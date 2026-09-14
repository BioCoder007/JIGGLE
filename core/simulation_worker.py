"""
Runs an OpenMM simulation on a background QThread and emits a signal
every `report_interval` steps with the current coordinates and a handful
of live metrics (RMSD, Rg, interface distance, energies).

Kept intentionally simple:
- Single prepared complex PDB in (receptor + binder already docked/posed),
  or a single protein/monomer (leave binder_chains empty).
- PDBFixer used to add missing atoms/hydrogens, and to strip
  crystallization additives (see strip_heterogens below).
- Implicit solvent (GBSA-OBC2) by default for speed, since the point is
  real-time visualization, not production-grade sampling. Explicit
  solvent is available as an option but will be much slower per frame.
- Receptor/binder split is done by chain ID for the interface-distance
  metric and for the live 3D view's chain colors.
- All live metrics (RMSD, Rg) and everything sent to the 3D viewer are
  restricted to the protein atoms (receptor + binder chains) — solvent is
  fully simulated but never rendered or included in the metrics. See the
  notes above `protein_idx` below for why this matters.
"""
import io
import os
import threading
import time
from collections import Counter

import numpy as np
from openmm import (
    CustomExternalForce,
    LangevinMiddleIntegrator,
    OpenMMException,
    Platform,
    app,
    unit,
)
from pdbfixer import PDBFixer
from PySide6.QtCore import QThread, Signal

from .analysis import (
    RunningRMSF,
    interface_distance,
    radius_of_gyration,
    rmsd,
    rmsd_fit_on,
)
from .forcefields import DEFAULT_FORCE_FIELD, DEFAULT_IMPLICIT, resolve


class SimulationWorker(QThread):
    frame_ready = Signal(dict)   # see payload keys below
    structure_ready = Signal(str)  # fixed-up PDB text, emitted once before stepping
    log = Signal(str)
    finished_ok = Signal()
    failed = Signal(str)

    def __init__(self, pdb_path, receptor_chains, binder_chains,
                 temperature_k=300.0, total_steps=1_000_000,
                 report_interval=50, timestep_fs=2.0,
                 solvent="implicit", output_dir="./md_output",
                 platform_name="auto", strip_heterogens=True,
                 force_field=DEFAULT_FORCE_FIELD, implicit_model=DEFAULT_IMPLICIT,
                 water_model=None, minimization_steps=0):
        super().__init__()
        self.pdb_path = pdb_path
        self.receptor_chains = set(receptor_chains)
        self.binder_chains = set(binder_chains)
        self.temperature_k = temperature_k
        self.total_steps = total_steps
        self.report_interval = report_interval
        self.timestep_fs = timestep_fs
        self.solvent = solvent
        self.output_dir = output_dir
        self.platform_name = platform_name
        self.strip_heterogens = strip_heterogens
        self.force_field = force_field
        self.implicit_model = implicit_model
        self.water_model = water_model
        # 0 means "run to OpenMM's own convergence tolerance", which is what
        # maxIterations=0 means to minimizeEnergy(). Any positive value is a
        # hard cap on iterations.
        self.minimization_steps = minimization_steps
        self._stop_requested = False
        # Set = running. Cleared = paused. Starts set so a run that is never
        # paused behaves exactly as before.
        self._resume = threading.Event()
        self._resume.set()

    def request_stop(self):
        self._stop_requested = True
        # Never leave the worker blocked on a pause it can no longer be
        # released from — Stop has to win over Pause.
        self._resume.set()

    def request_pause(self):
        self._resume.clear()

    def request_resume(self):
        self._resume.set()

    @property
    def is_paused(self):
        return not self._resume.is_set()

    def _pick_platform(self):
        if self.platform_name != "auto":
            try:
                return Platform.getPlatformByName(self.platform_name)
            except OpenMMException:
                self.log.emit(
                    f'Platform "{self.platform_name}" is not available on '
                    f"this OpenMM install — falling back to auto-detect."
                )
        for name in ("CUDA", "OpenCL", "CPU"):
            try:
                return Platform.getPlatformByName(name)
            except OpenMMException:
                # An unavailable platform is the normal case here, not an
                # error — this loop exists precisely to probe for one that
                # works, so there is nothing worth logging per miss.
                continue
        return Platform.getPlatformByName("Reference")

    def _chain_atom_indices(self, topology):
        """Recomputed AFTER the topology is finalized (i.e. after solvent
        has been added, if explicit) — doing this on the pre-solvation
        topology and reusing the indices post-solvation was a real bug:
        Modeller.addSolvent() builds a new Topology object, and while it
        happens to preserve solute atom order in practice, relying on
        stale indices from a different Topology object is fragile and
        was never actually verified correct."""
        receptor_idx, binder_idx = [], []
        for atom in topology.atoms():
            chain_id = atom.residue.chain.id
            if chain_id in self.receptor_chains:
                receptor_idx.append(atom.index)
            elif chain_id in self.binder_chains:
                binder_idx.append(atom.index)
        return np.array(receptor_idx, dtype=int), np.array(binder_idx, dtype=int)

    def run(self):
        try:
            self.log.emit("Fixing structure (PDBFixer)...")
            fixer = PDBFixer(filename=self.pdb_path)
            fixer.findMissingResidues()

            # Strip crystallization/purification additives (isopropanol,
            # glycerol, PEG, sulfate, DMSO, acetate, etc.) BEFORE adding
            # missing atoms/hydrogens. These show up constantly in real
            # PDB files (crystal structures especially) and have no force
            # field template in amber14-all.xml, which otherwise fails
            # with a cryptic "No template found for residue N (XXX)"
            # error deep inside createSystem(). Order matters: this must
            # run before findMissingAtoms/addMissingAtoms, or PDBFixer may
            # try to complete/model atoms onto residues we're about to
            # throw away anyway.
            if self.strip_heterogens:
                # keepWater depends on solvent mode: implicit solvent
                # (GBSA-OBC2) has NO water residue template at all — it
                # represents solvent effects implicitly, without any
                # explicit water molecules — so keeping crystallographic
                # waters here crashes createSystem() with "No template
                # found for residue N (HOH)" on ANY structure that has
                # ordered waters (extremely common in real crystal
                # structures). Found this via an actual end-to-end test
                # run, not by inspection. Explicit solvent gets its own
                # fresh water box from addSolvent() anyway, so crystal
                # waters aren't needed there either — stripping them
                # avoids potential clashes with the newly added box.
                keep_water = False
                removed = fixer.removeHeterogens(keepWater=keep_water)
                if removed:
                    counts = Counter(r.name for r in removed)
                    summary = ", ".join(f"{name} x{n}" for name, n in counts.items())
                    self.log.emit(
                        f"Stripped {len(removed)} non-standard residue(s) "
                        f"(crystallization additives / heterogens, "
                        f"including crystallographic water): {summary}"
                    )

            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(7.0)

            topology = fixer.topology
            positions = fixer.positions

            xml_files, water_model_name = resolve(
                self.force_field, self.solvent,
                implicit_name=self.implicit_model, water_name=self.water_model,
            )
            self.log.emit(
                f"Building system — force field: {self.force_field} "
                f"[{', '.join(xml_files)}], {self.solvent} solvent"
                + (f", water model '{water_model_name}'" if water_model_name else "")
            )
            try:
                forcefield = app.ForceField(*xml_files)
                if self.solvent == "implicit":
                    system = forcefield.createSystem(
                        topology, nonbondedMethod=app.NoCutoff,
                        constraints=app.HBonds,
                    )
                else:
                    modeller = app.Modeller(topology, positions)
                    # Passing model= is REQUIRED for anything that isn't a
                    # 3-site water. addSolvent defaults to 'tip3p', so a 4-site
                    # force field would be handed 3-site waters and fail deep
                    # inside createSystem with "No template found for residue
                    # N (HOH)" — a real trap, confirmed by direct test.
                    modeller.addSolvent(
                        forcefield, model=water_model_name,
                        padding=1.0 * unit.nanometer,
                    )
                    topology, positions = modeller.topology, modeller.positions
                    # addSolvent centres the water on the solute, but it does
                    # NOT move anything to the box's own coordinate frame — the
                    # whole system keeps whatever coordinates the input PDB had.
                    # That is fine for the physics (PME is periodic either way),
                    # but every file written with enforcePeriodicBox=True then
                    # wraps molecules into [0, L]. Water is thousands of small
                    # molecules so it tiles the box perfectly; the protein is
                    # ONE molecule, so it gets translated as a unit to wherever
                    # its centre happens to land — typically a corner, half
                    # outside the water. That is what makes topology.pdb look
                    # broken in PyMOL even though the simulation is correct.
                    #
                    # Translating the system so the solute sits at the box
                    # centre costs nothing and makes the written files look
                    # like what they physically are.
                    box_nm = np.array(
                        topology.getUnitCellDimensions().value_in_unit(unit.nanometer)
                    )
                    pos_nm = np.array(positions.value_in_unit(unit.nanometer))
                    solute_mask = np.array([
                        a.residue.name not in ("HOH", "WAT")
                        for a in topology.atoms()
                    ])
                    solute = pos_nm[solute_mask]
                    solute_centre = 0.5 * (solute.min(axis=0) + solute.max(axis=0))
                    shift = box_nm / 2.0 - solute_centre
                    positions = (pos_nm + shift) * unit.nanometer
                    self.log.emit(
                        f"Centred the solute in the {box_nm[0]:.1f} x {box_nm[1]:.1f} x "
                        f"{box_nm[2]:.1f} nm water box (shifted by "
                        f"{shift[0]:.1f}, {shift[1]:.1f}, {shift[2]:.1f} nm) so the "
                        f"written topology.pdb and trajectory are correctly imaged."
                    )
                    system = forcefield.createSystem(
                        topology, nonbondedMethod=app.PME,
                        nonbondedCutoff=1.0 * unit.nanometer,
                        constraints=app.HBonds,
                    )
            except ValueError as exc:
                if "No template found" in str(exc):
                    raise RuntimeError(
                        f"{exc}\n\n"
                        "This residue has no amber14 force field parameters. "
                        "If it's a crystallization additive (glycerol, PEG, "
                        "sulfate, DMSO, etc.) it should have been stripped "
                        "automatically — check the log above for what got "
                        "removed. If it's a real ligand/cofactor you want "
                        "kept (e.g. a small-molecule inhibitor, a metal not "
                        "covered by amber14), it needs its own force field "
                        "parameters (e.g. GAFF/OpenFF via openmmforcefields) "
                        "before OpenMM can build a system with it — this "
                        "app doesn't generate those automatically."
                    ) from exc
                raise

            # Recompute chain groups on the FINAL topology (post-solvation
            # if explicit) — see _chain_atom_indices docstring.
            receptor_idx, binder_idx = self._chain_atom_indices(topology)
            protein_idx = np.sort(np.concatenate([receptor_idx, binder_idx])) \
                if len(receptor_idx) or len(binder_idx) else None
            if protein_idx is None or len(protein_idx) == 0:
                raise RuntimeError(
                    "No atoms matched the receptor/binder chain IDs given — "
                    "check the 'Detected chains' line against what you typed "
                    "into the chain fields."
                )
            self.log.emit(
                f"Tracking {len(protein_idx)} protein atom(s) for live view "
                f"and metrics (receptor: {len(receptor_idx)}, binder: {len(binder_idx)})"
                + (" — solvent excluded from both." if self.solvent == "explicit" else "")
            )

            # --- Positional restraints during minimization/heating -----------
            # Standard MD equilibration practice, and a real gap in the round-5
            # heating fix: without this, ramping temperature on a freshly-
            # minimized structure can genuinely shock it out of shape — heavy
            # atoms are free to swing however far a few hundred K of thermal
            # energy pushes them before the system has settled. Tethering
            # heavy (non-hydrogen) atoms to their starting positions with a
            # moderate harmonic spring during minimization+heating prevents
            # this, then the restraint is fully released (k -> 0) before the
            # production run, so the reported dynamics are completely
            # unrestrained. This is a standard technique in essentially every
            # published MD equilibration protocol (GROMACS/AMBER tutorials,
            # etc.) — its absence here was a real bug, not a rendering issue.
            # Two things here were wrong and together produced
            # "Particle coordinate is NaN" on every explicit-solvent run:
            #
            # 1. The loop restrained EVERY heavy atom in the topology, which
            #    after solvation means every water oxygen too — 4,268 atoms on
            #    a system whose protein has 121. Tethering the solvent is not
            #    what "prevent thermal shock distorting the structure" means,
            #    and it makes the water a rigid lattice instead of a liquid.
            #
            # 2. "(x-x0)^2+..." is a NON-PERIODIC distance. Under PME, OpenMM
            #    wraps molecules back into the primary box as they drift out.
            #    The instant a restrained water wraps, x jumps by a full box
            #    vector and the restraint sees a ~5 nm displacement. Measured
            #    on a 5.1 nm box: 13,020 kJ/mol of spurious energy and a
            #    5,103 kJ/mol/nm force on a single atom. At 2 fs that is an
            #    instant blow-up. With thousands of restrained waters over
            #    1,500 heating steps it is not a question of whether one wraps.
            #
            # Fix: restrain protein heavy atoms only, and measure the
            # displacement with periodicdistance() so a wrap costs nothing.
            # periodicdistance needs real box vectors, so it is only used when
            # there is a periodic box — in implicit solvent the system's
            # default box is meaningless and the plain form is correct.
            if self.solvent == "explicit":
                restraint_expr = "0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2"
            else:
                restraint_expr = "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
            restraint = CustomExternalForce(restraint_expr)
            restraint.addGlobalParameter("k", 1000.0 * unit.kilojoule_per_mole / unit.nanometer**2)
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addPerParticleParameter("z0")
            n_restrained = 0
            atoms_by_index = {a.index: a for a in topology.atoms()}
            for i in protein_idx:
                atom = atoms_by_index[int(i)]
                if atom.element is None or atom.element.symbol == "H":
                    continue
                pos = positions[atom.index].value_in_unit(unit.nanometer)
                restraint.addParticle(atom.index, list(pos))
                n_restrained += 1
            system.addForce(restraint)
            self.log.emit(
                f"Restraining {n_restrained} protein heavy atom(s) to their starting "
                f"positions during minimization/heating (released before the main run) "
                f"to prevent thermal shock distorting the structure. Solvent is left "
                f"free"
                + (" and displacements are measured with the periodic minimum image."
                   if self.solvent == "explicit" else ".")
            )

            integrator = LangevinMiddleIntegrator(
                self.temperature_k * unit.kelvin,
                1.0 / unit.picosecond,
                self.timestep_fs * unit.femtosecond,
            )
            platform = self._pick_platform()
            simulation = app.Simulation(topology, system, integrator, platform)
            simulation.context.setPositions(positions)

            os.makedirs(self.output_dir, exist_ok=True)
            dcd_path = os.path.join(self.output_dir, "trajectory.dcd")
            top_path = os.path.join(self.output_dir, "topology.pdb")
            # --- Minimization -------------------------------------------------
            # Generated/docked/AI-designed complexes routinely have local
            # clashes and imperfect bond geometry that a shallow, hard-capped
            # minimization won't resolve. maxIterations=0 tells OpenMM to run
            # until convergence (its default tolerance) rather than stopping
            # at an arbitrary step count — this was previously capped at 500,
            # which was nowhere near enough for a rough input structure and
            # was a likely contributor to the instability seen below.
            if self.minimization_steps == 0:
                self.log.emit("Minimizing energy (running to convergence)...")
            else:
                self.log.emit(
                    f"Minimizing energy (hard cap: {self.minimization_steps} iterations — "
                    f"set 0 to run to convergence instead)..."
                )
            simulation.minimizeEnergy(maxIterations=self.minimization_steps)
            pe_after_min = simulation.context.getState(getEnergy=True) \
                .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            self.log.emit(f"Post-minimization potential energy: {pe_after_min:.1f} kJ/mol")

            # --- Gradual heating ------------------------------------------------
            # Assigning full-temperature velocities in one shot onto a
            # freshly-minimized (but not necessarily fully relaxed)
            # structure is a classic way to trigger an instant blow-up —
            # a large local strain released as kinetic energy all at once
            # can be enough to send the integrator unstable within a few
            # steps, after which every subsequent frame is garbage (this
            # matches "structure is distorted / plots are flat" reports).
            # Ramping from a low temperature over a short equilibration
            # phase gives the system a chance to dissipate that strain
            # gradually instead. This equilibration time is IN ADDITION to
            # total_steps — it's prep, not part of the requested run.
            n_heat_stages = 6
            heat_steps_per_stage = 250
            self.log.emit(
                f"Heating gradually to {self.temperature_k:.0f} K "
                f"({n_heat_stages} stages, not counted toward Total steps)..."
            )
            for stage in range(1, n_heat_stages + 1):
                if self._stop_requested:
                    self.log.emit("Stop requested during heating — halting.")
                    self.finished_ok.emit()
                    return
                # Pause is honoured during equilibration too. Without this a
                # user who hits Pause while the system is still heating sees
                # nothing happen until production begins, then watches it stop
                # immediately — which reads as a bug.
                if not self._resume.is_set():
                    self.log.emit("Paused (during heating).")
                    self._resume.wait()
                    if self._stop_requested:
                        self.log.emit("Stop requested during heating — halting.")
                        self.finished_ok.emit()
                        return
                    self.log.emit("Resumed.")
                stage_temp = self.temperature_k * stage / n_heat_stages
                integrator.setTemperature(stage_temp * unit.kelvin)
                simulation.context.setVelocitiesToTemperature(stage_temp * unit.kelvin)
                try:
                    simulation.step(heat_steps_per_stage)
                except Exception as exc:
                    if "NaN" in str(exc):
                        raise RuntimeError(
                            f"The integrator became unstable during heating "
                            f"(stage {stage}/{n_heat_stages}, {stage_temp:.0f} K): {exc}\n\n"
                            f"Usual causes, in order of likelihood:\n"
                            f"  - severe steric clashes in the input that minimization "
                            f"could not resolve (check the structure, or minimize it "
                            f"externally first),\n"
                            f"  - a timestep too large for the system — try 1.0 fs,\n"
                            f"  - an inappropriate force field / water pairing."
                        ) from exc
                    raise

                state_check = simulation.context.getState(getPositions=True)
                check_coords = state_check.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                if not np.all(np.isfinite(check_coords)):
                    raise RuntimeError(
                        f"Simulation became numerically unstable during heating "
                        f"(stage {stage}/{n_heat_stages}, {stage_temp:.0f} K) — "
                        f"positions turned into NaN/Inf. This usually means the "
                        f"input structure has severe clashes minimization "
                        f"couldn't resolve, or the timestep is too large. Try a "
                        f"smaller timestep, or inspect the input structure for "
                        f"steric clashes before re-running."
                    )

            state_heated = simulation.context.getState(getEnergy=True)
            pe_after_heat = state_heated.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            ke_after_heat = state_heated.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            self.log.emit(
                f"Post-heating: PE {pe_after_heat:.1f} kJ/mol, "
                f"KE {ke_after_heat:.1f} kJ/mol"
            )

            # Release the positional restraint now that heating is done —
            # everything from here on (the RMSD reference capture and the
            # entire production run) is fully unrestrained real dynamics.
            simulation.context.setParameter("k", 0.0)
            self.log.emit("Released positional restraints — starting unrestrained dynamics.")

            # The DCD reporter is registered HERE, not before minimization.
            # Registering it earlier (as this previously did) meant the six
            # 250-step heating stages were written to the trajectory as ~30
            # restrained frames, so the DCD and the live plots were offset
            # from each other and anyone opening the file in VMD/ChimeraX saw
            # half a second of tethered, non-physical motion before the real
            # run began. The trajectory now contains production frames only,
            # and frame N in the file is frame N on the plots.
            simulation.reporters.append(
                app.DCDReporter(dcd_path, self.report_interval)
            )
            self.log.emit(
                f"Writing production trajectory to {dcd_path} "
                f"(equilibration frames excluded; companion topology: {top_path})"
            )

            # RMSD reference is captured AFTER heating (not after raw
            # minimization) — a more physically meaningful "start" state,
            # and avoids measuring RMSD relative to a pre-equilibration
            # snapshot that the system immediately moves away from anyway.
            state0 = simulation.context.getState(
                getPositions=True, enforcePeriodicBox=(self.solvent == "explicit")
            )
            all_coords0 = state0.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            ref_coords = all_coords0[protein_idx]

            # Build a PROTEIN-ONLY topology/positions for the live 3D view
            # and the emitted reference structure. Sending solvent atoms to
            # the viewer was the direct cause of "explicit solvent renders
            # very badly" — thousands of extra water lines/points on top of
            # the actual structure of interest. Solvent is still fully
            # simulated under the hood; it's just never rendered or
            # included in RMSD/Rg.
            view_modeller = app.Modeller(topology, state0.getPositions())
            protein_idx_set = {int(i) for i in protein_idx}
            atoms_to_remove = [a for a in topology.atoms() if a.index not in protein_idx_set]
            view_modeller.delete(atoms_to_remove)
            view_topology, view_positions = view_modeller.topology, view_modeller.positions

            with open(top_path, "w") as f:
                # keepIds=True is essential — PDBFile.writeFile defaults to
                # keepIds=False, which silently renumbers residues AND
                # relabels chains sequentially (A, B, C, ...) regardless of
                # the topology's real chain IDs. Without this, the DCD's
                # companion topology.pdb would carry different chain labels
                # than what the user actually specified/sees in the app.
                app.PDBFile.writeFile(topology, state0.getPositions(), f, keepIds=True)  # full system, incl. solvent, for post-hoc analysis

            pdb_buffer = io.StringIO()
            # keepIds=True here is the critical fix for a real bug: without
            # it, the chain IDs sent to the 3D viewer get silently
            # relabeled to sequential A/B/C — breaking the receptor/binder
            # color matching (any chain the user typed no longer matches
            # anything, since the labels the viewer sees aren't the labels
            # the user typed). This was caught from an actual run's log:
            # "[3D viewer] structure loaded OK. chains=A,B,C;unmatched=B,C"
            # when the real chains were H, L, A.
            app.PDBFile.writeFile(view_topology, view_positions, pdb_buffer, keepIds=True)

            # OpenMM's PDB writer NEVER emits HELIX/SHEET header records —
            # confirmed directly (write a trivial topology, grep the output:
            # neither keyword appears). 3Dmol explicitly parses HELIX/SHEET
            # when present and trusts them as ground truth; when absent, it
            # falls back to its own geometric hydrogen-bond-distance-based
            # secondary-structure guesser, which is far less reliable and
            # was the actual cause of the "tangled/distorted" cartoon
            # reports — the physics was fine (RMSD stayed small and sane),
            # but the cartoon was being drawn from a bad guess at which
            # residues are helix/sheet/coil, not from real per-atom error.
            # This affected every structure sent through this pipeline
            # (never just one file), while "Preview Structure" — which
            # loads the ORIGINAL uploaded file directly — was always fine,
            # because real deposited PDB files usually carry real
            # author-curated HELIX/SHEET records that 3Dmol reads directly.
            # Fix: carry the original file's HELIX/SHEET lines forward into
            # what we send to the viewer, since chain IDs and residue
            # numbers are preserved (keepIds=True, and PDBFixer doesn't
            # renumber existing residues). This only helps files that HAD
            # such records to begin with — AI-generated/docked outputs
            # that never had them will still fall back to the geometric
            # guess, same as they always would.
            helix_sheet_lines = []
            try:
                with open(self.pdb_path, errors="ignore") as f:
                    for line in f:
                        if line.startswith(("HELIX", "SHEET")):
                            helix_sheet_lines.append(line.rstrip("\n"))
            except OSError:
                pass

            view_pdb_text = pdb_buffer.getvalue()
            if helix_sheet_lines:
                lines = view_pdb_text.split("\n")
                insert_at = 1 if lines and lines[0].startswith("REMARK") else 0
                lines[insert_at:insert_at] = helix_sheet_lines
                view_pdb_text = "\n".join(lines)
                self.log.emit(
                    f"Carried forward {len(helix_sheet_lines)} HELIX/SHEET "
                    f"record(s) from the original file so the 3D view uses "
                    f"real secondary structure instead of a geometric guess."
                )
            else:
                self.log.emit(
                    "No HELIX/SHEET records in the original file — the 3D "
                    "view's cartoon shape is a geometric guess (3Dmol's own "
                    "hydrogen-bond-based detector), not ground truth. This "
                    "is normal for generated/docked structures that never "
                    "had this annotation."
                )

            self.structure_ready.emit(view_pdb_text)

            all_masses = np.array([
                system.getParticleMass(i).value_in_unit(unit.dalton)
                for i in range(system.getNumParticles())
            ])
            protein_masses = all_masses[protein_idx]

            # Per-residue RMSF needs the CA atoms' positions WITHIN the
            # protein-only coordinate slice, not their indices in the full
            # system, plus a label per residue for the x axis.
            protein_pos_of = {int(g): i for i, g in enumerate(protein_idx)}
            ca_local_idx, rmsf_labels, rmsf_resids = [], [], []
            for atom in topology.atoms():
                if atom.name != "CA" or atom.index not in protein_pos_of:
                    continue
                ca_local_idx.append(protein_pos_of[atom.index])
                res = atom.residue
                # The x axis is a CONTIGUOUS 1..N index, deliberately not the
                # PDB residue number. Residue numbering restarts at 1 for every
                # chain, so in any receptor/binder complex the raw ids run
                # 1..24, 1..12 — plotting against that sends the line jumping
                # back to the left half way across, drawing a criss-cross of
                # straight segments over the real trace. The true chain and
                # residue number are preserved in the label and written to
                # rmsf.csv, so nothing is lost.
                rmsf_resids.append(len(rmsf_resids) + 1)
                rmsf_labels.append(f"{res.chain.id}:{res.name}{res.id}")
            ca_local_idx = np.array(ca_local_idx, dtype=int)
            # x positions where the chain changes, for a divider line on the plot
            chain_breaks = [i + 1 for i in range(1, len(rmsf_labels))
                            if rmsf_labels[i].split(":")[0] != rmsf_labels[i - 1].split(":")[0]]

            # Receptor/binder positions WITHIN the protein-only slice, so the
            # split RMSDs can be computed on the same array the view uses.
            rec_local = np.array([protein_pos_of[int(i)] for i in receptor_idx
                                  if int(i) in protein_pos_of], dtype=int)
            bnd_local = np.array([protein_pos_of[int(i)] for i in binder_idx
                                  if int(i) in protein_pos_of], dtype=int)
            split_rmsd = len(rec_local) >= 3 and len(bnd_local) >= 3
            if split_rmsd:
                self.log.emit(
                    "Two-body system — reporting receptor RMSD, binder internal "
                    "RMSD, and binder RMSD fitted on the receptor (pose drift) "
                    "alongside the whole-complex value."
                )

            rmsf_tracker = None
            if len(ca_local_idx) >= 3:
                rmsf_tracker = RunningRMSF(ref_coords[ca_local_idx])
                self.log.emit(
                    f"Tracking per-residue RMSF over {len(ca_local_idx)} CA atoms "
                    f"(accumulated incrementally, so it refines as the run proceeds)."
                )
            else:
                self.log.emit(
                    "Fewer than 3 CA atoms found — skipping the RMSF plot."
                )

            n_chunks = self.total_steps // self.report_interval
            self.log.emit(f"Running {self.total_steps} steps "
                           f"({n_chunks} reporting frames)...")

            t_start = time.time()
            for chunk in range(n_chunks):
                if self._stop_requested:
                    self.log.emit("Stop requested — halting simulation.")
                    break

                # Pause: block here rather than spinning. The integrator is
                # simply not advanced, so the physics is untouched — this is a
                # true pause, not a slow-down, and resuming continues from the
                # exact same state. Woken by request_resume() or request_stop().
                if not self._resume.is_set():
                    self.log.emit("Paused.")
                    self._resume.wait()
                    if self._stop_requested:
                        self.log.emit("Stop requested — halting simulation.")
                        break
                    self.log.emit("Resumed.")

                simulation.step(self.report_interval)

                state = simulation.context.getState(
                    getPositions=True, getEnergy=True,
                    enforcePeriodicBox=(self.solvent == "explicit"),
                )
                all_coords = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)

                if not np.all(np.isfinite(all_coords)):
                    raise RuntimeError(
                        f"Simulation became numerically unstable at step "
                        f"{(chunk + 1) * self.report_interval} — positions "
                        f"turned into NaN/Inf. Stopping rather than render "
                        f"garbage frames. Try a smaller timestep or check "
                        f"the input structure for clashes."
                    )

                coords = all_coords[protein_idx]  # protein-only, for view + metrics
                pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)

                current_rmsd = rmsd(coords, ref_coords, align=True)
                current_rg = radius_of_gyration(coords, protein_masses)
                iface_dist = None
                if len(receptor_idx) and len(binder_idx):
                    # receptor_idx/binder_idx are indices into the FULL
                    # (solvent-included) coordinate array, not the
                    # protein-only slice, so use all_coords here.
                    iface_dist = interface_distance(all_coords, receptor_idx, binder_idx)

                step_no = (chunk + 1) * self.report_interval
                payload = {
                    "step": step_no,
                    "time_ps": step_no * self.timestep_fs / 1000.0,
                    "positions_nm": coords,           # protein-only Nx3 numpy array
                    "rmsd_nm": current_rmsd,
                    "rg_nm": current_rg,
                    "interface_distance_nm": iface_dist,
                    "potential_energy_kjmol": pe,
                    "kinetic_energy_kjmol": ke,
                }

                if split_rmsd:
                    payload["rmsd_receptor_nm"] = rmsd(
                        coords[rec_local], ref_coords[rec_local], align=True)
                    payload["rmsd_binder_self_nm"] = rmsd(
                        coords[bnd_local], ref_coords[bnd_local], align=True)
                    # The pose-stability number: hold the receptor still, then
                    # measure how far the binder has moved in the site.
                    payload["rmsd_binder_on_receptor_nm"] = rmsd_fit_on(
                        coords, ref_coords, rec_local, bnd_local)

                if rmsf_tracker is not None:
                    rmsf_tracker.add_frame(coords[ca_local_idx])
                    payload["rmsf_nm"] = rmsf_tracker.values()
                    payload["rmsf_resids"] = rmsf_resids
                    payload["rmsf_labels"] = rmsf_labels
                    payload["rmsf_chain_breaks"] = chain_breaks
                    payload["rmsf_n_frames"] = rmsf_tracker.n_frames
                self.frame_ready.emit(payload)

            elapsed = time.time() - t_start
            self.log.emit(f"Done. {elapsed:.1f}s wall time.")
            self.finished_ok.emit()

        except Exception as exc:
            # Deliberately broad. This runs on a worker thread, so anything
            # not caught here dies silently and the UI hangs on "running"
            # forever. Every failure must reach the user.
            self.failed.emit(str(exc))
