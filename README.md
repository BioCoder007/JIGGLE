# Jiggle

*"Everything that living things do can be understood in terms of the jigglings and wigglings of atoms."* — Richard Feynman

A desktop app that runs a molecular dynamics simulation and renders the structure moving **while it runs** — not playback of a finished trajectory.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Lint](https://github.com/BioCoder007/JIGGLE/actions/workflows/lint.yml/badge.svg)](https://github.com/BioCoder007/JIGGLE/actions/workflows/lint.yml)

Point it at a PDB, choose a force field, hit Start. PDBFixer repairs the
structure, OpenMM minimises it, heats it gradually under positional
restraints, then releases into production dynamics — and you watch the whole
thing in a 3D viewer with live plots updating beside it.

---

## Install

```bash
conda create -n jiggle python=3.11 -y
conda activate jiggle
conda install -c conda-forge openmm pdbfixer -y
pip install -r requirements.txt
```

Install OpenMM through conda-forge rather than pip — the pip wheels are
unreliable for CUDA. Check what platforms you actually have:

```bash
python -m openmm.testInstallation
```

## Run

```bash
python main.py
```

---

## What you can configure

**Force fields** — AMBER14 (ff14SB), AMBER19 (ff19SB), AMBER99SB-ILDN,
AMBER99SB, AMBER10, AMBER03, CHARMM36.

**Solvent** — implicit (OBC2, OBC1, GBn2, GBn, HCT) or explicit
(TIP3P-FB, TIP3P, SPC/E, TIP4P-Ew, matched to the force field).

Every force-field/solvent pairing in the dropdowns was verified by building
an actual `System` from it, so the UI cannot offer a combination that fails
at runtime. Water models are declared with the name `Modeller.addSolvent`
needs, which matters: 4-site waters fail with a misleading *"No template
found for residue N (HOH)"* if the model name is not passed.

**Run length** — enter nanoseconds or steps; the two stay in sync as you
change the timestep. Minimisation iterations, temperature, timestep and
reporting interval are all exposed. Platform is CUDA / OpenCL / CPU / auto.

## Live metrics

| Plot | What it tells you |
|---|---|
| **RMSD** | Structural drift from the starting conformation |
| **Radius of gyration** | Compactness — expansion or collapse |
| **Potential energy** | Whether the system is behaving |
| **Per-residue RMSF** | Which regions are flexible and which are rigid |
| **Interface distance** | Receptor–binder centroid separation |

### RMSD is split four ways for two-body systems

When a binder chain is set, the RMSD panel shows:

- **complex** — everything, aligned on everything
- **receptor** — is the receptor itself stable?
- **binder (internal)** — is the binder changing shape?
- **binder (fit on receptor)** — superpose on the receptor, then measure the
  binder. **This is the pose-stability number.**

That last one matters because a single combined RMSD averages a binder's
drift across every receptor atom that did not move. In a test where a binder
slid exactly 0.500 nm while the receptor stayed put, the combined trace
reported 0.120 nm; fitting on the receptor recovered the full 0.500 nm. A
centroid interface distance misses it too, since it cannot see rotation or
sliding that preserves centroid separation.

Per-residue RMSF is accumulated incrementally from running sums rather than
stored coordinates, so it refines live and uses memory proportional to the
number of residues, not the length of the run.

## Other features

- **Pause / resume** mid-run. The integrator is not advanced while paused, so
  this is a true pause and resuming continues from the identical state.
- **Save any live frame** straight to PDB, with step, time, RMSD and Rg in the
  header.
- **Automatic chain detection** with a receptor/binder guess you can override.
- **Preview Structure** loads the raw file instantly, before any processing.
- **Monomer support** — leave the binder field empty.
- Receptor/binder colour pickers that re-render immediately.

## Outputs

Written to the output directory on completion, or on demand via
*Save plots + data*:

```
metrics.csv        every series, per frame (including the split RMSDs)
rmsf.csv           per-residue RMSF with chain and residue labels
rmsd.png  radius_of_gyration.png  potential_energy.png  rmsf.png  interface_distance.png
trajectory.dcd     production frames only
topology.pdb       matching topology, solute centred in the box
```

`trajectory.dcd` and `topology.pdb` load together in VMD, PyMOL, ChimeraX,
MDAnalysis or mdtraj, and feed straight into MM/GBSA tooling such as
`gmx_MMPBSA`.

Two details worth knowing:

- The DCD contains **production frames only**. The reporter is registered
  after the restraints are released, so equilibration is not written and
  frame *N* in the file is frame *N* on the plots.
- For explicit solvent, the solute is translated to the box centre before the
  run. `addSolvent` centres water on the solute but leaves everything in the
  input file's coordinate frame; periodic imaging then wraps the protein — one
  molecule — to a corner while the water tiles the box, which makes a
  perfectly valid simulation look broken in PyMOL.

---

## How it works

`core/simulation_worker.py` runs OpenMM on a `QThread`. Every
`report_interval` steps it pulls positions and energies from the context and
emits a `frame_ready` signal. The main thread converts nm to Angstrom and
pushes a flat coordinate array into the embedded 3Dmol.js viewer, and appends
to the pyqtgraph plots.

The first structure sent to the viewer is the PDBFixer-processed,
post-minimisation topology — not your input file. PDBFixer adds hydrogens and
missing atoms, so the viewer's atom list has to match what OpenMM is
integrating or every later coordinate push silently misaligns. The viewer
parses with `keepH: true`; without it 3Dmol drops every hydrogen at parse
time and the atom counts diverge.

Only protein atoms are sent to the viewer and used for metrics. Solvent is
fully simulated but excluded from both — a live cartoon of 12,000 water
molecules is neither useful nor fast, and Rg over a water box measures the
box.

## Limitations

- **Interface distance is a centroid separation**, not a reaction coordinate.
  A real free-energy profile along dissociation needs umbrella sampling or
  metadynamics — a separate, much longer job.
- **No binding free energy.** A complex trajectory tells you about pose
  stability, not affinity. For binding free energy use MM/GBSA or FEP on the
  output.
- **NVT only.** No barostat, so no NPT density equilibration.
- **Explicit solvent is slow to watch.** Lower `Steps / frame` and expect a
  slideshow unless you are on a strong GPU.
- **Standard residues only.** Non-canonical amino acids, cyclised or stapled
  peptides, glycans and most cofactors need parameters this app does not
  generate — see `openmmforcefields` for those.
- **Secondary structure for generated structures.** HELIX/SHEET records are
  carried forward from the input file when present. RFdiffusion and docking
  outputs usually have none, so 3Dmol falls back to its geometric guess and
  cartoons look looser.
- No checkpoint/restart, no run queue, no trajectory scrubbing yet.

## File layout

```
main.py
requirements.txt
pyproject.toml            ruff configuration
core/
  analysis.py             Kabsch alignment, RMSD, fit-on RMSD, Rg, running RMSF
  forcefields.py          verified force field / solvent catalogue
  simulation_worker.py    OpenMM on a QThread, streams frames
ui/
  main_window.py          PySide6 layout, plots, exports
  viewer_template.html    3Dmol.js live viewer
  vendor/3Dmol-min.js     vendored, BSD licensed - no CDN dependency
```

## Contributing

Issues and PRs welcome — see [CONTRIBUTING.md](CONTRIBUTING.md). Run
`ruff check .` before opening a PR; CI runs the same check with a pinned ruff
version.

## License

MIT — see [LICENSE](LICENSE).
