# Contributing

Thanks for considering a contribution. This is an early-stage tool — the
core loop (OpenMM worker → live 3D view → live plots) works, but there's
plenty of room to extend it.

## Setup

```bash
git clone https://github.com/<your-org>/md-realtime-gui.git
cd md-realtime-gui
conda create -n mdgui python=3.11 -y
conda activate mdgui
conda install -c conda-forge openmm pdbfixer -y
pip install -r requirements.txt
```

## Good first issues

- DCD trajectory writing alongside the live view (see README "Known
  limitations")
- Pause / resume controls
- Hiding solvent atoms in the 3D view when running explicit solvent
- Swapping the crude interface-distance metric for a proper collective
  variable
- Multi-run queue (batch several complexes without relaunching the app)

## Pull requests

- Keep PRs scoped to one change.
- If you touch `core/simulation_worker.py`, note any performance impact
  on live reporting (it runs every N steps, so overhead compounds).
- If you touch `ui/viewer_template.html`, test with both a small
  (~2-3k atom, implicit solvent) and larger (explicit solvent) system —
  3Dmol coordinate updates scale with atom count.
- Add yourself to a CONTRIBUTORS section in the README if you'd like.

## Reporting bugs

Use the bug report issue template. Include your OS, GPU (if any),
OpenMM platform (`CUDA`/`OpenCL`/`CPU`), and the PDB you tested with if
you can share it.
