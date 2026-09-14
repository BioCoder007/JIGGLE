# Citing Jiggle

Jiggle is a front end. It runs **OpenMM**, prepares structures with
**PDBFixer**, renders with **3Dmol.js**, and applies published force fields,
solvent models and algorithms that other people developed. Cite Jiggle for
what Jiggle did, and cite the underlying methods for the science.

The rule of thumb: **cite what your results depend on.** If a reader could
not reproduce your numbers without knowing you used it, it belongs in your
methods.

---

## 1. Jiggle

```bibtex
@software{safeer_jiggle_2026,
  author  = {Safeer, Amman},
  title   = {{Jiggle: a live molecular dynamics viewer}},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/BioCoder007/JIGGLE},
  license = {MIT}
}
```

Safeer, A. (2026). *Jiggle: a live molecular dynamics viewer* (Version 0.1.0)
[Computer software]. https://github.com/BioCoder007/JIGGLE

> Add the Zenodo DOI here once the first release is archived, and to the
> `doi:` field in `CITATION.cff`. A DOI is what makes this citable in a
> journal rather than a URL that can rot.

---

## 2. Always cite these

Every Jiggle run uses all three.

**OpenMM** — the simulation engine. Jiggle implements none of the physics.

> Eastman, P., Swails, J., Chodera, J. D., McGibbon, R. T., Zhao, Y.,
> Beauchamp, K. A., Wang, L.-P., Simmonett, A. C., Harrigan, M. P.,
> Stern, C. D., Wiewiora, R. P., Brooks, B. R., & Pande, V. S. (2017).
> OpenMM 7: Rapid development of high performance algorithms for molecular
> dynamics. *PLOS Computational Biology*, 13(7), e1005659.
> https://doi.org/10.1371/journal.pcbi.1005659

**PDBFixer** — structure repair: missing atoms and residues, hydrogens,
heterogen removal. Part of the OpenMM project; covered by the citation above.

**3Dmol.js** — the live 3D viewer. Cite if any figure came from Jiggle's
viewer.

> Rego, N., & Koes, D. (2015). 3Dmol.js: molecular visualization with WebGL.
> *Bioinformatics*, 31(8), 1322-1324.
> https://doi.org/10.1093/bioinformatics/btu829

---

## 3. Cite the force field you selected

**AMBER14 / ff14SB** (Jiggle default)

> Maier, J. A., Martinez, C., Kasavajhala, K., Wickstrom, L., Hauser, K. E.,
> & Simmerling, C. (2015). ff14SB: Improving the accuracy of protein side
> chain and backbone parameters from ff99SB. *Journal of Chemical Theory and
> Computation*, 11(8), 3696-3713. https://doi.org/10.1021/acs.jctc.5b00255

**AMBER19 / ff19SB**

> Tian, C., Kasavajhala, K., Belfon, K. A. A., Raguette, L., Huang, H.,
> Migues, A. N., Bickel, J., Wang, Y., Pincay, J., Wu, Q., & Simmerling, C.
> (2020). ff19SB: Amino-acid-specific protein backbone parameters trained
> against quantum mechanics energy surfaces in solution. *Journal of Chemical
> Theory and Computation*, 16(1), 528-552.
> https://doi.org/10.1021/acs.jctc.9b00591

**AMBER99SB-ILDN**

> Lindorff-Larsen, K., Piana, S., Palmo, K., Maragakis, P., Klepeis, J. L.,
> Dror, R. O., & Shaw, D. E. (2010). Improved side-chain torsion potentials
> for the Amber ff99SB protein force field. *Proteins*, 78(8), 1950-1958.
> https://doi.org/10.1002/prot.22711

**AMBER99SB**

> Hornak, V., Abel, R., Okur, A., Strockbine, B., Roitberg, A., & Simmerling,
> C. (2006). Comparison of multiple Amber force fields and development of
> improved protein backbone parameters. *Proteins*, 65(3), 712-725.
> https://doi.org/10.1002/prot.21123

**AMBER03**

> Duan, Y., Wu, C., Chowdhury, S., Lee, M. C., Xiong, G., Zhang, W., Yang, R.,
> Cieplak, P., Luo, R., Lee, T., Caldwell, J., Wang, J., & Kollman, P. (2003).
> A point-charge force field for molecular mechanics simulations of proteins.
> *Journal of Computational Chemistry*, 24(16), 1999-2012.
> https://doi.org/10.1002/jcc.10349

**CHARMM36**

> Best, R. B., Zhu, X., Shim, J., Lopes, P. E. M., Mittal, J., Feig, M., &
> MacKerell, A. D., Jr. (2012). Optimization of the additive CHARMM all-atom
> protein force field targeting improved sampling of the backbone phi, psi and
> side-chain chi1 and chi2 dihedral angles. *Journal of Chemical Theory and
> Computation*, 8(9), 3257-3273. https://doi.org/10.1021/ct300400x

---

## 4. Cite the solvent model you selected

### Implicit (generalised Born)

**OBC2 (Jiggle default) and OBC1**

> Onufriev, A., Bashford, D., & Case, D. A. (2004). Exploring protein native
> states and large-scale conformational changes with a modified generalized
> Born model. *Proteins*, 55(2), 383-394. https://doi.org/10.1002/prot.20033

**GBn**

> Mongan, J., Simmerling, C., McCammon, J. A., Case, D. A., & Onufriev, A.
> (2007). Generalized Born model with a simple, robust molecular volume
> correction. *Journal of Chemical Theory and Computation*, 3(1), 156-169.
> https://doi.org/10.1021/ct600085e

**GBn2**

> Nguyen, H., Roe, D. R., & Simmerling, C. (2013). Improved generalized Born
> solvent model parameters for protein simulations. *Journal of Chemical
> Theory and Computation*, 9(4), 2020-2034.
> https://doi.org/10.1021/ct3010485

**HCT**

> Hawkins, G. D., Cramer, C. J., & Truhlar, D. G. (1996). Parametrized models
> of aqueous free energies of solvation based on pairwise descreening of
> solute atomic charges from a dielectric medium. *The Journal of Physical
> Chemistry*, 100(51), 19824-19839. https://doi.org/10.1021/jp961710n

### Explicit water

**TIP3P** — and, as mTIP3P, the CHARMM water model

> Jorgensen, W. L., Chandrasekhar, J., Madura, J. D., Impey, R. W., & Klein,
> M. L. (1983). Comparison of simple potential functions for simulating
> liquid water. *The Journal of Chemical Physics*, 79(2), 926-935.
> https://doi.org/10.1063/1.445869

**TIP3P-FB** (Jiggle's default explicit water for the AMBER14/19 sets)

> Wang, L.-P., Martinez, T. J., & Pande, V. S. (2014). Building force fields:
> An automatic, systematic, and reproducible approach. *The Journal of
> Physical Chemistry Letters*, 5(11), 1885-1891.
> https://doi.org/10.1021/jz500737m

**SPC/E**

> Berendsen, H. J. C., Grigera, J. R., & Straatsma, T. P. (1987). The missing
> term in effective pair potentials. *The Journal of Physical Chemistry*,
> 91(24), 6269-6271. https://doi.org/10.1021/j100308a038

**TIP4P-Ew**

> Horn, H. W., Swope, W. C., Pitera, J. W., Madura, J. D., Dick, T. J.,
> Hura, G. L., & Head-Gordon, T. (2004). Development of an improved four-site
> water model for biomolecular simulations: TIP4P-Ew. *The Journal of Chemical
> Physics*, 120(20), 9665-9678. https://doi.org/10.1063/1.1683075

---

## 5. Algorithms, for methods sections that go into detail

**Particle Mesh Ewald** — used for all explicit-solvent electrostatics.

> Essmann, U., Perera, L., Berkowitz, M. L., Darden, T., Lee, H., &
> Pedersen, L. G. (1995). A smooth particle mesh Ewald method. *The Journal of
> Chemical Physics*, 103(19), 8577-8593. https://doi.org/10.1063/1.470117

**Kabsch superposition** — underlies every RMSD and RMSF value Jiggle reports,
including the fit-on-receptor pose-drift metric.

> Kabsch, W. (1976). A solution for the best rotation to relate two sets of
> vectors. *Acta Crystallographica Section A*, 32(5), 922-923.
> https://doi.org/10.1107/S0567739476001873

**NumPy** — all analysis arithmetic.

> Harris, C. R., Millman, K. J., van der Walt, S. J., *et al.* (2020). Array
> programming with NumPy. *Nature*, 585, 357-362.
> https://doi.org/10.1038/s41586-020-2649-2

Also used, and worth acknowledging rather than formally citing: **PySide6/Qt**
for the interface and **pyqtgraph** for the plots.

---

## 6. A worked example

For an implicit-solvent run on a receptor-binder complex with the defaults,
the citation set is:

Jiggle · OpenMM (Eastman 2017) · ff14SB (Maier 2015) · OBC2 (Onufriev 2004)
· 3Dmol.js if you used a viewer figure (Rego & Koes 2015)

Suggested methods text:

> Molecular dynamics was performed in Jiggle v0.1.0 (Safeer, 2026), a
> real-time front end for OpenMM 8 (Eastman *et al.*, 2017). Structures were
> prepared with PDBFixer: heterogens and crystallographic water were removed,
> missing atoms and residues rebuilt, and hydrogens added at pH 7.0. The
> system was parameterised with AMBER ff14SB (Maier *et al.*, 2015) and the
> OBC2 generalised Born implicit solvent model (Onufriev *et al.*, 2004),
> with hydrogen bond lengths constrained. Energy was minimised to
> convergence, then the system was heated to 300 K over six stages with
> harmonic positional restraints (1000 kJ mol⁻¹ nm⁻²) on protein heavy atoms,
> released before production. Production dynamics used a LangevinMiddle
> integrator at 300 K with a 1 ps⁻¹ friction coefficient and a 2 fs timestep,
> for N ns. Binder pose stability was assessed as the binder RMSD computed
> after superposition on the receptor, which separates displacement within
> the binding site from global motion and from the binder's own internal
> rearrangement.

Swap in the explicit-solvent sentence where relevant:

> The complex was solvated in a cubic box with 1.0 nm padding using the
> TIP3P-FB water model (Wang *et al.*, 2014), with long-range electrostatics
> treated by Particle Mesh Ewald (Essmann *et al.*, 1995) and a 1.0 nm
> real-space cutoff.

---

## 7. Note on scope

Jiggle's own contribution is the real-time interface and the analysis layer:
the four-way RMSD decomposition, the incremental per-residue RMSF, and the
live coupling of all of it to a running simulation. It is not a simulation
engine and should not be described as one. Citing it alongside OpenMM rather
than instead of OpenMM is both accurate and the norm for tools of this kind —
the same convention VMD, ChimeraX, MDAnalysis and gmx_MMPBSA are cited under.
