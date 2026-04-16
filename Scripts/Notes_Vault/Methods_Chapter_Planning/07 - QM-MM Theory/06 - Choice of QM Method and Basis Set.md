# Choice of QM Method and Basis Set

**Topic:** [[07 - QM-MM Theory]]
**Tags:** #qmmm #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

**QM Methods:**
- **Hartree-Fock (HF):** Mean field electron-electron interactions. Fast but misses electron correlation. Underestimates dispersion, can give poor dipole moments.
- **DFT:** Workhorse of modern QC. Electron density as fundamental variable. Good accuracy/cost balance. Choice of functional matters — hybrid functionals like B3LYP common, dispersion-corrected functionals increasingly preferred.
- **Post-HF (MP2, CCSD(T)):** More accurate but expensive. Used for benchmarking.

**Basis sets:** Molecular orbitals as linear combinations of atom-centred Gaussian functions. Larger = more accurate but more expensive.
- Basis set size: minimal (STO-3G) → split valence (6-31G) → triple zeta (6-311G)
- Polarisation functions (*): allow electron density to distort from spherical symmetry. Important for accurate dipole moments.
- Diffuse functions (+): extend far from nucleus. Important for anions, excited states, outer electron density.

For dipole moment calculations: choice of basis set particularly important — dipole moment sensitive to accuracy of electron density. Directly relevant to SCEE work.

---

## Textbook Notes

### Dill & Bromberg
*Start with accessible conceptual treatment*

#dill-bromberg

---

### Allen & Tildesley
*Simulation specific detail*

#allen-tildesley

---

### Frenkel & Smit
*Rigorous treatment*

#frenkel-smit

---

### Other Sources

#paper

---

## My Thoughts & Observations

> *Your own reactions, connections, questions, disagreements. This is where your academic voice develops.*

---

## Key Equations & Figures

> *Flag equations and figures as you encounter them. Link to equation cards in the Equations folder.*

- [ ] 

---

## Connections to Other Sections

- [[]]

---

## Synthesis Status

- [ ] Notes gathered
- [ ] Thoughts developed
- [ ] Ready for synthesis
- [ ] Draft written
