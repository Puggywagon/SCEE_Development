# Barostats

**Topic:** [[03 - Molecular Dynamics]]
**Tags:** #molecular-dynamics #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

Barostat controls pressure by allowing volume to fluctuate — gives NPT ensemble. Box dimensions are dynamical variables evolving alongside atomic positions.

**Berendsen barostat:** Rescales box dimensions and atomic positions. Good for equilibration but does not correctly sample NPT ensemble.

**Parrinello-Rahman barostat:** Production standard. Extended degrees of freedom for box vectors. Rigorously samples NPT ensemble including correct volume fluctuations. Often paired with Nosé-Hoover thermostat.

**Practical note for SCEE:** Equilibrate in NPT to get correct density at target state point, then switch to NVT for production configuration generation.

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
