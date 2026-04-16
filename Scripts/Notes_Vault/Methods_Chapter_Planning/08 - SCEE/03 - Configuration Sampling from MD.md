# Configuration Sampling from MD

**Topic:** [[08 - SCEE]]
**Tags:** #scee #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

Single QM/MM calculation gives dipole moment of one molecule in one instantaneous configuration. Dipole moment is a thermodynamic property — should represent average over full ensemble.

MD trajectory provides configurations sampling canonical ensemble at target T and P. Extract configurations at intervals sufficient to ensure statistical independence — too close together risks correlated structures.

**Number of configurations:** Balance between statistical convergence and computational cost. Monitor convergence of average dipole moment as configurations added.

**Size of MM region:** Must be large enough to capture full electrostatic influence of liquid on central molecule. Common: include all molecules within defined cutoff radius. Check results converged with respect to this radius.

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
