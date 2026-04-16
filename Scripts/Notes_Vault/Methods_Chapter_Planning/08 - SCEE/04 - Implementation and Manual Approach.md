# Implementation and Manual Approach

**Topic:** [[08 - SCEE]]
**Tags:** #scee #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

For each configuration from MD trajectory:
1. Identify central QM molecule, extract MM environment within cutoff
2. Prepare Gaussian/ONIOM input file with correct QM/MM specification
3. Run QM/MM calculation, extract ESP/RESP charges from wavefunction
4. Update MM charges for surrounding molecules
5. Repeat inner loop to convergence
6. Extract dipole moment from final wavefunction
7. Repeat for every configuration in ensemble
8. Average over all configurations

Manual approach makes methodology transparent and easy to validate. But practical burden significant — potentially hundreds of configurations, each requiring multiple QM/MM iterations. Scope for human error considerable. Time investment substantial even for single state point or molecule.

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
