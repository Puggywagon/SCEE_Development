# The ONIOM Formalism

**Topic:** [[07 - QM-MM Theory]]
**Tags:** #qmmm #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

ONIOM — Our own N-layered Integrated molecular Orbital and molecular Mechanics. Developed by Morokuma. Implemented in Gaussian. Flexible framework for combining different levels of theory.

**Subtractive energy scheme:**
E(ONIOM) = E(high, model) + E(low, real) - E(low, model)

Where high = QM method, low = MM method, real = entire system, model = QM region alone. Takes MM energy of whole system, adds QM energy of inner region, subtracts MM energy of inner region to avoid double counting. Coupling handled implicitly through subtraction.

For SCEE work: in-house ONIOM scripts automate input preparation and output processing across many configurations from MD trajectory.

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
