# Multicanonical MC

**Topic:** [[04 - Monte Carlo]]
**Tags:** #monte-carlo #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

When system has free energy barrier between states (e.g. liquid/vapour near coexistence), standard MC spends almost all time in one phase.

Multicanonical MC modifies sampling weights so simulation visits all energies with equal probability — effectively removing the free energy barrier. Iteratively determines weights that flatten the energy histogram.

Raw multicanonical data reweighted using histogram reweighting to recover correct thermodynamic properties at any state point.

Multicanonical sampling + histogram reweighting: complementary tools. Multicanonical ensures you visit relevant regions of configuration space, histogram reweighting extracts thermodynamic information efficiently.

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
