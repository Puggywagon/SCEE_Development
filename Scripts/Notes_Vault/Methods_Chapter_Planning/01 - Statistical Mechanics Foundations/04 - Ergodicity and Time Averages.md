# 04 - Ergodicity and Time Averages

**Topic:** [[01 - Statistical Mechanics Foundations]]
**Tags:** #statmech #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

In practice you cannot sum over all microstates — there are far too many. In MD you follow the system through time and average over the trajectory. The question is: is a **time average** the same as an **ensemble average**?

The **ergodic hypothesis** says yes — provided you run long enough, the system will visit all accessible microstates, and the time average will equal the ensemble average.

In practice this is an assumption that can break down:
- If the system gets trapped in one region of configuration space
- Poor sampling or slow relaxation
- Phase transitions or metastable states

Checking for ergodicity problems is an important part of validating a simulation. It is closely related to the equilibration problem — a system that has not properly equilibrated has not sampled the full ensemble.

*Note: The Boltzmann-Ehrenfest connection — Boltzmann's entropy formula S = k ln W is on his tombstone. Both Boltzmann (1906) and Ehrenfest (1933) died having worked on statistical mechanics — the dark joke of the field.*

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

- [ ] [[Boltzmann Entropy Formula]] — S = k ln W
- [ ] [[Liouville Equation]] — phase space and ergodicity
- [ ] Time average vs ensemble average — schematic or formal equivalence

---

## Connections to Other Sections

- [[03 - Partition Functions and Thermodynamic Properties]] — ergodicity justifies computing ensemble averages from trajectories
- [[04 - Energy Minimisation and Equilibration]] — poor equilibration leads to ergodicity problems
- [[03 - Molecular Dynamics]] — MD relies on ergodicity for time averaging

---

## Synthesis Status

- [ ] Notes gathered
- [ ] Thoughts developed
- [ ] Ready for synthesis
- [ ] Draft written
