# 02 - Statistical Ensembles

**Topic:** [[01 - Statistical Mechanics Foundations]]
**Tags:** #statmech #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

An ensemble defines which microstates are accessible and how probable each one is, based on what macroscopic quantities you hold fixed. Your choice of ensemble is a physical statement about what constraints your real system operates under.

**The main ensembles:**

- **Microcanonical (NVE)** — N, V, and E fixed. Every accessible microstate has equal probability. Most fundamental but least practical.
- **Canonical (NVT)** — N and V fixed, system exchanges energy with heat bath at temperature T. Workhorse of most simulations.
- **Grand Canonical (μVT)** — T and V fixed, particles can be exchanged with a reservoir. Important for phase equilibria.
- **Isothermal-Isobaric (NPT)** — T and P fixed, volume can fluctuate. Most physically realistic for condensed phase simulations.

When you use a thermostat in MD or specify an ensemble in MC, you are choosing which of these frameworks applies.

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

- [ ] [[Boltzmann Factor]] — probability weighting for canonical ensemble
- [ ] Schematic comparing ensemble constraints — useful figure to find

---

## Connections to Other Sections

- [[01 - Classical Statistical Mechanics]] — ensembles follow from microstate probability
- [[03 - Partition Functions and Thermodynamic Properties]] — partition function is defined per ensemble
- [[03 - Molecular Dynamics]] — thermostat choice determines ensemble
- [[04 - Monte Carlo]] — ensemble sampling is central to MC

---

## Synthesis Status

- [ ] Notes gathered
- [ ] Thoughts developed
- [ ] Ready for synthesis
- [ ] Draft written
