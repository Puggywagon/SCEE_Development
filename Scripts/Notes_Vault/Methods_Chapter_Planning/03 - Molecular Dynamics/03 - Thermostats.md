# Thermostats

**Topic:** [[03 - Molecular Dynamics]]
**Tags:** #molecular-dynamics #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

Plain Newtonian MD conserves total energy — samples NVE ensemble. Thermostat couples simulation to heat bath at target temperature for NVT.

**Velocity rescaling:** Periodically scale all velocities. Simple but doesn't strictly sample canonical ensemble.

**Berendsen thermostat:** Exponentially relaxes temperature towards target with time constant. Stable, good for equilibration. Suppresses temperature fluctuations — does not correctly sample canonical ensemble. Not recommended for production runs.

**Nosé-Hoover thermostat:** Gold standard for production NVT. Extended degree of freedom representing heat bath. Rigorously generates canonical ensemble including correct temperature fluctuations.

**Bussi velocity rescaling:** Combines Berendsen stability with correct canonical ensemble sampling via stochastic element. Popular modern choice.

For SCEE: Nosé-Hoover or Bussi rescaling appropriate for production MD runs.

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
