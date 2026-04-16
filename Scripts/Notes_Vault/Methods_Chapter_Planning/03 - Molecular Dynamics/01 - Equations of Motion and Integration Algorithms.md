# Equations of Motion and Integration Algorithms

**Topic:** [[03 - Molecular Dynamics]]
**Tags:** #molecular-dynamics #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

MD propagates Newton's equations of motion forward in time. Forces come from the negative gradient of the potential energy — force field gives forces, integrator moves atoms forward in time.

**Verlet/Leapfrog algorithm:** Approximates continuous trajectory as series of small steps. Leapfrog updates positions and velocities in staggered fashion — velocities at half-integer timesteps, positions at integer timesteps. Good energy conservation.

**Timestep choice:** Too large → overshoot, atoms overlap, simulation crashes. Too small → computationally wasteful. Limited by fastest motions — typically bond stretching involving hydrogen. ~1-2 fs standard for atomistic simulations with constrained bonds. Rule of thumb: timestep ~ one tenth of period of fastest motion.

Energy conservation over long simulations is a key diagnostic.

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
