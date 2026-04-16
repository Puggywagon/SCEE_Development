# The Self-Consistency Condition

**Topic:** [[08 - SCEE]]
**Tags:** #scee #concept
**Status:** #to-read

---

## AI Refresher

> *Conceptual scaffold — starting point before reading. Expand and refine from textbooks.*

Calculation converged when MM point charges used as input to QM/MM calculation are consistent with charges derived from resulting QM wavefunction.

**Iterative cycle:**
1. Start with initial MM charges (typically gas phase)
2. Perform QM/MM electrostatic embedding calculation on central molecule
3. Derive new charges from resulting QM wavefunction (ESP/RESP charges)
4. Use updated charges as MM charges for next iteration
5. Repeat until charges no longer change within convergence threshold

Charges derived from QM wavefunction fitted to reproduce molecular electrostatic potential — ESP or RESP charges.

Conceptually analogous to SCF procedure in Hartree-Fock — molecular orbitals iterated until consistent with the mean field they generate.

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
