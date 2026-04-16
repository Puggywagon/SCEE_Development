# Equation — Lennard-Jones Potential

**Topic:** See connections below
**Tags:** #force-fields #equation
**Status:** #to-verify

---

## Equation

$$
u(r) = 4\varepsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]
$$

---

## Term Breakdown

| Symbol | Meaning | Units |
|--------|---------|-------|
| $u(r)$ | Pair interaction energy as function of separation | J or reduced units |
| $r$ | Separation between particle centres | m or reduced units |
| $\varepsilon$ | Well depth — controls strength of interaction | J or reduced units |
| $\sigma$ | Collision diameter — sets length scale | m or reduced units |
| $r^{-12}$ | Repulsive term — Pauli repulsion | — |
| $r^{-6}$ | Attractive term — London dispersion | — |

---

## Physical Interpretation

The workhorse pair potential of molecular simulation. Captures both the short range Pauli repulsion (r⁻¹²) and the longer range London dispersion attraction (r⁻⁶). Two parameters: epsilon (well depth) and sigma (collision diameter). Despite its simplicity captures qualitative behaviour of noble gases and simple molecules remarkably well.

---

## Sources

| Book / Paper | Chapter/Section | Notes |
|------|----------------|-------|
| Dill & Bromberg | Intermolecular forces chapters |  |
| Allen & Tildesley | Chapter 1 | Extensive treatment |
| Frenkel & Smit | Chapter 2 |  |

---

## Connections

- [[Mie Potential]]
- [[01 - Pair Potentials]]
- [[05 - The Zeno Line — Definition and Properties]]
