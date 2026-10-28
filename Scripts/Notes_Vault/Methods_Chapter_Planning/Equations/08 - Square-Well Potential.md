# Equation — Square-Well Potential

**Topic:** See connections below
**Tags:** #force-fields #equation #zenoline
**Status:** #to-verify

---

## Equation

$$
u(r) = \begin{cases} +\infty & \text{for } 0 < r < \sigma \\ -\varepsilon & \text{for } \sigma < r < \lambda\sigma \\ 0 & \text{for } \lambda\sigma < r \end{cases}
$$

---

## Term Breakdown

| Symbol | Meaning | Units |
|--------|---------|-------|
| $u(r)$ | Pair interaction energy | J or reduced units |
| $r$ | Separation between particle centres | m or reduced units |
| $\sigma$ | Hard sphere diameter | m or reduced units |
| $\varepsilon$ | Well depth | J or reduced units |
| $\lambda$ | Well width parameter (well extends to $\lambda\sigma$) | Dimensionless |

---

## Physical Interpretation

Highly simplified pair potential — hard core repulsion at short range, constant attractive well of fixed width and depth beyond that. Analytically tractable — exact virial expansion to third coefficient. Captures essential physics of short range attraction and hard core repulsion in simplest possible form. Central to Zeno line project alongside Mie potential.

---

## Sources

| Book / Paper | Chapter/Section | Notes |
|------|----------------|-------|
| Dill & Bromberg |  |  |
| Allen & Tildesley |  |  |
| Frenkel & Smit |  |  |
| Paterson, Bannerman, Lue (2024) | J. Chem. Phys. 160, 154503 | Group paper — key reference |

---

## Connections

- [[Lennard-Jones Potential]]
- [[01 - Pair Potentials]]
- [[03 - The Virial Expansion]]
