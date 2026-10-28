# Equation — Mie Potential

**Topic:** See connections below
**Tags:** #force-fields #equation #zenoline
**Status:** #to-verify

---

## Equation

$$
u(r) = \frac{n\varepsilon}{n-m}\left(\frac{n}{m}\right)^{m/(n-m)} \left[\left(\frac{\sigma}{r}\right)^{n} - \left(\frac{\sigma}{r}\right)^{m}\right]
$$

---

## Term Breakdown

| Symbol | Meaning | Units |
|--------|---------|-------|
| $u(r)$ | Pair interaction energy | J or reduced units |
| $r$ | Separation between particle centres | m or reduced units |
| $\varepsilon$ | Interaction energy (well depth) | J or reduced units |
| $\sigma$ | Particle diameter | m or reduced units |
| $n$ | Repulsive exponent (often varied, LJ uses n=12) | Dimensionless |
| $m$ | Attractive exponent (typically m=6 for physical basis) | Dimensionless |

---

## Physical Interpretation

Generalisation of the Lennard-Jones potential with adjustable exponents. Setting n=12, m=6 recovers the LJ 12-6 potential. The repulsive exponent n controls the steepness of the repulsive wall. Central to the Zeno line project — varying n and examining effect on Zeno line linearity.

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
- [[05 - The Zeno Line — Definition and Properties]]
- [[06 - The Zeno Line and Corresponding States]]
