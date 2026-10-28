# Equation — Truncated and Shifted LJ Potential

**Topic:** See connections below
**Tags:** #force-fields #equation
**Status:** #to-verify

---

## Equation

$$
u_{\text{shift}}(r) = \begin{cases} u_{\text{LJ}}(r) - u_{\text{LJ}}(r_c) & r \leq r_c \\ 0 & r > r_c \end{cases}
$$

---

## Term Breakdown

| Symbol | Meaning | Units |
|--------|---------|-------|
| $u_{\text{shift}}(r)$ | Truncated and shifted potential | J or reduced units |
| $u_{\text{LJ}}(r)$ | Standard Lennard-Jones potential | J or reduced units |
| $r_c$ | Cutoff distance | m or reduced units |

---

## Physical Interpretation

Shifts the LJ potential so it goes smoothly to zero at the cutoff, avoiding the discontinuity that arises from simple truncation. The discontinuity in simple truncation introduces systematic errors in energy and can cause numerical instability. Shifting removes this but changes the potential slightly — tail corrections needed to recover the full LJ thermodynamics.

---

## Sources

| Book / Paper | Chapter/Section | Notes |
|------|----------------|-------|
| Allen & Tildesley | Chapter 2 — force truncation | Standard reference |
| Frenkel & Smit |  |  |

---

## Connections

- [[Lennard-Jones Potential]]
- [[03 - Cutoffs]]
- [[04 - Long-Range Corrections]]
