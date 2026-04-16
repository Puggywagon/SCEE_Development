# Equation — Metropolis Acceptance Criterion

**Topic:** See connections below
**Tags:** #monte-carlo #equation
**Status:** #to-verify

---

## Equation

$$
P_{acc}(m \rightarrow n) = \min\left(1, e^{-\beta \Delta U}\right) = \min\left(1, e^{-(U_n - U_m)/k_BT}\right)
$$

---

## Term Breakdown

| Symbol | Meaning | Units |
|--------|---------|-------|
| $P_{acc}$ | Acceptance probability for proposed move from state m to n | Dimensionless |
| $\Delta U$ | Change in potential energy $U_n - U_m$ | J |
| $U_m, U_n$ | Potential energy of current and proposed states | J |
| $\beta$ | Inverse temperature $1/k_BT$ | J⁻¹ |
| $k_B$ | Boltzmann constant | J K⁻¹ |
| $T$ | Temperature | K |

---

## Physical Interpretation

If the proposed move lowers the energy (ΔU < 0), it is always accepted (probability = 1). If it raises the energy (ΔU > 0), it is accepted with probability given by the Boltzmann factor of the energy difference. This asymmetry drives the system towards low energy configurations while allowing occasional uphill moves to escape local minima. Satisfies detailed balance and therefore generates the correct canonical ensemble.

---

## Sources

| Book / Paper | Chapter/Section | Notes |
|------|----------------|-------|
| Dill & Bromberg |  |  |
| Allen & Tildesley | Chapter 4 | Classic treatment |
| Frenkel & Smit | Chapter 3 | Rigorous derivation including detailed balance |

---

## Connections

- [[Boltzmann Factor]]
- [[02 - Metropolis Algorithm and Detailed Balance]]
- [[03 - Ensemble Sampling]]
