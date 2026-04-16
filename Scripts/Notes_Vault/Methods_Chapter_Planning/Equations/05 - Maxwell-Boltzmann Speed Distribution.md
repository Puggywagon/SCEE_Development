# Equation — Maxwell-Boltzmann Speed Distribution

**Topic:** See connections below
**Tags:** #molecular-dynamics #statmech #equation
**Status:** #to-verify

---

## Equation

$$
f(v) = 4\pi \left(\frac{m}{2\pi k_B T}\right)^{3/2} v^2 \exp\left(-\frac{mv^2}{2k_B T}\right)
$$

---

## Term Breakdown

| Symbol | Meaning | Units |
|--------|---------|-------|
| $f(v)$ | Probability density of speed v | s m⁻¹ |
| $m$ | Particle mass | kg |
| $k_B$ | Boltzmann constant | J K⁻¹ |
| $T$ | Temperature | K |
| $v$ | Speed | m s⁻¹ |

---

## Physical Interpretation

The probability distribution of molecular speeds in a gas at thermal equilibrium. Used to assign initial velocities in MD simulations. Follows from the equipartition theorem — each quadratic degree of freedom carries average kinetic energy ½k_BT.

---

## Sources

| Book / Paper | Chapter/Section | Notes |
|------|----------------|-------|
| Dill & Bromberg | Kinetic theory chapters | Very intuitive treatment |
| Allen & Tildesley | Chapter 3 — initial conditions |  |
| Frenkel & Smit |  |  |

---

## Connections

- [[02 - Initial Conditions]]
- [[Boltzmann Factor]]
