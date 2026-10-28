# Equation — ONIOM Energy Expression

**Topic:** See connections below
**Tags:** #qmmm #scee #equation
**Status:** #to-verify

---

## Equation

$$
E_{\text{ONIOM}} = E_{\text{high, model}} + E_{\text{low, real}} - E_{\text{low, model}}
$$

---

## Term Breakdown

| Symbol | Meaning | Units |
|--------|---------|-------|
| $E_{\text{high, model}}$ | QM energy of the inner (model/QM) region | J or Hartree |
| $E_{\text{low, real}}$ | MM energy of the entire system | J or kcal/mol |
| $E_{\text{low, model}}$ | MM energy of the inner region alone — subtracted to avoid double counting | J or kcal/mol |

---

## Physical Interpretation

The ONIOM subtractive energy scheme. Rather than computing coupling between regions explicitly, ONIOM constructs total energy by taking MM energy of whole system, adding QM energy of inner region, and subtracting MM energy of inner region to avoid double counting. Coupling handled implicitly through subtraction. Elegant and flexible.

---

## Sources

| Book / Paper | Chapter/Section | Notes |
|------|----------------|-------|
| Dill & Bromberg |  | Not covered — QM/MM not in scope |
| Allen & Tildesley |  | Not covered |
| Frenkel & Smit |  | Not covered |
| Gaussian documentation | ONIOM method | Primary reference |
| Morokuma group papers |  | Original development |

---

## Connections

- [[05 - The ONIOM Formalism]]
- [[04 - Electrostatic Embedding]]
- [[08 - SCEE]]
