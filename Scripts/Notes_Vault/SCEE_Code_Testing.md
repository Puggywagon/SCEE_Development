# SCEE Pipeline — Testing Protocol

A staged checkout of [[SCEE_Code_Overview]]. Each stage says what to run, what should appear, and what to
look at to decide it worked. Work down the list — a failure at stage *n* makes everything after it
meaningless, so there is no point skipping ahead.

> **Before starting:** the pipeline is a single script with no resume capability, so re-running from the top
> repeats the MD. Two ways to work stage by stage:
>
> - **Full-run mode** — run `Run_SCEE.py` on the reduced settings below and inspect at each checkpoint. Copy
>   the whole test directory after each successful stage (`cp -r Test Test_after_stage4`) so a failure later
>   doesn't cost you the earlier work.
> - **Interactive mode** — drive the classes directly from a Python session started *inside* the system
>   directory. Snippets are given per stage. This is the faster way to test one stage repeatedly.

---

## Tier 0 — Test system and reduced settings

Two tiers are worth running, and they answer different questions.

| Tier | System | Purpose | Cost |
|---|---|---|---|
| **A — plumbing** | any small molecule, heavily reduced settings | Does every stage run, hand off the right files, and produce output of the right shape? | minutes to a few hours |
| **B — validation** | ethanol, production settings | Does it reproduce a number we already trust? | full production cost |

**Tier A numbers are not physically meaningful.** With 5 configurations and a small box the dipole is
statistically worthless; the point is that the machinery works. Do not read anything into the value.

Suggested Tier A `Settings.yml` overrides:

```yaml
State_Conditions:
    Temperature: [298]
    Pressure: [1]
    Replicas: 1
Advanced_Settings:
    configuration:
        box_length_nm: 3.0        # ~200 molecules instead of ~1800
        cluster_radius_nm: 1.2    # must stay below half the equilibrated box length
    sampling:
        n_configurations: 5       # 5 configs x 3 charge sets = 15 Gaussian jobs
    filtering:
        min_kept_configs: 3       # otherwise the floor warning fires constantly
```

Production length in `template_2.mdp` also needs shortening to match — with `b=1000` and `dt=20`, five
configurations need roughly 1100 ps of production, not 5000.

- [ ] Reduced `Settings.yml` and `template_2.mdp` prepared in a scratch directory
- [ ] `Results/` subdirectory created (nothing creates it automatically)
- [ ] Topology, `.itp` and `.mdp` files copied in

---

## Stage 1 — Settings and derived quantities

Cheapest possible check, and it needs no MD at all.

```python
from Settings_Loader import load_settings
s = load_settings()
print(s.N, s.initial_molecules, s.molecule.calculated_dielectric)
print(s.is_mixture, s.build_solvent_gro, s.build_solute_gro)
```

- [ ] Loads without a `KeyError` (a missing YAML key fails here, not three hours in)
- [ ] `calculated_dielectric` = ε − n² + 1, checked by hand on a calculator
- [ ] `initial_molecules` is consistent with N = ρV/M × N_A. Worked example: propylamine at 6.2 nm gives
      V = 238 nm³ and N ≈ 1745, so `initial_molecules` ≈ 1830 after the 5 % buffer and the −1
- [ ] `cluster_radius_nm` is **less than half** the expected equilibrated box length — otherwise the
      minimum-image convention cuts a cluster that isn't a sphere. At the defaults, 2.8 nm against a 6.2 nm
      box leaves very little margin

**Fails if:** derived counts are orders of magnitude out (usually a units error in density or molar mass).

---

## Stage 2 — Molecule generation

```python
import Gro_Builder
gb = Gro_Builder.Gro_Builder()
gb.AA_Solvent('NCCC')
```

- [ ] `Solvent.gro` and `Solvent.png` written
- [ ] The `.png` depicts the intended molecule — this is the check the interactive prompt was written for
- [ ] Atom count in line 2 of the `.gro` matches the molecular formula
- [ ] Coordinates are in nm (values of order 0.1–1.0, not 1–10)
- [ ] Box line reads `0.0 0.0 0.0` — correct at this stage, `editconf` sets it next
- [ ] **Atom order matches `Solvent.itp` line for line.** This is the single most important manual check in
      the whole pipeline and nothing verifies it automatically. A mismatch produces plausible-looking
      numbers that are wrong

For the united-atom path, additionally:

- [ ] Aliphatic hydrogens on sp³ carbons are absent
- [ ] Order is heteroatoms → carbons → polar H → explicit H
- [ ] Atom count matches the UA `.itp`

**Fails if:** RDKit cannot embed the molecule (raises for the UA path, silently returns −1 in the AA path —
worth checking the `.gro` is not empty).

---

## Stage 3 — Vacuum equilibration and box filling

```python
import Gro_Simulations
MD = Gro_Simulations.Gro_Simulations()
MD.ntmpi, MD.ntomp, MD.pin, MD.pinoffset = 1, 8, 'on', 0
MD.run_md1(s, role='pure')
```

- [ ] `em.gro`, `nvt_vacuum2.gro`, `nvt_vacuum2.edr` present
- [ ] `nvt_vacuum2.gro` contains **velocity columns** — later stages assume nine fields per atom line
- [ ] The molecule is intact and not stretched across the box (visual check in VMD or PyMOL)
- [ ] `out.gro` exists, and its atom count equals `n_atoms_central + n_atoms_bath × initial_molecules`
- [ ] The bath count line was appended to `Pure.top` **exactly once** — re-running this stage appends again

Count check:

```bash
head -2 out.gro            # total atoms
grep -c "Solv" out.gro     # or count by residue
tail -3 Pure.top
```

**Fails if:** `insert-molecules` placed fewer molecules than requested. Nothing checks this, and it surfaces
much later as an atom-count mismatch. Compare the `out.gro` atom count against the number written into
`Pure.top` before going further.

---

## Stage 4 — Condensed-phase equilibration

```python
MD.run_md2(s, role='pure')
```

- [ ] `em1.gro`, `nvt_eq.gro`, `npt.gro` all present
- [ ] NPT density is within a few percent of the experimental `density` in the settings — **the first real
      physics check in the pipeline.** A density that is 10 % out means the topology or the force field is
      wrong, and every dipole downstream will be wrong with it
- [ ] Equilibrated box length is roughly `box_length_nm`, and comfortably more than twice
      `cluster_radius_nm`
- [ ] Energy and density traces are flat over the last part of the run, not still drifting

```bash
gmx energy -f npt.edr -s npt.tpr   # select Density, then Potential
```

---

## Stage 5 — `oniom.inp` generation

This file is small enough to check entirely by hand, and doing so once is worth a great deal.

```python
import Oniom_Generation
og = Oniom_Generation.Oniom_Generation()
og.Gen_File('oniom.inp', s.advanced.sampling.n_configurations,
            s.molecule.system_title, s.advanced.configuration.cluster_radius_nm)
total, qmax = og.QM_Inputs(s, role='pure', Oniom='oniom.inp')
og.MM_Inputs(s, Oniom='oniom.inp')
og.Counting_Molecules('oniom.inp', s.initial_molecules)
```

Open the file against [[SCEE_Code_Overview]] §7 and check:

- [ ] Header: configuration count, system title, cutoff radius, all as set
- [ ] QM atom count matches the `.itp`; central atom index points at the atom you expect (largest |q|, or the
      water oxygen)
- [ ] Element symbols are right — particularly any atom whose GROMACS name starts `CA`, `NA`, `CL`, `BR`
      or `SI`, since those are resolved by mass and are the ones that can silently resolve wrong
- [ ] Dummy sites appear as `Bq` with the flag set to 1
- [ ] σ values are the `.itp` values × 10 (nm → Å); ε values are unchanged from the `.itp` (kJ/mol)
- [ ] Any polar hydrogen with zero LJ in the `.itp` has picked up the hard-coded 0.2673 / 0.0418 floor
- [ ] Charges in column q1 equal the `.itp` charges (since `qr1 = 1.0`); q2 and q3 are those values × 1.35
      and × 1.7
- [ ] Each charge column sums to approximately zero across the molecule
- [ ] MM block matches the QM block for a pure system, except for the missing dummy column
- [ ] Final line equals `initial_molecules`, and matches what went into `Pure.top`

**Fails if:** an atom type is missing from the `.top` include tree — this raises a clear `ValueError` naming
the type.

---

## Stage 6 — Production and frame extraction

```python
import os
MD.run_md3(s, role='pure', HOMEDIR=os.getcwd(), T=298.0, p=1.0)
```

- [ ] `Prod.mdp` was written with the literal `TEMPERATURE` and `PRESSURE` tokens actually replaced
- [ ] `Pure_QMMM_md3.{tpr,trr,xtc,edr,gro}` present
- [ ] `conf_*.gro` files appear; count them and note the **first index**

```bash
ls conf_*.gro | head -3
ls conf_*.gro | wc -l
```

- [ ] There are at least `n_configurations` files numbered from `conf_1` upward. `trjconv -sep` numbers from
      zero, so confirm the run is long enough that `conf_1 … conf_N` all exist
- [ ] Every `conf_*.gro` has the same atom count as `out.gro`
- [ ] Molecules are whole, not split across the periodic boundary (`pbc=whole` should have handled it —
      check one frame visually)

---

## Stage 7 — Cluster extraction (single configuration first)

Run this on **one** configuration before letting it loose on the full set.

```python
from Shell_Oniom import parse_oniom_inp, process_configuration
setup = parse_oniom_inp('oniom.inp')
process_configuration(setup, 1)
```

- [ ] Exactly six files appear for that configuration: `q1`, `q2`, `q3` and their `_chg` counterparts
- [ ] The printed molecule count is physically sensible. Expected ≈ (4/3)πr³ × ρN_A/M — for propylamine at
      r = 2.8 nm that is roughly 670 bath molecules, i.e. an ONIOM job of ~8700 atoms
- [ ] The total molecule count (QM + bath) is **odd** — if `REMOVED` was printed, the even-count workaround
      fired correctly
- [ ] Coordinates are in Å and centred on the QM molecule: the central atom's line should be near
      `0.000 0.000 0.000`, and no coordinate should exceed roughly 10 × the cutoff in nm
- [ ] No coordinate is wildly out of range, which would mean the minimum-image correction failed
- [ ] The `q1`, `q2`, `q3` files are identical except for the charge fields
- [ ] `_chg` files contain only the bath atoms, one line each, and their charges sum to approximately zero
- [ ] The `VDW` block at the foot lists every QM and MM atom type once, with σ ≈ the `oniom.inp` value ×
      0.561 and ε ≈ the `oniom.inp` value ÷ 4.184
- [ ] Connectivity indices run 1 … n_atoms with no gaps and no repeats

**Fails if:** the atom count in `conf_1.gro` disagrees with `oniom.inp` — this raises immediately and
usually means the insertion shortfall from stage 3.

---

## Stage 8 — Reference Gaussian jobs

```python
import Gaussian_Calculations
G = Gaussian_Calculations.Gaussian_Calculations(s)
mu_vac = G.init0()
pcm1 = G.init1()
pcm2 = G.init2()
```

- [ ] `Vacuum/`, `PCM1/`, `PCM2/` each contain a `.dat` and a `.out`
- [ ] Every `.out` ends with `Normal termination` (twice — these are two-step `--Link1--` jobs)
- [ ] `Dipole_Vacuum.csv`, `Dipole_PCM1.csv`, `Dipole_PCM2.csv` written
- [ ] μ_Vacuum is close to the experimental gas dipole. Not identical — it is a calculated value at the v0
      level — but the same order and within a few tenths of a Debye
- [ ] Ordering is μ_Vacuum < PCM2 < PCM1. A continuum always polarises the molecule upward, and PCM2 uses a
      lower ε than PCM1, so a different ordering means something is wrong with the ε values
- [ ] The PCM1/PCM2 ratio is modest, of order 1.0–1.1. A large ratio propagates straight into Δμ
- [ ] The geometry in the `.dat` files looks like the molecule (`gro_to_dat` recentres and converts to Å —
      a units slip here shows up as absurd bond lengths)

**Fails if:** `g09root` is wrong — the generated `run_.sh` fails to source the profile and the `.out` files
are empty or missing. Inspect `run_.sh` directly; it is a plain bash script.

---

## Stage 9 — SCEE jobs and the self-consistency fit

Run stage 7 for all configurations first (`ShellOniom().process('oniom.inp')`), then:

```python
G.natom = total     # from stage 5
G.qmax = qmax
G.init3()
scee_df = G.init4()
print(scee_df)
```

- [ ] `SCEE/` contains three `.dat`/`.out` pairs per configuration for `init3`, and three `_v1` pairs for
      `init4`
- [ ] All terminate normally; time one job before launching the batch, and multiply
- [ ] The throttle works: `pgrep -c g09` during the run should never exceed `max_jobs`
- [ ] `Dipole_SCEE.csv` has one row per configuration with `dipole_l`, `dipole_m`, `dipole_h`, `muL_SCEE`
- [ ] Within each row, `dipole_l < dipole_m < dipole_h` — stronger embedding charges must polarise the
      molecule more. Any row breaking this ordering indicates an SCF that converged somewhere odd
- [ ] `muL_SCEE` is larger than μ_Vacuum and of a plausible magnitude for the liquid
- [ ] Row count equals the number of configurations. A shortfall means some jobs produced no parseable
      dipole and were dropped silently by the finite-value filter — check for crashed jobs before accepting
      the result

Worth instrumenting once: how often the quadratic fit falls back to the linear branch. Frequent fallbacks
mean the three `qr` values are too close together to define the curvature.

---

## Stage 10 — Analysis and collation

- [ ] `Dipole_Calculations.csv` written per leaf, with `delta_mu` and `mu_liquid` columns
- [ ] Spot-check one row by hand: `delta_mu = muL_SCEE × (PCM1/PCM2) − mu_Vacuum`, and
      `mu_liquid = delta_mu + Gas_Dip`
- [ ] `Density.txt` and `Dipole2.txt` present, and the simulated density again matches experiment
- [ ] ΔH_vap is positive and of a sensible magnitude for the molecule
- [ ] `Results/dipole_stats.csv` contains `raw`, `pass1` … `filtered` rows for the state point. **Check it
      was not appended to from a previous run** — see the filename-casing item in [[SCEE_Code_Tasks]]
- [ ] The distribution plot renders, with kept and removed populations distinguishable
- [ ] `Results.csv` has one row per state point; `Per_Replica.csv` one row per replica
- [ ] `N_configs_kept` is a sensible fraction of `N_configs_total`. Heavy filtering means either a genuinely
      broad distribution or a problem upstream
- [ ] `below_floor` is `False`
- [ ] The settings snapshot in `Results/` matches what was actually run

---

## Stage 11 — Validation run

The only stage that tests whether the answer is *right* rather than merely produced.

- [ ] Full production settings, ethanol, one replica at 298 K and 1 bar
- [ ] μ_liquid agrees with the previously obtained value for the same force field
- [ ] Simulated density agrees with experiment
- [ ] The filtered configuration count is comparable to previous runs
- [ ] A second replica gives a mean within the quoted standard error of the first

If this reproduces, the pipeline is behaving. If Tier A passed and this does not, the problem is physical
or statistical rather than a plumbing fault, which is a much narrower place to look.

---

## Quick reference — what each stage hands to the next

| Stage | Hands forward |
|---|---|
| 2 | `Solvent.gro` |
| 3 | `nvt_vacuum2.gro` (all reference QM geometry), `out.gro`, topology molecule count |
| 4 | `npt.gro` |
| 5 | `oniom.inp`, plus `Total_Atoms` and `qmax` needed by stage 9 |
| 6 | `conf_*.gro`, `Pure_QMMM_md3.trr/.edr` |
| 7 | six `.inp` files per configuration |
| 8 | μ_Vacuum, PCM1, PCM2 |
| 9 | `muL_SCEE` per configuration |
| 10 | `Results.csv`, `Per_Replica.csv` |

Break any one of these links and the failure appears at the next stage, not the one that caused it — which
is the reason for testing in this order.
