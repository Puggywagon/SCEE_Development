# Section 4: Analysis and Data Saving

---

## Collective Goal

This section takes the per-configuration dipole values from Section 3 and applies the SCEE self-consistency procedure to calculate the liquid-phase dipole moment. It then aggregates results across all replicas, applies IQR-based outlier filtering, plots the distribution of liquid dipole moments, computes summary statistics, and saves everything to CSV. It also collects additional simulation observables (density, model dipole, model dielectric constant) for the final results table.

---

## Script Files Involved

| Version | Files |
|---------|-------|
| V2 | `Run_SCEE.py` (all analysis logic inline) |
| V3 | `Run_SCEE.py` (orchestration), `Simulations_Analysis.py` (analysis class), `Gaussian_Calculations.py` (SCEE fitting in `init4`) |

---

## Imports and Inputs

### Imports

**V2** — all in `Run_SCEE.py`:
```python
import numpy, pandas, matplotlib.pyplot, glob, os, argparse
```

**V3:**
```python
# Run_SCEE.py
import numpy, pandas, matplotlib.pyplot, glob, os

# Simulations_Analysis.py
import gromacs, numpy, os, pandas, matplotlib.pyplot, argparse
```

### Inputs to This Section

These values flow in from earlier sections:

| Variable | Source | Description |
|----------|--------|-------------|
| `df1['dipole_l/m/h']` | `Gaussian_Calculations.get_multipole_statistics()` | Per-config dipole at three charge scaling levels |
| `num1`, `num2`, `num3` | Computed in `Run_SCEE` (V2) / `init4` (V3) | The three actual scaled charge values used in each Gaussian run |
| `mu_Vacuum` | `Vacuum.get_multipole_statistics()` / `Gauss.init0()` | Self-consistent vacuum dipole |
| `PCM1` | `PCM1.get_multipole_statistics()` / `Gauss.init1()` | Liquid dipole under full dielectric screening |
| `PCM2` | `PCM2.get_multipole_statistics()` / `Gauss.init2()` | Liquid dipole under electronic-only screening |
| `gas_dipole` | User input | Experimental gas-phase dipole moment (Debye) |
| `dipole_model` | `get_dipole_model()` | Force-field vacuum dipole |
| `Model_Dipole`, `Epsilon` | `get_dipole_model_liquid()` | Force-field liquid dipole and model dielectric constant |

### User Inputs

**V2** — inline:
| Variable | Description | Example (Aniline) |
|----------|-------------|-------------------|
| `gas_dipole` | Experimental gas-phase dipole moment (D) | `1.49` |
| `IQR_K` | IQR multiplier for outlier detection | `2.0` (hardcoded in analysis block) |

**V3** — from `Settings.yml` (partially — `gas_dipole` not yet wired from YAML in current WIP):
| YAML Key | Description |
|----------|-------------|
| `State_Conditions.Di_Const` | Experimental dielectric constant |
| `State_Conditions.Ref_Ind` | Experimental refractive index |

---

## How This Section is Called in `Run_SCEE`

**V2 — all inline in the main script body:**
```python
# Per-configuration SCEE self-consistency loop
for index, row in df1.iterrows():
    # polynomial fit and root-finding to get muL_SCEE
    ...
df1['muL_SCEE'] = muL_list

# Per-configuration liquid dipole
for i in df1['muL_SCEE']:
    delta_mu = (mu_SCEE * (PCM1[0] / PCM2[0])) - mu_Vacuum[0]
    mu_liquid = delta_mu + gas_dipole

# Aggregation across replicas
IQR_K = 2.0
write_filtered_data(k=IQR_K)
plot_dipoledist(central, k=IQR_K)
plot_dipoledist_filtered(central, k=IQR_K)
# ... build collection DataFrame and save CSV
```

**V3 — split between `Run_SCEE.py` and `Simulations_Analysis.py`:**
```python
# SCEE fitting now happens inside Gauss.init4(), returns SCEE series

# Per-configuration liquid dipole (still in Run_SCEE.py)
for i in df1['muL_SCEE']:
    delta_mu = (mu_SCEE * (PCM1 / PCM2)) - mu_Vacuum
    mu_liquid = delta_mu + gas_dipole

# Aggregation
Analysis.write_filtered_data(k=IQR_K)
Analysis.plot_dipoledist(solvent, k=IQR_K)
Analysis.plot_dipoledist_filtered(solvent, k=IQR_K)
```

---

## How the Relevant Functions Work

### SCEE Self-Consistency Fitting

This is the mathematical core of the method. For each configuration, three dipole values are available from Gaussian (`dipole_l`, `dipole_m`, `dipole_h`), each computed at a different charge scaling level (`num1`, `num2`, `num3`). The self-consistent dipole is the value where the input charge (used to polarise the solvent shell) equals the output dipole (computed by QM).

The fitting procedure:
1. Fit a quadratic polynomial through the three points `(num_i, dipole_i)`:
   ```
   dipole = a * charge² + b * charge + c
   ```
2. The self-consistent solution satisfies `dipole = charge`, so solve:
   ```
   a * muL² + (b-1) * muL + c = 0
   ```
   using the quadratic formula (taking the physically meaningful root).
3. If the leading coefficient `a` is negligibly small (< 1×10⁻⁶), fall back to a linear fit and solve:
   ```
   muL = -c / (b - 1)
   ```
4. For the vacuum case (no charge scaling), the same procedure is applied to get `muG`, the self-consistent vacuum dipole.

**In V2:** This fitting is done inline in `Run_SCEE.py` for all four quantities (Vacuum, PCM1, PCM2, SCEE).

**In V3:** The SCEE fitting is done inside `Gaussian_Calculations.init4()`. The Vacuum, PCM1, and PCM2 fittings remain in `Run_SCEE.py` as part of the readback blocks for `init0/1/2`. Note that for single-molecule calculations (Vacuum, PCM1, PCM2), the three `dipole_l/m/h` values will be identical (the three charge-scaled runs only differ for the cluster SCEE case), so the fitting there effectively just recovers the single dipole value.

### Liquid Dipole Calculation

Once `muL_SCEE` is obtained for each configuration, the liquid-phase dipole is calculated:

```python
delta_mu = (mu_SCEE * (PCM1 / PCM2)) - mu_Vacuum
mu_liquid = delta_mu + gas_dipole
```

**Physical interpretation:**
- `mu_SCEE` is the QM dipole of the solute embedded in a cluster of solvent molecules with scaled charges.
- The `PCM1 / PCM2` ratio scales this by the ratio of the full dielectric response to the electronic-only response, correcting for the incomplete dielectric environment of the finite cluster.
- Subtracting `mu_Vacuum` removes the self-consistent vacuum contribution.
- Adding `gas_dipole` (the experimental gas-phase value) anchors the result to experiment.

This per-configuration `mu_liquid` is the primary observable of the simulation.

### `get_dipole_model_liquid()` (both versions, called during analysis)

Calls the GROMACS dipole tool on the production trajectory to get the force-field liquid dipole and model dielectric constant:
- `Model_Dipole` — the average dipole moment of the model in the liquid phase (from the trajectory).
- `Epsilon` — the model dielectric constant computed from dipole fluctuations.

These are used in the collection DataFrame and in the optional dielectric constant correction formula.

**V3 note:** Currently points to `replica_1` with a hardcoded path and glob pattern (`*K/*Bar`). The return statement has a bug — `Epsilon` is capitalised but `epsilon` is the local variable name.

### `Sim_Density()` (V3 only, `Simulations_Analysis`)

Extracts the average simulation density from the production run energy file using `gromacs.tools.Energy`. This is added to the per-replica results as `df['Simulation Density']` and averaged across replicas in the final collection. The V2 pipeline did not collect this.

**V3 note:** Currently has hardcoded path to `replica_1` and accepts `system_title, i, T, p` as arguments that are not used internally.

### IQR Outlier Filtering — `write_filtered_data()` and `_iqr_masks()`

Outlier detection uses the interquartile range (IQR) method:
```
lower bound = Q1 - k * IQR
upper bound = Q3 + k * IQR
```
where `k = 2.0` (hardcoded). Data points outside these bounds are flagged as outliers.

`write_filtered_data()` reads each `Dipole_Calculations.csv` across all replica/T/p directories and writes a corresponding `Dipole_Calculations_filtered.csv` where outlier rows have their `muL_SCEE`, `delta_mu`, and `mu_liquid` values set to `NaN` (rather than removing the rows entirely, preserving the config index for reference). A minimum of 5 finite data points is required to perform the filtering — if a file has fewer, it is skipped.

### Distribution Plotting — `plot_dipoledist()` and `plot_dipoledist_filtered()`

Both functions follow the same core plotting pattern:
1. Collect all `mu_liquid` values from across all replica CSV files.
2. Apply IQR outlier detection.
3. Plot a histogram of kept values (purple) and visible outliers (orange).
4. Overlay a fitted Gaussian (normal distribution) curve using the mean and standard deviation of the kept data (teal).
5. Save to both `.png` and `.pdf`.

**`plot_dipoledist`** uses simple axis limits based on the data range.

**`plot_dipoledist_filtered`** applies additional axis intelligence:
- Computes where the Gaussian curve drops below 0.1% of its peak value to set axis limits.
- Detects gaps in the histogram (bins of zero count wider than 3 bins) and trims the axis to avoid long empty tails.
- Extends the axis to show mild outliers within 2σ of the mean.

**Known plotting issues (not being fixed in current work):** The axis limit logic in `plot_dipoledist_filtered` can produce odd results in some cases. These will be addressed in a future update.

**V2 vs V3 output path difference:**
- V2 saves to `../../../../../../../Research/SCEE_Project/Results/Aromatics_Chapter/figures/` — a hardcoded absolute path relative to the run directory.
- V3 saves to `../figures/` — a relative path, more portable.

### Summary Statistics Collection

After the per-replica loop, a `collection` DataFrame is assembled with one row per molecule summarising the full run. This includes:

**Simulation parameters:** solvent name, temperatures, pressures, replicas, configuration counts (total and filtered).

**Charge scaling:** `qr1/2/3` ratios, `num1/2/3` actual scaled charges.

**Dipole moments:** experimental gas dipole, vacuum model dipole, liquid model dipole, `mu_Vacuum` (SCEE), `PCM1`, `PCM2`, `muL_SCEE` (mean), `delta_mu` (mean), `mu_liquid` (mean, unfiltered and filtered).

**Statistical measures:** Standard deviation and standard error computed from per-replica means (using `ddof=1` for sample standard deviation). Both unfiltered and filtered versions are computed.

**Dielectric constant correction (optional):** If `mu_liquid > Model_Dipole`:
```
epsilon_corrected = n² + (mu_liquid / Model_Dipole) * (Model_Epsilon + 1)
```
where `n` is the experimental refractive index. If the condition is not met, this is set to NaN.

**V3 addition:** Simulation density and model dielectric constant are now included in the collection, which V2 did not capture.

---

## Outputs

| File | Description |
|------|-------------|
| `Dipole_Calculations.csv` | Per-configuration results: `config`, `dipole_l/m/h`, `muL_SCEE`, `delta_mu`, `mu_liquid` (one per replica/T/p directory) |
| `Dipole_Calculations_filtered.csv` | Same as above with outlier rows NaN-ed out |
| `Dipole_Vacuum.csv` | Vacuum self-consistent dipole results |
| `Dipole_PCM1.csv` | PCM1 self-consistent dipole results |
| `Dipole_PCM2.csv` | PCM2 self-consistent dipole results |
| `Dipole_SCEE.csv` | Per-configuration SCEE dipole results |
| `Dipole.txt` | Raw GROMACS vacuum dipole output |
| `Dipole2.txt` | Raw GROMACS liquid dipole output |
| `Density.txt` (V3) | Raw GROMACS density output |
| `dipole_distribution_{Molecule}.png/.pdf` | Unfiltered dipole distribution plot |
| `dipole_distribution_{Molecule}_filtered.png/.pdf` | Filtered dipole distribution plot |
| `dipole_stats_{Molecule}.csv` | Summary statistics for unfiltered distribution |
| `dipole_stats_{Molecule}_filtered.csv` | Summary statistics for filtered distribution |
| `Results_{solvent}.csv` (V3) / `Results_{central}.csv` (V2) | Full summary results table |

---

## Changes from V2 to V3 and Why

| Change | Reason |
|--------|--------|
| All analysis functions moved from `Run_SCEE.py` into `Simulations_Analysis.py` class | Better organisation. V2 had analysis logic, plotting, and data collection all inline in the main script, making it difficult to reuse or modify independently. |
| SCEE self-consistency fitting moved into `Gaussian_Calculations.init4()` | Keeps the fitting immediately after the calculation that produces the data. |
| `Sim_Density()` added (V3) | Captures simulation density as an additional observable for validation against experiment. |
| Output paths changed from absolute to relative (`../figures/`) | Makes the script portable across machines and users. |
| `read_filtered_data()` and related functions moved to `Simulations_Analysis` | Consistent with the refactoring of all analysis into the class. |
| `combine_group()` added to `Oniom_Generation.py` (V3, misplaced) | This function computes a weighted mean and combined standard error across replica groups. It appears intended for `Simulations_Analysis.py` but has been placed in `Oniom_Generation.py` — likely a development artefact that will be moved. |
| `collection` DataFrame extended with density and epsilon | Additional observables for comparison with experimental data. |

---

## Other Notes

- **`gas_dipole` not yet in YAML (V3 WIP):** The experimental gas-phase dipole moment is a critical input for the `mu_liquid` formula but is not yet read from `Settings.yml` in V3. It must be added as a user-facing YAML key.
- **Replica mean vs configuration mean:** Standard deviation and standard error are computed from per-replica means (not from the full pool of configurations). This is the statistically correct approach — individual configurations within a replica are correlated, so treating them as independent samples would underestimate the true uncertainty.
- **IQR_K = 2.0 hardcoded:** The outlier threshold is fixed in both V2 and V3. A tighter value removes more configurations; a looser value retains more. For Aniline this was set to 2.0. This may need to be molecule-dependent and could be promoted to a YAML setting.
- **Dielectric constant correction formula:** The correction `n² + (mu_liquid / Model_Dipole) * (Model_Epsilon + 1)` is described in the comments as something to ask Miguel/Leo about — the physical interpretation of this correction and when it should be applied (only when `mu_liquid > Model_Dipole`) is not fully documented in the code.
- **Plot issues flagged for future work:** Both plotting functions can produce axis limits that don't look ideal for all distributions. The filtered plot in particular uses multi-step heuristics for axis limits (Gaussian tails, gap detection, mild outlier extension) that interact in ways that occasionally produce sub-optimal results. These are known and will be updated.
