# Section 2: Using the Vacuum Phase

---

## Collective Goal

This section prepares the isolated solute molecule for use in the quantum chemical calculations that form the core of the SCEE method. It covers: equilibrating the single solute molecule in vacuo using GROMACS MD; extracting its vacuum-phase geometry; and running three Gaussian09 calculations — a vacuum geometry optimisation/single-point (Vacuum), and two polarisable continuum model (PCM) single-points using the experimental dielectric constant (PCM1) and the optical (calculated) dielectric constant (PCM2). Together these three calculations provide the reference dipole moments needed to correct the SCEE liquid-phase dipole for the dielectric environment. The vacuum-phase model dipole moment is also extracted from GROMACS at this stage for use in the charge scaling later.

---

## Script Files Involved

| Version | Files |
|---------|-------|
| V2 | `Pre_Eq_Simulations.py`, `Vacuum.py`, `PCM1.py`, `PCM2.py`, `dat_generation.py` |
| V3 | `Gro_Simulations.py`, `Gaussian_Calculations.py`, `Atoms.py`, `Simulations_Analysis.py` |

---

## Imports and Inputs

### Imports

**V2:**
```python
# Pre_Eq_Simulations.py / MD classes
import gromacs, subprocess, numpy, os, socket, time

# Vacuum.py / PCM1.py / PCM2.py
import glob, gromacs, pandas, re, os
from subprocess import check_call
from dat_generation import Dat_Generation

# dat_generation.py
import pandas, re, numpy
```

**V3:**
```python
# Gro_Simulations.py
import gromacs, subprocess, numpy, os, socket, time

# Gaussian_Calculations.py
import glob, gromacs, pandas, re, os, numpy
import Atoms
import Simulations_Analysis

# Atoms.py
import re

# Simulations_Analysis.py
import gromacs, numpy, os, pandas, argparse
```

### User Inputs

**V2** — inline variables in `Run_SCEE.py`:
| Variable | Description | Example (Aniline) |
|----------|-------------|-------------------|
| `central` | Name of the solute molecule | `'Aniline'` |
| `L` | Target simulation box side length (nm) | `6.2` |
| `sol_keyword` | Gaussian solvent keyword for PCM | `'Aniline'` |
| `exp_diconst` | Experimental dielectric constant | `6.890` |
| `ref_ind` | Experimental refractive index | `1.586` |

**V3** — from `Settings.yml`:
| YAML Key | Description |
|----------|-------------|
| `State_Conditions.Di_Const` | Experimental dielectric constant |
| `State_Conditions.Ref_Ind` | Experimental refractive index |
| `State_Conditions.g09root` | Path to Gaussian09 installation |
| `State_Conditions.GAUSS_SCRDIR` | Path to Gaussian scratch directory |
| `State_Conditions.max_jobs` | Maximum concurrent Gaussian + MD jobs |
| `Advanced_Settings.configuration.box_length_nm` | Box length L |

`sol_keyword` and `exp_diconst` are referenced in V3 `Run_SCEE.py` but are not yet read from the YAML — this is a known gap in the current WIP version.

### Input Files (Template Folder)
| File | Description |
|------|-------------|
| `minim_vacuum.mdp` | GROMACS MDP parameters for vacuum energy minimisation |
| `nvt_vacuum.mdp` | GROMACS MDP parameters for vacuum NVT equilibration |
| `{central}_AA.top` (V2) / `{Topology_File}` (V3) | GROMACS topology file for the solute |
| `{central}_AA.gro` (V2) / `{Gro_File}` (V3) | Starting geometry from Section 1 (or user-provided) |

---

## How This Section is Called in `Run_SCEE`

**V2:**
```python
Pre_Eq_Simulations.Pre_Eq_Solute(central, L)
dipole_model = Pre_Eq_Simulations.get_dipole_model()

Vacuum.init_v0(system_title, sol_keyword)
Vacuum.run_gaussian(step=0)
df4 = Vacuum.get_multipole_statistics()

PCM1.init_v0(system_title, sol_keyword, exp_diconst)
PCM1.run_gaussian(step=0)
df2 = PCM1.get_multipole_statistics()

PCM2.init_v0(system_title, sol_keyword, cal_diconst)
PCM2.run_gaussian(step=0)
df3 = PCM2.get_multipole_statistics()
```
Note: in the V2 Aniline run, the `Pre_Eq_Solute`, `Vacuum.init_v0`, `PCM1.init_v0`, and `PCM2.init_v0` calls are all commented out, meaning the script reads back results from already-completed runs rather than re-running them. This is a common workflow pattern — run the heavy computations once, then comment them out and just read the outputs.

**V3:**
```python
if Box_Build == 'Yes':
    MD.run_md1(Gro_File, Topology_File, L, initial_molecules, solresnametop, Mixture='No', MD='Vacuum')
    dipole_model = Analysis.get_dipole_model()

mu_Vacuum = Gauss.init0(sol_keyword, Gaus='Vacuum')
PCM1 = Gauss.init1(sol_keyword, exp_diconst, Gaus='PCM1')
PCM2 = Gauss.init2(sol_keyword, cal_diconst, Gaus='PCM2')
```
In V3, `run_md1` is only called when `Box_Build == 'Yes'` (i.e., when the user has not provided a pre-equilibrated box). The Gaussian calls are all wrapped into methods of the single `Gaussian_Calculations` class.

---

## How the Relevant Functions Work

### MD Steps — `Pre_Eq_Simulations.Pre_Eq_Solute()` (V2) / `Gro_Simulations.run_md1()` (V3)

These functions equilibrate the single solute molecule in vacuum before it is inserted into the solvent box.

**Sequence of steps:**
1. The box size is set to `L + 1.7` nm using `gromacs.editconf`. The extra 1.7 nm is added as a buffer to avoid the molecule interacting with its periodic image.
2. A vacuum energy minimisation is run using the `minim_vacuum.mdp` parameters (`gromacs.grompp` + `gromacs.mdrun`), producing `em.gro`.
3. A vacuum NVT equilibration is run using `nvt_vacuum.mdp`, producing `nvt_vacuum2.gro`. This is the file that the Gaussian input generation reads from.

**V2 threading control:** `Pre_Eq_Simulations` contains a `_wait_for_md_slot()` method that checks the number of running `mdrun` and `g09` processes before launching any new MD job, keeping the total below `max_heavy_jobs`. This throttling was configured based on hostname (`Tower3` vs office PC).

**V3 change:** The `_wait_for_md_slot()` throttle has been removed from `Gro_Simulations`. Job throttling in V3 is only handled within the Gaussian bash runner scripts.

Additionally, in V3 `run_md1` calls `insert_Molecules` and `Write_to_top` at the end, combining the vacuum prep and molecule insertion into a single method call. In V2 these were separate calls.

### `Simulations_Analysis.get_dipole_model()` (both versions)

After the vacuum NVT run, the GROMACS dipole tool is used to extract the model's vacuum dipole moment from the `nvt_vacuum2` trajectory:
1. Calls `gromacs.tools.Dipoles` on the `nvt_vacuum2` edr/tpr files.
2. Writes the raw output to `Dipole.txt`.
3. Parses line 8 of the output to extract the average dipole moment value.
4. Returns `dipole_model` — the force-field dipole moment of the isolated molecule.

This value is stored and later used in the charge scaling calculation: `ratio = Model_Dipole_liquid / qmax`, which scales the three charge sets proportionally to the liquid-phase enhancement.

### Coordinate Conversion — `dat_generation.Dat_Generation.gro_to_dat()` (V2) / `Gaussian_Calculations.gro_to_dat()` (V3)

Before a Gaussian calculation can be set up, the coordinates from `nvt_vacuum2.gro` must be converted from GROMACS format (nm, GROMACS atom names) to Gaussian format (Å, element symbols), centred on the molecular centre of mass.

**Steps:**
1. Reads `nvt_vacuum2.gro`, skipping the header and box line.
2. Parses into a DataFrame with columns: Residue, atom name, atom number, x, y, z, velocities.
3. Computes the mean x, y, z positions and subtracts them to centre the molecule at the origin.
4. Multiplies coordinates by 10 to convert nm → Å.
5. Maps GROMACS atom names to element symbols.

**V2 atom mapping:** Done via a long `if/elif` chain in `dat_generation.gro_to_dat_atoms()`, covering ~30 hardcoded GROMACS atom type names (e.g. `'CA'` → `'C'`, `'NH'` → `'N'`, `'MW'` → `'Bq'`).

**V3 atom mapping:** Delegated to `Atoms.gaussian_symbol_from_atomname()`, which uses a more robust resolution strategy: dummy prefix matching → NA ambiguity resolution by mass → two-letter element check → one-letter element check → mass-based fallback. This is more generalisable and less brittle for novel atom types.

**Note (V3 WIP):** In `Gaussian_Calculations.gro_to_dat()`, the call `Atoms.Atom_Types(Gros, Masses=0)` passes `Masses=0` rather than the actual mass list, which means the mass-based fallback in `Atoms` cannot function correctly. This is a bug in the current WIP version.

### Vacuum Calculation — `Vacuum.init_v0()` (V2) / `Gaussian_Calculations.init0()` (V3)

Generates and runs the vacuum-phase Gaussian input, then reads back the dipole moment.

**V2 — Two-step job via `--Link1--`:**
- Step 1: Geometry optimisation at B3LYP/6-31+G(d,p) with `fopt=tight`
- Step 2: Geometry optimisation again at B3LYP/6-31+G(d,p) starting from the checkpoint, also with `fopt=tight`

The checkpoint file (`{system_title}_vacuum.chk`) carries the optimised geometry from step 1 into step 2 via `geom=checkpoint`.

**V3 — Single-step job:**
- Single geometry optimisation at B3LYP/6-31+G(d,p) with `fopt=tight`

The `--Link1--` second step has been removed. This is a meaningful method change — the two-step approach in V2 was likely used to ensure convergence, with the second pass refining from the checkpoint. The V3 simplification may reflect confidence that a single pass is sufficient for the molecules studied.

In both versions, the `.dat` file is written to a subdirectory (`opt-test_vacuum/` in V2, `Vacuum/` in V3), a bash script is written to run Gaussian with job throttling, and the script is executed via `check_call`.

**Multipole readback:** After Gaussian completes, `get_multipole_statistics()` (V2) / `get_multipole_statistics(dips_out, rundir)` (V3) reads the output `.out` file using regex patterns to extract the dipole, quadrupole, and octapole moments. Only the total dipole magnitude (`multipole['total dipole']`) is used at this stage. The result is returned as `mu_Vacuum`.

### PCM1 Calculation — `PCM1.init_v0()` (V2) / `Gaussian_Calculations.init1()` (V3)

Runs a two-step Gaussian job using the Polarisable Continuum Model with the **experimental** dielectric constant (`exp_diconst`). This represents the molecule's dipole in the full dielectric environment.

**Both versions — Two-step job via `--Link1--`:**
- Step 1: B3LYP/cc-pVTZ single-point with `scrf=(pcm,solvent={sol_keyword},read)`, reading `eps={exp_diconst}`
- Step 2: B3LYP/aug-cc-pVTZ single-point from checkpoint at the same geometry, same PCM settings

The two-step basis set escalation (cc-pVTZ → aug-cc-pVTZ) allows a cheaper initial SCF to converge before the more expensive augmented basis set single-point is computed. The `geom=checkpoint` keyword in step 2 reuses the geometry from step 1's checkpoint file. The dipole from step 2 (aug-cc-pVTZ) is the value used.

Returns `PCM1` — the liquid-phase dipole moment under full dielectric screening.

### PCM2 Calculation — `PCM2.init_v0()` (V2) / `Gaussian_Calculations.init2()` (V3)

Structurally identical to PCM1 but uses the **optical (calculated) dielectric constant** (`cal_diconst`) instead of the experimental one.

The optical dielectric constant is calculated as:
```
cal_diconst = exp_diconst - ref_ind² + 1
```
This approximation separates the electronic polarisation contribution (captured by the refractive index) from the total dielectric response. PCM2 therefore represents the dipole in an environment that responds only to fast electronic polarisation, not slower orientational/nuclear degrees of freedom.

Returns `PCM2` — the liquid-phase dipole under electronic-only dielectric screening.

### The PCM1/PCM2 Correction Ratio

The ratio `PCM1 / PCM2` is used later in the analysis to correct the SCEE dipole:
```
delta_mu = mu_SCEE * (PCM1 / PCM2) - mu_Vacuum
```
This scales the SCEE result to account for the difference between the dielectric environment experienced by the molecule in the cluster QM calculation versus the true bulk liquid environment.

### Gaussian Job Runner — `run_gaussian()` (all classes, both versions)

All Gaussian classes share the same pattern for running jobs:
1. Write a bash script that loops over `.dat` files in the run directory.
2. For each file, check the combined count of running `g09` and `mdrun` processes.
3. If the count is below `MAXJOBS`, launch the Gaussian job in the background with `g09 < datfile > outfile &`.
4. After all jobs are launched, `wait` for completion.

**V3 enhancement:** The bash script supports multiple scratch directories as a round-robin array (`SCRATCH_DIRS`), cycling through them across jobs. This was a single fixed path in V2.

---

## Outputs

| File | Description |
|------|-------------|
| `em.gro` | Energy-minimised solute in vacuum |
| `nvt_vacuum2.gro` / `nvt_vacuum2.edr` / `nvt_vacuum2.trr` | Vacuum NVT trajectory and outputs |
| `Dipole.txt` | Raw GROMACS dipole output for the vacuum model |
| `Vacuum/Vacuum.dat` (V3) / `opt-test_vacuum/nvt_vacuum.dat` (V2) | Gaussian input for vacuum calculation |
| `Vacuum/Vacuum.out` (V3) / `opt-test_vacuum/nvt_vacuum.out` (V2) | Gaussian output for vacuum calculation |
| `Dipole_Vacuum.csv` | Vacuum dipole moments (`dipole_l`, `dipole_m`, `dipole_h`, `muG`) |
| `PCM1/PCM1.dat` (V3) / `opt_b3lypaug-cc-pvtz.../PCM1.dat` (V2) | Gaussian input for PCM1 calculation |
| `PCM1/PCM1.out` (V3) | Gaussian output for PCM1 calculation |
| `Dipole_PCM1.csv` | PCM1 dipole moments |
| `PCM2/PCM2.dat` (V3) | Gaussian input for PCM2 calculation |
| `PCM2/PCM2.out` (V3) | Gaussian output for PCM2 calculation |
| `Dipole_PCM2.csv` | PCM2 dipole moments |
| `run_.sh` | Bash script used to execute Gaussian jobs |

---

## Changes from V2 to V3 and Why

| Change | Reason |
|--------|--------|
| `dat_generation.py` absorbed into `Gaussian_Calculations.py` | Consolidation — the coordinate conversion was only ever used as a pre-step to Gaussian input generation, so keeping it as a separate class added unnecessary indirection. |
| `Vacuum.py`, `PCM1.py`, `PCM2.py` merged into `Gaussian_Calculations.py` as `init0`, `init1`, `init2` | All three classes were structurally near-identical (same `run_gaussian`, `read_multipoles`, `get_multipole_statistics` pattern). Merging eliminates significant code duplication. |
| `Atoms.py` introduced for atom name resolution | Replaces the brittle hardcoded `if/elif` chains in V2. More generalisable to new molecules. |
| `_wait_for_md_slot()` removed from `Gro_Simulations` | Job throttling for MD is less critical than for Gaussian — the Gaussian scripts throttle themselves. |
| Vacuum calculation simplified from two-step to single-step | Removes the redundant second `fopt` pass from `--Link1--`. |
| Run directories renamed (`opt-test_vacuum/` → `Vacuum/`, `opt_b3lypaug-cc-pvtz_spb3lyp/` → `PCM1/`, `PCM2/`) | Cleaner, more descriptive naming. |
| Hardware config moved to `Settings.yml` | V2 used hardcoded hostname checks and absolute paths. V3 reads these from YAML, making the script portable across machines. |
| `Pre_Eq_Simulations.Pre_Eq_Solute` absorbed into `Gro_Simulations.run_md1` | The vacuum equilibration and molecule insertion are now a single logical unit. |

---

## Other Notes

- **`get_multipole_statistics` bug in V3 (WIP):** The method loops `Q` from 1 to 3 (intended to separately read the `_q1`, `_q2`, `_q3` output files) but the glob pattern inside does not use `Q` — it always uses the same `dips_out` string. For `init0/1/2` this doesn't matter because there is only one output file each (Vacuum/PCM1/PCM2 are single jobs), but the three-column DataFrame structure (`dipole_l`, `dipole_m`, `dipole_h`) will contain identical values rather than distinct results. This is a carry-over from the SCEE `get_multipole_statistics` structure which does differentiate by Q.
- **Variable name bugs in `init1` and `init2` (V3 WIP):** `init1` calls `get_multipole_statistics` into `df4` but then references `df2`; `init2` similarly references `df3`. These are unfixed carry-over names from the V2 split-class structure.
- **`HO`, `HN`, `HW` Lennard-Jones override:** In both versions, hydrogen atoms bonded to heteroatoms (`HO`, `HN`, `HW`) have their σ and ε hardcoded to 0.2673 Å and 0.0418 kJ/mol respectively, rather than reading from the force field. This is because OPLS-AA assigns zero LJ parameters to these hydrogens, which would cause issues in the ONIOM calculation. The hardcoded values provide a small but non-zero interaction.
- The `sol_keyword` used for Gaussian's PCM must be a keyword that Gaussian recognises. For Aniline, `'Aniline'` is a valid Gaussian solvent keyword. For molecules without an exact match, an approximate keyword must be used — for example, benzyl alcohol was used as a substitute for phenol in earlier runs.
