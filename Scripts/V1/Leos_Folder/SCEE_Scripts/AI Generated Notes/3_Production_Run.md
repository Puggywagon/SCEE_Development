# Section 3: The Production Run

---

## Collective Goal

This section takes the equilibrated single-molecule geometry from Section 2, builds and equilibrates a full solvent box at the target density, runs the production MD simulation, and then processes the trajectory into per-configuration ONIOM inputs for Gaussian. It also covers the ONIOM input configuration file generation and the SCEE Gaussian calculations themselves (both the ONIOM geometry optimisation step and the single-point charge embedding step). The outputs of this section are per-configuration dipole moments at three charge scaling levels (`dipole_l`, `dipole_m`, `dipole_h`) that are passed to Section 4 for analysis.

---

## Script Files Involved

| Version | Files |
|---------|-------|
| V2 | `Pre_Eq_Simulations.py`, `MD_Simulations.py`, `Oniom_Generation.py`, `SCEE.py`, `Shell_Oniom.f90` |
| V3 | `Gro_Simulations.py`, `Oniom_Generation.py`, `Gaussian_Calculations.py`, `Atoms.py`, `Shell_Oniom.f90` |

---

## Imports and Inputs

### Imports

**V2:**
```python
# Pre_Eq_Simulations.py / MD_Simulations.py
import gromacs, subprocess, numpy, os, socket, time

# Oniom_Generation.py
import pandas, re

# SCEE.py
import glob, gromacs, pandas, re, os, numpy
from subprocess import check_call
```

**V3:**
```python
# Gro_Simulations.py
import gromacs, subprocess, numpy, os, socket, time

# Oniom_Generation.py
import pandas, re
import Atoms

# Gaussian_Calculations.py
import glob, gromacs, pandas, re, os, numpy
import Atoms
import Simulations_Analysis
```

### User Inputs

**V2** — inline variables in `Run_SCEE.py`:
| Variable | Description | Example (Aniline) |
|----------|-------------|-------------------|
| `density` | Experimental density (g/cm³) | `1.022` |
| `mol_mass` | Molar mass (g/mol) | `93.13` |
| `L` | Target box side length (nm) | `6.2` |
| `initial_molecules` | Calculated number of solvent molecules | computed |
| `solresnametop` | Residue name as it appears in the topology `[molecules]` section | `'Ani'` |
| `Configurations` | Number of trajectory frames to extract | `200` |
| `Cut_Off_Radius` | Cluster radius for ONIOM shell selection (nm) | `2.8` |
| `qr1`, `qr2`, `qr3` | Charge scaling ratios | `1.0`, `1.35`, `1.7` |
| `replicas` | Range of replica indices | `range(3, 6, 1)` |
| `T_list`, `P_list` | Temperature and pressure lists | `['298.0']`, `['1.0']` |

**V3** — from `Settings.yml`:
| YAML Key | Description |
|----------|-------------|
| `State_Conditions.Temperature` | List of temperatures |
| `State_Conditions.Pressure` | List of pressures |
| `State_Conditions.Replicas` | Number of replicas |
| `Advanced_Settings.configuration.box_length_nm` | Box side length |
| `Advanced_Settings.configuration.cluster_radius_nm` | ONIOM cluster cutoff radius |
| `Advanced_Settings.sampling.n_configurations` | Number of frames to extract |
| `Advanced_Settings.electrostatics.charge_scaling.qr1/qr2/qr3` | Charge scaling ratios |
| `No_Gro.density` / `Build_Gro.density` | Experimental density |
| `No_Gro.mol_mass` / `Build_Gro.mol_mass` | Molar mass |
| `No_Gro.solresnametop` / `Build_Gro.solresnametop` | Residue name in topology |

### Input Files (Template Folder)
| File | Description |
|------|-------------|
| `minim.mdp` | GROMACS MDP for box energy minimisation |
| `nvt.mdp` | GROMACS MDP for box NVT equilibration |
| `npt.mdp` | GROMACS MDP for box NPT equilibration |
| `template_2.mdp` | Template MDP for production run (c-rescale barostat); `TEMPERATURE` and `PRESSURE` are placeholders replaced at runtime. Note: `template_1.mdp` uses a Parrinello-Rahman barostat — used for molecules requiring a different barostat for stability |
| `{system_title}.top` (V2) / `{Topology_File}` (V3) | Full system topology file |
| `oplsaaff.itp` (V2) | OPLS-AA force field file (required when `split='yes'` and `dummy='no'`) |
| `{central}_AA.itp` (V2) | Solute topology include file |
| `{solvent}_UA.itp` (V2) | Solvent topology include file |
| `nvt_vacuum2.gro` | Equilibrated single molecule from Section 2 |
| `oniom.inp` | ONIOM configuration file (generated in this section) |
| `Shell_Oniom.f90` | Fortran source for trajectory processing |

---

## How This Section is Called in `Run_SCEE`

**V2:**
```python
# Box building
Pre_Eq_Simulations.insert_Molecules(central, solvent, system_title, initial_molecules)
Pre_Eq_Simulations.Write_to_top(system_title, initial_molecules, solresnametop)
Pre_Eq_Simulations.Pre_Eq_System(system_title)

# ONIOM input generation
Oniom_Generation.Gen_File(Oniom, Configurations, system_title, Cut_Off_Radius)
# (QM and MM inputs via split itp logic)
Oniom_Generation.QM_Inputs(Solute, Oniom, dummy, split, qr1, qr2, qr3, qmax)
Oniom_Generation.MM_Inputs(Solvent, Oniom, split, qr1, qr2, qr3, qmax)
Oniom_Generation.Counting_Molecules(Solvent, Oniom, initial_molecules)

# Production MD (inside replica/T/p loop)
MD_Simulations.create_mdpfile(HOMEDIR, 'junk.mdp', T, p)
MD_Simulations.run_md('junk.mdp', HOMEDIR, system_title)
MD_Simulations.process_trajectory(system_title)
MD_Simulations.process_gro()

# SCEE Gaussian calculations (inside replica loop)
Molecule.process_gro(exe=HOMEDIR+'/shellO', inp=HOMEDIR+'/oniom.inp')
Molecule.init_v0()
Molecule.run_gaussian(step=0)
Molecule.init_v1()
Molecule.run_gaussian(step=1)
```
In the V2 Aniline run, all of the above calls are commented out and the script reads back results from completed runs.

**V3:**
```python
# Box building (inside Box_Build == 'Yes' block, already includes vacuum prep from Section 2)
# insert_Molecules and Write_to_top are called at the end of run_md1

# Box equilibration (always called regardless of mode)
MD.run_md2(Topology_File, Mixture='No', MD='Box')

# ONIOM input generation
Topologys = Oniom_Generation.expand_topology_with_itps(Topology_File)
# ... written to Concat_Top.top ...
Oniom_Generation.Gen_File(Oniom, Configurations, solvent, Cut_Off_Radius)
Total_Atoms, qmax = Oniom_Generation.QM_Inputs(Topology, Oniom, qr1, qr2, qr3)
Oniom_Generation.MM_Inputs(Topology, Oniom, qr1, qr2, qr3)
Oniom_Generation.Counting_Molecules(Oniom, initial_molecules)

# Production MD + SCEE (inside replica/T/p loop)
MD.run_md3(mdpfile, HOMEDIR, system_title, T, p, Mixture='No', MD='Production')

# SCEE Gaussian calculations (inside replica directory loop)
Gauss.process_gro(exe=HOMEDIR+'/shellO', inp=HOMEDIR+'/oniom.inp')
Gauss.init3(Gaus='SCEE_V0')
SCEE = Gauss.init4(Gaus='SCEE_V1')
```
**Note (V3 WIP):** `system_title` is not currently defined anywhere in V3 `Run_SCEE.py` — this is an open bug. `sol_keyword` and `gas_dipole` are also referenced but not yet read from the YAML.

---

## How the Relevant Functions Work

### Molecule Count Calculation (both versions)

Before building the box, the number of solvent molecules required to fill a cubic box of side `L` nm at the experimental density is calculated:

```python
N = (density * (L * 1e-7)**3 / mol_mass) * 6.022e23
initial_molecules = int(ceil(N + 0.05*N)) - 1
```

The 5% excess accounts for the fact that `gmx insert-molecules` may not always successfully insert all requested molecules due to steric clashes. One molecule is subtracted because the box already contains the solute molecule from `nvt_vacuum2.gro`.

### Molecule Insertion — `Pre_Eq_Simulations.insert_Molecules()` (V2) / `Gro_Simulations.insert_Molecules()` (V3)

Uses `gromacs.insert_molecules` to pack `initial_molecules` copies of the solvent geometry (`{solvent}_UA.gro` in V2, `nvt_vacuum2.gro` in V3) around the equilibrated solute in `nvt_vacuum2.gro`, producing `out.gro`.

**V2:** Inserts copies of the solvent UA `.gro` file — keeping solute and solvent representations distinct.

**V3:** Inserts copies of `nvt_vacuum2.gro` into itself — this reflects the pure liquid case where solute and solvent are the same molecule. The `Mixture` parameter is passed through but not yet used internally.

### Topology Update — `Write_to_top()` (both versions)

Appends the solvent molecule count to the `[molecules]` section of the topology file:
```
{solresnametop}    {initial_molecules}
```
This tells GROMACS how many solvent molecules are in the system alongside the single solute molecule already listed in the topology.

### Box Equilibration — `Pre_Eq_Simulations.Pre_Eq_System()` (V2) / `Gro_Simulations.run_md2()` (V3)

Runs a three-stage equilibration of the full solvent box:
1. **Energy minimisation** (`minim.mdp`): removes bad contacts from the insertion step → `em1.gro`
2. **NVT equilibration** (`nvt.mdp`): brings the system to the target temperature → `nvt_eq.gro`
3. **NPT equilibration** (`npt.mdp`): equilibrates volume/density at the target pressure → `{system_title}.gro` (V2) / `system.gro` (V3)

The output of the NPT step is the starting point for all production runs.

**V3 note:** `run_md2` is called regardless of `Box_Build` mode, even when the user provides an existing `.gro` file. This is intentional — as noted in the script comment: "Just because they gave us a box I don't trust that they have minimised it well."

### Production MD — `MD_Simulations.run_md()` (V2) / `Gro_Simulations.run_md3()` (V3)

Runs the production NpT MD simulation for each replica/temperature/pressure combination:
1. Creates a customised MDP file from `template_2.mdp` by replacing the `TEMPERATURE` and `PRESSURE` placeholders.
2. Runs `gromacs.grompp` + `gromacs.mdrun` to produce the trajectory (`.xtc`), energy (`.edr`), and coordinate (`.gro`) files, all named `{system_title}_QMMM_md3.*`.
3. Calls `process_trajectory` to extract frames.
4. Calls `process_gro` to run Shell_Oniom.

The run operates from within the replica/T/p directory, with `HOMEDIR` used to locate the topology and starting structure.

**Template MDP note:** `template_2.mdp` uses the c-rescale barostat. `template_1.mdp` (Parrinello-Rahman) is available for molecules that are unstable with c-rescale. The choice of which template to use is currently a manual decision.

### Trajectory Processing — `process_trajectory()` (both versions)

After the production run, `gromacs.tools.Trjconv` is used to extract individual frames from the trajectory:
- Source: `{system_title}_QMMM_md3.xtc`
- Parameters: `pbc=whole` (makes molecules whole across PBC), starting from 1000 ps (`b=1000`) with 20 ps spacing (`dt=20`)
- Output: `conf_0.gro`, `conf_1.gro`, ... `conf_N.gro` — one file per extracted frame

The `b=1000` parameter skips the first 1000 ps to allow for equilibration of the production trajectory. With `dt=20` and a typical run of ~5 ns, this yields approximately 200 configurations.

### Shell_Oniom Processing — `process_gro()` + `Shell_Oniom.f90`

This is the critical step that converts the extracted `.gro` frames into ONIOM-formatted Gaussian input files.

`process_gro()` compiles and runs the Fortran program:
```bash
gfortran -o Shell_Oniom Shell_Oniom.f90 && ./Shell_Oniom
```
In V2 this is a standalone call in `MD_Simulations`. In V3 it is embedded at the end of `run_md3`. In the SCEE calculation loop, `Gauss.process_gro()` in V3 additionally copies `oniom.inp`, the compiled `shellO` executable, and `Shell_Oniom.f90` into the replica directory before running.

**Shell_Oniom.f90 — detailed operation:**

For each of the `N_conf` configurations (read from `oniom.inp`):

1. Opens the `oniom.inp` file and reads the control parameters:
   - Number of configurations
   - Output file name stem (the solvent name)
   - Cutoff radius for the solvent shell (in nm)
   - Solute block: atom count, central atom index, then per-atom: element symbol, GROMACS type, σ, ε, three charge values (q1/q2/q3), dummy flag
   - Solvent block: atom count, central atom index, then per-atom: element symbol, GROMACS type, σ, ε, three charge values
   - Total solvent molecule count

2. Opens `conf_N.gro` and reads all atomic coordinates plus box dimensions.

3. **Unit and unit-cell conversion:**
   - Translates all coordinates so the solute central atom is at the origin.
   - Applies minimum image convention (PBC) to handle atoms that cross box boundaries.
   - Converts nm → Å (multiplies by 10).
   - Converts σ from GROMACS units using factor 0.561231 (Amber distance convention).
   - Converts ε from kJ/mol to kcal/mol (divides by 4.184).

4. **Solvent shell selection:** For each solvent molecule, the distance from its central atom to the solute central atom is computed (with PBC). If the distance ≤ `Rcut1`, the molecule is flagged for inclusion.

5. **Even/odd correction:** If the number of selected solvent molecules (plus solute) is even, the furthest selected molecule is removed to make the count odd. This corrects for a known Gaussian bug with even numbers of atoms in certain ONIOM calculations.

6. **Outputs six files per configuration:**
   - `{solvent}_cN_q1.inp`, `_q2.inp`, `_q3.inp` — full ONIOM geometry files at three charge scaling levels. Solute atoms are tagged `H` (high level, QM), solvent atoms tagged `L` (low level, MM). Dummy atoms are excluded from coordinates and connectivity.
   - `{solvent}_cN_q1_chg.inp`, `_q2_chg.inp`, `_q3_chg.inp` — charge-only files for the solvent shell, containing only xyz coordinates and the per-atom charge for each scaling level. Used as point charge embeddings in the V1 single-point calculation.

7. **Writes connectivity section** for the ONIOM input (one line per atom listing its index), then the VDW parameter section listing `VDW {type} {sigma} {epsilon}` for all solute and solvent atom types.

### ONIOM Input Configuration — `Oniom_Generation` (both versions)

Before `Shell_Oniom.f90` runs, `oniom.inp` must be built by `Oniom_Generation`. This file tells the Fortran code how to interpret the topology and what parameters to use.

**`Gen_File()`:** Writes the header block:
```
{N_configurations}
{solvent_name}        ← used as the output file stem
{cutoff_radius}
```
**`QM_Inputs()`:** Writes the solute (QM region) block. Reads atom names, masses, charges, and LJ parameters from the topology. The three scaled charge sets are computed as:
- For the atom with maximum absolute charge (`wi = 1`): `q = qmax * qr`
- For other atoms: `q = (qi/qmax) * qmax * qr = qi * qr` — i.e. all charges are scaled by the same ratio
- Special handling for water: OW gets negative charges, MW (dummy) gets zero charges

The central atom (`Central_Atom`) is identified as the atom with the largest absolute charge (or the OW atom index for water). This is the atom that Shell_Oniom uses as the reference point for distance calculations and coordinate translation.

**V2 vs V3 key difference in `QM_Inputs`:**
- V2: Takes separate `Solute` (from `{central}_AA.itp`) and looks up LJ parameters via the `oplsaaff.itp` force field file. The `dummy` and `split` flags control which columns the topology is expected to have.
- V3: Takes the full concatenated topology string (from `expand_topology_with_itps`) and reads both charges (from `[atoms]`) and LJ parameters (from `[atomtypes]`) from the same source. No separate force field file lookup needed. `Atoms.Atom_Types()` handles the element symbol resolution.

**`MM_Inputs()`:** Writes the solvent (MM region) block using the same logic as `QM_Inputs` but reading the solvent molecule's atoms. In V3 this reads the second `[atoms]` block from the concatenated topology. Note that in V3 the charge sign convention for water is **reversed** relative to `QM_Inputs` — MM OW gets zero charges and MM MW gets the negative scaled charges, consistent with the TIP4P model where the charge sits on the dummy atom.

**`expand_topology_with_itps()` (V3 only):** Recursively resolves all `#include` statements in the topology file, building a single concatenated text string. Special handling: `forcefield.itp` itself is not included, but `ffnonbonded.itp` within it is — this is where the `[atomtypes]` σ and ε values live. A `seen` set prevents circular includes.

**`Counting_Molecules()`:** Appends the total solvent molecule count to `oniom.inp`.

### `get_multipole_statistics()` and `qmax` Usage

`qmax` (the maximum absolute partial charge in the QM molecule) and `Total_Atoms` are returned from `Oniom_Generation.QM_Inputs()` and stored in `Run_SCEE`. They are used to:
- Set `Gauss.natom = Total_Atoms` so the SCEE single-point knows how many atoms to extract from the Gaussian Z-matrix output.
- Compute the charge ratio: `ratio = Model_Dipole_liquid / qmax`, from which the three actual charge values passed to Gaussian are derived: `num1 = qr1 * ratio * qmax`, etc.

### SCEE Calculations — `SCEE.init_v0/v1()` (V2) / `Gaussian_Calculations.init3/init4()` (V3)

**Step 1 — ONIOM geometry optimisation (`init_v0` / `init3`):**

For each of the ~200 `_c*_q?.inp` files (three per configuration, one per charge scaling level), a Gaussian ONIOM input is assembled:
- Route card: `oniom({method_v0}/{basis_v0}:amber=(print,SoftFirst,LastEquiv))=embedcharge`
- Options: `opt=quadmacro iop(1/98=66,1/19=7)`
- Charge/spin: `0 1 0 1 0 1` (neutral, singlet for both QM and MM regions)

The `_c*_q?.inp` geometry from Shell_Oniom is appended to this header to form the complete `.dat` file. All files are written to the `SCEE/` subdirectory (V3) or `opt-test-SCEE/` (commented out, V2). Gaussian runs are submitted via the bash throttling script.

**Step 2 — Single-point charge embedding (`init_v1` / `init4`):**

After the ONIOM optimisation, the relaxed QM-region geometry is extracted from each `.out` file using a regex on the Z-Matrix orientation block. For each configuration, a new single-point calculation is set up:
- Route card: `{method_v1}/{basis_v1} charge` — runs at the higher basis set level with the solvent point charges embedded
- The QM atom coordinates come from the Z-matrix of the V0 output
- The solvent charges come from the corresponding `_chg.inp` file

`Atoms_Dict.Atom_Symbol()` (V3) / `self.atom_dict[atomic_number]` (V2) converts the atomic number from the Z-matrix output back to an element symbol.

**In V3, charge scaling is computed inside `init4`:**
```python
Model_Dipole = Analysis.get_dipole_model_liquid(system_title)
ratio = Model_Dipole / qmax
num1 = qr1 * ratio * qmax
num2 = qr2 * ratio * qmax
num3 = qr3 * ratio * qmax
```
This is a significant structural change — in V2 this calculation lived in `Run_SCEE.py`. Moving it into `init4` creates a dependency on `Analysis.get_dipole_model_liquid()` being callable from within `Gaussian_Calculations`, which requires the production trajectory to already exist.

`get_multipole_statistics()` then reads back all `_v1.out` files, extracts dipole moments for each configuration at each charge level, and organises them into a DataFrame with columns `config`, `dipole_l`, `dipole_m`, `dipole_h`. Configurations where any of the three calculations failed (NaN values) are dropped. The resulting `muL_SCEE` per configuration is computed via the polynomial self-consistency fitting described in Section 4.

---

## Outputs

| File | Description |
|------|-------------|
| `out.gro` | Initial box after molecule insertion |
| `em1.gro` | Box after energy minimisation |
| `nvt_eq.gro` | Box after NVT equilibration |
| `{system_title}.gro` (V2) / `system.gro` (V3) | Box after NPT equilibration — production starting point |
| `{system_title}_QMMM_md3.xtc/.edr/.gro/.trr/.tpr` | Production MD trajectory and associated files |
| `conf_0.gro` ... `conf_N.gro` | Extracted trajectory frames |
| `oniom.inp` | ONIOM configuration file for Shell_Oniom |
| `Concat_Top.top` (V3 only) | Concatenated topology with all includes resolved |
| `shellO` (V3) | Compiled Shell_Oniom executable |
| `{solvent}_cN_q1/q2/q3.inp` | ONIOM geometry inputs for each config at three charge levels |
| `{solvent}_cN_q1/q2/q3_chg.inp` | Point charge files for each config at three charge levels |
| `SCEE/{solvent}_cN_q?.dat` | Gaussian ONIOM input files (V0) |
| `SCEE/{solvent}_cN_q?_v1.dat` | Gaussian single-point input files (V1) |
| `SCEE/{solvent}_cN_q?.out` | Gaussian ONIOM outputs (V0) |
| `SCEE/{solvent}_cN_q?_v1.out` | Gaussian single-point outputs (V1) |
| `Dipole_SCEE.csv` | Per-configuration `dipole_l`, `dipole_m`, `dipole_h`, `muL_SCEE` |

---

## Changes from V2 to V3 and Why

| Change | Reason |
|--------|--------|
| `Pre_Eq_Simulations`, `MD_Simulations`, `SCEE` merged into `Gro_Simulations` and `Gaussian_Calculations` | Eliminates near-duplicate class structures. MD steps are now numbered `run_md1/2/3` reflecting their natural sequence. |
| `run_md2` always runs even when box is provided | Defensive practice — user-provided boxes may not be well-equilibrated. |
| ONIOM generation no longer requires separate `.itp` files | `expand_topology_with_itps` resolves the full topology in memory, removing the dependency on `split`, `dummy`, and `oplsaaff.itp` flags. |
| `Atoms.py` replaces hardcoded atom name maps in `Oniom_Generation` | More robust resolution, especially for two-letter elements (Cl, Na, etc.) noted as a known issue in V2 code comments. |
| `atomtypes` parser handles both 7- and 8-column formats | V2 assumed a fixed column structure. V3 detects which format by checking if column 4 is `'A'` or `'D'`. |
| `Gen_File` second argument changed from `system_title` to `solvent` | The output file stem in `oniom.inp` should be the solvent name (used by Shell_Oniom as the file prefix), not the full system title. |
| Charge scaling moved into `init4` from `Run_SCEE` | Keeps the charge computation closer to where it's used, though creates an internal dependency on `get_dipole_model_liquid` being available. |
| Production run output renamed from `{system_title}.gro` to `system.gro` | Consistent naming independent of molecule name. |

---

## Other Notes

- **Template MDP files:** Two templates exist — `template_1.mdp` (Parrinello-Rahman barostat) and `template_2.mdp` (c-rescale barostat). The script always uses `template_2.mdp`. Switching between them requires a manual code change. The choice matters for stability: some molecules run stably only with one or the other. For the Aniline run, `template_2.mdp` was used.
- **Phenol case (Other Notes):** During the Phenol run, the NPT pre-equilibration step from `Pre_Eq_System` / `run_md2` caused the production run to fail. The workaround was to merge the NPT step with the production run into a single extended simulation. The reason for this behaviour is not fully understood. This is a specific case where the standard pipeline was modified.
- **Even/odd molecule correction in Shell_Oniom:** The removal of the furthest solvent molecule when the cluster count is even is a workaround for a Gaussian bug. The exact nature of this bug (likely related to spin state or orbital counting in ONIOM with AMBER embedding) is not documented in the code.
- **`Mixture` parameter:** Both `Gro_Simulations` methods and `Shell_Oniom` have infrastructure for mixture support (e.g., different solute/solvent molecules) but this is not yet active. The parameter is accepted but not used internally.
- **`process_gro` signature mismatch (V3 WIP):** `Gro_Simulations.run_md3` calls `self.process_gro(Mixture)` but `process_gro` takes no arguments. This is a bug in the current WIP version.
- **`init4` variable scope issue (V3 WIP):** `system_title`, `qmax`, `qr1`, `qr2`, `qr3` are referenced inside `Gaussian_Calculations.init4()` but are not passed as arguments or defined as class attributes — they are expected to be in the calling scope. This will fail when called as a class method.
