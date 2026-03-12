# Section 1: Building the Molecule

---

## Collective Goal

This section covers the generation of initial molecular geometry files (`.gro` format) for the solute and solvent molecules directly from their SMILES string representations. The aim is to produce a starting structure that is chemically reasonable and formatted correctly for use in subsequent GROMACS simulations. In V2 this also included the generation of united-atom (UA) representations; in V3 this has been simplified to all-atom (AA) structures only.

---

## Script Files Involved

| Version | File |
|---------|------|
| V2 | `Gro_Builder.py` |
| V3 | `Gro_Builder.py` |

---

## Imports and Inputs

### Imports
```python
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
```

### User Inputs (from `Run_SCEE.py`)

**V2** — provided inline as hardcoded variables:
| Variable | Description | Example (Aniline run) |
|----------|-------------|----------------------|
| `mol` | SMILES string of the solute/central molecule | `'c1ccccc1N'` |
| `resname` | Residue name for the solute in the `.gro` file | `'ANI'` |
| `system_title` | Name used for the output file | `'Aniline'` |
| `solmol` | SMILES string of the solvent molecule | `'c1ccccc1N'` (same as solute for a pure liquid) |
| `solresname` | Residue name for the solvent in the `.gro` file | `'ANI'` |

**V3** — provided via `Settings.yml` under the `Build_Gro` mode block:
| YAML Key | Description |
|----------|-------------|
| `solmol` | SMILES string |
| `solvent` | Name used for output file (equivalent to `system_title`) |
| `solresname` | Residue name in the `.gro` file |

### Input Files
None. The structure is generated entirely from the SMILES string using RDKit.

---

## How This Section is Called in `Run_SCEE`

**V2:**
```python
Gro_Builder.AA_Structure(mol, system_title, resname)
Gro_Builder.UA_Structure(solmol, system_title, solresname)
```
Both calls are made unconditionally at the start of the run, regardless of whether `.gro` files already exist.

**V3:**
`Gro_Builder.AA_Structure` is only called if the user has selected `Mode: Build_Gro` in `Settings.yml`. The returned filename is stored and passed forward:
```python
Gro_File = Gro_Builder.AA_Structure(solmol, solvent, solresname)
```
After generation, the script pauses (in interactive use) to allow the user to visually confirm the structure is correct before continuing. In the current WIP version this confirmation is hardcoded to `'Yes'`.

---

## How the Relevant Functions Work

### `Gro_Builder.AA_Structure(mol, molecule, resname)`

1. Takes the SMILES string and creates an RDKit molecule object via `Chem.MolFromSmiles()`.
2. Adds explicit hydrogen atoms with `Chem.AddHs()`.
3. Generates an initial 3D conformation using `AllChem.EmbedMolecule()`, which uses distance geometry.
4. Optimises the geometry using the Universal Force Field (UFF) via `AllChem.UFFOptimizeMolecule()`.
5. Generates a 2D depiction image of the molecule using `Chem.Draw.MolToImage()`. In V2 this image is saved to disk as `{system_title}.png`; in V3 the save is commented out (noted as not working in the developer's environment — likely a display/backend issue).
6. Writes a GROMACS `.gro` file with the following format:
   - Header line: molecule name
   - Atom count line
   - One line per atom: residue number, residue name, element symbol, atom index, x, y, z coordinates (in nm — RDKit positions are in Å and are divided by 10)
   - Box size line: all zeros (no box at this stage)

**V3 change:** `AA_Structure` now returns `Gro_File` (the filename string), which is used downstream. V2 did not return anything.

### `Gro_Builder.UA_Structure(solmol, system_title, solresname)` *(V2 only)*

This method existed in V2 to generate a united-atom representation of the solvent molecule, where non-polar C–H groups are represented implicitly. It used the same RDKit embedding and UFF optimisation pipeline as `AA_Structure`, with a `skips` counter intended to track omitted hydrogen atoms — though in the V2 code this counter was always zero and never incremented, meaning the UA output was functionally identical to the AA output. The method also called `Functional_Group()` to classify the molecule, though this classification was not actually used to modify the output.

### `Gro_Builder.Functional_Group(solmol)` *(V2 only)*

Classified the molecule by functional group using SMARTS pattern matching via RDKit. The priority order of detection was: Ester > Aldehyde > Amine (secondary/tertiary/primary) > Thiol > Ether > Alcohol (secondary/tertiary/primary) > Diol > Amino Acid > Phenol > Aniline > Benzene > Carbons only. This function was called within `UA_Structure` but its return value was never used to alter the `.gro` output. It appears to have been intended for future logic to properly handle UA representations of different functional groups.

---

## Outputs

| File | Description |
|------|-------------|
| `{molecule}.gro` (V3) / `{system_title}_AA.gro` (V2) | All-atom GROMACS coordinate file for the solute |
| `{system_title}_UA.gro` (V2 only) | United-atom GROMACS coordinate file for the solvent |
| `{system_title}.png` (V2 only) | 2D depiction image of the molecule |

---

## Changes from V2 to V3 and Why

| Change | Reason |
|--------|--------|
| `UA_Structure` and `Functional_Group` methods removed | The UA method was not functionally producing a true united-atom structure (the `skips` counter was never used), making it effectively redundant. V3 simplifies to AA only. |
| `AA_Structure` now returns the filename | Allows `Run_SCEE` to capture and pass the `.gro` filename to subsequent MD steps without hardcoding it. |
| Structure generation is now mode-gated (`Build_Gro`) | V2 always generated `.gro` files regardless of need. V3 allows users to provide pre-existing files, so generation is only triggered when explicitly requested. |
| Output filename changed from `{system_title}_AA.gro` to `{molecule}.gro` | Cleaner naming when solute and solvent are the same molecule (pure liquid case). |
| Image saving commented out | Known environment compatibility issue with matplotlib display backend. |

---

## Other Notes

- The `.gro` files generated here have a box size of zero (all box dimensions set to `0.0`). A proper simulation box is assigned later in `Gro_Simulations.run_md1()` using `gromacs.editconf`. This means these files cannot be used directly with GROMACS without that box assignment step.
- The topology file (`.top`) for the molecule is **not** generated by this script — it must be provided by the user separately. Maintaining consistency between the RDKit-generated `.gro` geometry and the topology file is the user's responsibility, as noted in the `Run_SCEE` user-facing print statements.
- RDKit's UFF optimisation gives a reasonable starting geometry but is not a high-level optimisation. The geometry will be refined further during the vacuum phase energy minimisation steps.
- For the Aniline run, the solute and solvent SMILES (`'c1ccccc1N'`) are identical since it is a pure liquid simulation.
