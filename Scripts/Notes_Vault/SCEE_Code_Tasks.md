# SCEE Pipeline — Outstanding Work

Companion to [[SCEE_Code_Overview]], alongside [[SCEE_Code_Testing]]. Scoped to the **code**, not the wider project.

Everything here was derived from reading the source as it currently stands, so some items may already be
fixed in the working copy — treat unticked boxes as "worth confirming" rather than "definitely broken".

Grouped by what the work actually is: things that could affect a number, things that finish an
unfinished feature, things that make the code survivable by someone else, and things that need a
decision rather than a commit.

---

## 1. Correctness — could change a result

- [ ] **`dipole_stats.csv` filename casing.** `Run_SCEE.py` deletes `Results/Dipole_Stats.csv`; the analysis
      writes `Results/dipole_stats.csv`. On Linux these are different files, so filtering statistics
      accumulate across runs rather than resetting. Affects the audit trail, not the dipoles themselves.
- [ ] **`natom` includes dummy sites.** `Gauss.natom = Total_Atoms` counts dummies, but the ONIOM coordinate
      block excludes them, and `init4` reads `natom` rows of the Z-matrix orientation. For any model with
      virtual sites this reads past the QM molecule into the bath. Blocking for the dummy-hydrogen UA work.
- [ ] **Frame numbering off-by-one.** `trjconv -sep` numbers from `conf_0.gro`; `Shell_Oniom` reads
      `conf_1` upward. Confirm the production run yields `n_configurations` usable frames after the skip.
- [ ] **Insertion shortfall is never checked.** `gmx insert-molecules` can place fewer molecules than
      requested; both the topology line and `oniom.inp` are written from the *requested* count. Currently
      surfaces as a downstream `grompp` or `Shell_Oniom` exception — an explicit check on the output `.gro`
      would fail faster and more clearly.
- [ ] **Even-count workaround.** Confirm the Gaussian bug that motivates dropping the furthest molecule is
      still present in the version being used; if so, note the version in the code comment.

## 2. Unfinished features

- [ ] **Mixture branch.** Currently builds the solute and runs MD1/MD2/MD3, then stops — no ONIOM
      generation, Gaussian calculation or analysis. Decide scope before building (see §4).
- [ ] **`g09root` auto-detection.** Raises `NotImplementedError`; the path is machine-specific in the YAML.
- [ ] **Dielectric constant** is a post-hoc correction rather than part of the iteration.
- [ ] **Optimum model dipole** only partially integrated — first calculation is in, the loop is not.
- [ ] **Restore the interactive confirmation prompts** behind a `--batch` style flag; they are currently
      hard-coded to `'Yes'`, which is right for production and wrong for a first run on a new molecule.
- [ ] **Create `Results/` if missing** — nothing does, and the run fails late.

## 3. Robustness and handover

- [ ] **Topology ↔ RDKit atom ordering** is the one genuinely manual step and the most likely source of
      silent error for a new user. Either automate (SMILES-based topology generation) or add a validation
      check that compares `.itp` atom order against the `.gro` before anything runs.
- [ ] **Hard-coded polar-hydrogen LJ floor** (σ = 0.2673 Å, ε = 0.0418 kJ/mol) applied when `HO`/`HN`/`HW`
      carry zero parameters. Needs a source or a justification in the code, and arguably a settings key.
- [ ] **Restart / resume.** A crashed batch currently means re-running configurations that already have
      valid `.out` files. Skipping completed ones would save a lot of wall time.
- [ ] **Duplicate parse block** in `get_multipole_statistics_scee` — the same glob-and-read runs twice,
      doubling file reads over ~600 outputs. Harmless, easy.
- [ ] **Path asymmetry in `run_md3`** — pure uses an absolute topology path with a relative run name,
      mixture does the reverse.
- [ ] **The `pgrep` throttle** counts every `g09` and `mdrun` on the machine, including unrelated jobs.
- [ ] **Regression test on a known system.** One molecule with published numbers, run end to end, results
      committed — so a future change that breaks a dipole is visible immediately.
- [ ] **Requirements file / environment record** (GROMACS version, GromacsWrapper, RDKit, pyedr, Gaussian
      revision) so the pipeline can be rebuilt after submission.

## 4. Decisions — these are the ones for Leo

These are not bugs; they are scope questions where the answer changes what gets built and in what order.

- [ ] **What is the deliverable?** A reproducible archive that documents the thesis work, or a tool the
      group is expected to run on new molecules after submission? Almost everything in §3 is optional under
      the first reading and mandatory under the second.
- [ ] **How far does mixture support need to go before submission?** Full solute-in-solvent SCEE, or MD
      only with the QM side left as future work?
- [ ] **Should the dielectric correction move inside the iteration**, or is the post-hoc correction the
      defensible version for the thesis?
- [ ] **Is the polar-hydrogen LJ floor defensible as-is**, and where did those values originate?
- [ ] **How much of the automation belongs in the methods chapter** versus a repository README — i.e. how
      much of [[SCEE_Code_Overview]] is duplicating text he has already read.

---

## Suggested priority

If the code needs to keep producing thesis numbers while being tidied:

1. §1 items — quick, and two of them touch data you may already have generated.
2. The §4 scope decision — it determines whether §3 is worth doing at all.
3. §2 items that block remaining planned runs (dummy-hydrogen UA work in particular).
4. §3, in whatever order the answer to §4 implies.
