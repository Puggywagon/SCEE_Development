#!/usr/bin/python3
import math
import os
from dataclasses import dataclass, field
from typing import List, Tuple
 
 
# Constants
SIGMA_TO_AMBER = 0.561231   # σ_nm * 10 (Å) → R_min/2 in Å,
                            # i.e. R_min/2 = σ * 2^(1/6) / 2
KJ_TO_KCAL = 4.184          # 1 kcal = 4.184 kJ; ε is divided by this
NM_TO_ANG = 10.0            # nm → Å for coordinate output
 
 
# ---------------------------------------------------------------------------
# Data classes — represent parsed oniom.inp content
# ---------------------------------------------------------------------------
 
@dataclass
class AtomSpec:
    """A single atom's spec line from oniom.inp."""
    elem: str            # 1-2 char element symbol (e.g. "C", "O", "Bq")
    type_: str           # atom type tag (e.g. "OH", "CA", "MW")
    sigma: float         # converted to Amber units (R_min/2 in Å)
    epsilon: float       # converted to kcal/mol
    q1: str              # charge as string (preserves exact formatting)
    q2: str
    q3: str
    is_dummy: bool = False  # only meaningful for QM atoms
 
    def head(self, charge_str: str, layer_flag: str) -> str:
        """Build the 'Elem-Type-q-flag' header used in coordinate lines.
        layer_flag is "  0" for QM atoms and " -1" for MM atoms."""
        return (f"{self.elem.strip()}-"
                f"{self.type_.strip()}-"
                f"{charge_str.strip()}"
                f"{layer_flag}")
 
 
@dataclass
class OniomSetup:
    """Parsed contents of oniom.inp (everything except per-configuration data)."""
    n_conf: int
    outfile_prefix: str
    rcut: float                      # cutoff radius in nm
    qm_atoms: List[AtomSpec]         # QM (central) molecule atoms
    qm_central: int                  # 1-based index of central QM atom
    mm_atoms: List[AtomSpec]         # MM (bath) molecule atom specs (per molecule)
    mm_central: int                  # 1-based index of central MM atom
    n_mm_mols: int                   # number of MM molecules in the box
 
    @property
    def n_dummies(self) -> int:
        """Count of dummy QM atoms (excluded from coordinate/connectivity output)."""
        return sum(1 for a in self.qm_atoms if a.is_dummy)
 
 
# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
 
def _parse_atom_line(line: str, has_dummy_flag: bool) -> AtomSpec:
    """Parse one atom spec line from oniom.inp.
 
    Fortran format strings:
        QM: (a2,x,a2,2f7.4,3a9,x,i1)
        MM: (a2,x,a2,2f7.4,3a9)
 
    Python is more forgiving: we parse by whitespace splitting. This means
    we don't enforce the Fortran's column-width discipline — values just
    need to be space-separated.
    """
    tokens = line.split()
    if has_dummy_flag:
        if len(tokens) != 8:
            raise ValueError(
                f"Expected 8 tokens in QM atom line, got {len(tokens)}: {line!r}"
            )
        elem, type_, sigma_s, eps_s, q1, q2, q3, dum = tokens
        is_dummy = (int(dum) == 1)
    else:
        if len(tokens) != 7:
            raise ValueError(
                f"Expected 7 tokens in MM atom line, got {len(tokens)}: {line!r}"
            )
        elem, type_, sigma_s, eps_s, q1, q2, q3 = tokens
        is_dummy = False
 
    sigma = float(sigma_s) * SIGMA_TO_AMBER
    epsilon = float(eps_s) / KJ_TO_KCAL
 
    return AtomSpec(
        elem=elem,
        type_=type_,
        sigma=sigma,
        epsilon=epsilon,
        q1=q1,
        q2=q2,
        q3=q3,
        is_dummy=is_dummy,
    )
 
 
def parse_oniom_inp(path: str = "oniom.inp") -> OniomSetup:
    """Parse oniom.inp into an OniomSetup object."""
    with open(path) as f:
        # Read into a list, stripping right-hand whitespace but preserving
        # original lines so blank-line sentinels are visible.
        lines = [line.rstrip("\n") for line in f]
 
    # Drop trailing empty lines (file may end with a newline)
    while lines and lines[-1].strip() == "":
        lines.pop()
 
    # Use a generator-like cursor so we can advance through the file.
    idx = 0
 
    def next_nonblank() -> str:
        nonlocal idx
        while idx < len(lines) and lines[idx].strip() == "":
            idx += 1
        if idx >= len(lines):
            raise ValueError("Unexpected end of oniom.inp")
        line = lines[idx]
        idx += 1
        return line
 
    # Header
    n_conf = int(next_nonblank().split()[0])
    outfile_prefix = next_nonblank().split()[0]
    rcut = float(next_nonblank().split()[0])
 
    # QM block: header line "N_atoms central_atom", then N_atoms lines
    qm_header = next_nonblank().split()
    so_atoms = int(qm_header[0])
    qm_central = int(qm_header[1])
 
    qm_atom_specs: List[AtomSpec] = []
    for _ in range(so_atoms):
        line = next_nonblank()
        qm_atom_specs.append(_parse_atom_line(line, has_dummy_flag=True))
 
    # MM block: header line "N_atoms central_atom", then N_atoms lines
    mm_header = next_nonblank().split()
    sv_atoms = int(mm_header[0])
    mm_central = int(mm_header[1])
 
    mm_atom_specs: List[AtomSpec] = []
    for _ in range(sv_atoms):
        line = next_nonblank()
        mm_atom_specs.append(_parse_atom_line(line, has_dummy_flag=False))
 
    # Trailing: number of MM molecules
    n_mm_mols = int(next_nonblank().split()[0])
 
    return OniomSetup(
        n_conf=n_conf,
        outfile_prefix=outfile_prefix,
        rcut=rcut,
        qm_atoms=qm_atom_specs,
        qm_central=qm_central,
        mm_atoms=mm_atom_specs,
        mm_central=mm_central,
        n_mm_mols=n_mm_mols,
    )
 
 
# ---------------------------------------------------------------------------
# .gro reading
# ---------------------------------------------------------------------------
 
def parse_gro(path: str) -> Tuple[List[Tuple[float, float, float]], Tuple[float, float, float]]:
    """Read a .gro file and return (coords, box).
 
    coords is a list of (x, y, z) tuples for every atom, in file order.
    box is (Lx, Ly, Lz).
 
    Fortran read format for atom lines: (i5,2a5,i5,3f8.3,a)
    We use the same column-based parsing for the coordinate fields, since
    that's the only reliable way to read .gro files (whitespace splitting
    can fail when atoms have multi-digit indices that bleed into adjacent
    fields).
    """
    with open(path) as f:
        lines = f.readlines()
 
    # Line 0: title (ignored)
    # Line 1: atom count
    n_atoms = int(lines[1].strip())
 
    coords: List[Tuple[float, float, float]] = []
    for i in range(2, 2 + n_atoms):
        line = lines[i]
        # .gro coordinate format: residue_nr (5), residue_name (5),
        # atom_name (5), atom_nr (5), x (8.3f), y (8.3f), z (8.3f).
        # Columns are 0-indexed: x starts at 20, y at 28, z at 36.
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        coords.append((x, y, z))
 
    # Box line (after all atoms)
    box_line = lines[2 + n_atoms]
    parts = box_line.split()
    Lx = float(parts[0])
    Ly = float(parts[1])
    Lz = float(parts[2])
 
    return coords, (Lx, Ly, Lz)
 
 
# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
 
def _minimum_image(d: float, L: float) -> float:
    """Apply minimum-image PBC to a single displacement.
    Matches the Fortran's if/elif/else structure exactly."""
    if d > 0.5 * L:
        return d - L
    elif -d > 0.5 * L:
        return d + L
    else:
        return d
 
 
def _centred_position(atom_pos: Tuple[float, float, float],
                      centre_pos: Tuple[float, float, float],
                      box: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """Return atom_pos relative to centre_pos with minimum-image PBC,
    converted from nm to Angstrom."""
    Lx, Ly, Lz = box
    dx = _minimum_image(atom_pos[0] - centre_pos[0], Lx) * NM_TO_ANG
    dy = _minimum_image(atom_pos[1] - centre_pos[1], Ly) * NM_TO_ANG
    dz = _minimum_image(atom_pos[2] - centre_pos[2], Lz) * NM_TO_ANG
    return (dx, dy, dz)
 
 
def _displacement_nm(atom_pos: Tuple[float, float, float],
                     centre_pos: Tuple[float, float, float],
                     box: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """Return PBC-corrected displacement in nm (without nm→Å conversion).
    Used for distance calculation."""
    Lx, Ly, Lz = box
    return (
        _minimum_image(atom_pos[0] - centre_pos[0], Lx),
        _minimum_image(atom_pos[1] - centre_pos[1], Ly),
        _minimum_image(atom_pos[2] - centre_pos[2], Lz),
    )
 
 
# ---------------------------------------------------------------------------
# Per-configuration processing
# ---------------------------------------------------------------------------
 
def _format_coord_line(head: str, x: float, y: float, z: float, layer_tag: str) -> str:
    """Format one ONIOM coordinate line.
    Matches Fortran: (a17,x,3f8.3,a2)
    e.g.  'O -OH-  -0.64600  0   0.000   0.000   0.000 H'
    """
    return f"{head:17s} {x:8.3f}{y:8.3f}{z:8.3f}{layer_tag:2s}\n"
 
 
def _format_chg_line(x: float, y: float, z: float, q: str) -> str:
    """Format one charge-only line.
    Matches Fortran: (3f8.3,x,a10) — charge is right-justified in a 10-char field.
    """
    return f"{x:8.3f}{y:8.3f}{z:8.3f} {q:>10s}\n"
 
 
def process_configuration(setup: OniomSetup, conf_num: int) -> None:
    """Process one configuration: read conf_{N}.gro, write the 6 output files."""
    datafile = f"conf_{conf_num}.gro"
    coords, box = parse_gro(datafile)
 
    n_qm = len(setup.qm_atoms)
    n_mm_per_mol = len(setup.mm_atoms)
    expected_n_atoms = n_qm + n_mm_per_mol * setup.n_mm_mols
 
    if len(coords) != expected_n_atoms:
        raise ValueError(
            f"Atom count mismatch in {datafile}: "
            f"expected {expected_n_atoms} (={n_qm} QM + {n_mm_per_mol}*{setup.n_mm_mols} MM), "
            f"got {len(coords)}"
        )
 
    # Split into QM and MM coordinates
    qm_coords = coords[:n_qm]
    mm_coords = coords[n_qm:]  # length = n_mm_per_mol * n_mm_mols
 
    # Centre: 1-based to 0-based for indexing
    qm_centre = qm_coords[setup.qm_central - 1]
 
    # Output filenames
    prefix = setup.outfile_prefix
    out_main = [
        f"{prefix}_c{conf_num}_q1.inp",
        f"{prefix}_c{conf_num}_q2.inp",
        f"{prefix}_c{conf_num}_q3.inp",
    ]
    out_chg = [
        f"{prefix}_c{conf_num}_q1_chg.inp",
        f"{prefix}_c{conf_num}_q2_chg.inp",
        f"{prefix}_c{conf_num}_q3_chg.inp",
    ]
 
    # Open all six output files
    f_main = [open(p, "w") for p in out_main]
    f_chg = [open(p, "w") for p in out_chg]
 
    try:
        # --- Write QM atom coordinates (excluding dummies) ---
        for i, atom in enumerate(setup.qm_atoms):
            new_pos = _centred_position(qm_coords[i], qm_centre, box)
            if atom.is_dummy:
                continue
            qs = [atom.q1, atom.q2, atom.q3]
            for k in range(3):
                head = atom.head(qs[k], layer_flag="  0")
                f_main[k].write(_format_coord_line(head, *new_pos, layer_tag=" H"))
 
        # --- Find MM molecules within cutoff ---
        flag = [False] * setup.n_mm_mols
        cnt = 1  # count includes the central QM molecule
        largest = 0.0
        furthest_idx = -1
 
        for mol_idx in range(setup.n_mm_mols):
            # Index of this molecule's central atom in mm_coords
            j = mol_idx * n_mm_per_mol + (setup.mm_central - 1)
            mm_centre = mm_coords[j]
            dx, dy, dz = _displacement_nm(mm_centre, qm_centre, box)
            rrr = math.sqrt(dx * dx + dy * dy + dz * dz)
 
            if rrr <= setup.rcut:
                cnt += 1
                flag[mol_idx] = True
                if rrr > largest:
                    largest = rrr
                    furthest_idx = mol_idx
 
        # --- Even-count workaround (Gaussian bug) ---
        # Print the largest distance and the (1-based) index of the furthest
        # molecule before any removal — matches Fortran ordering.
        print(f"{largest:.6f}   {furthest_idx + 1}")
        if cnt % 2 == 0 and furthest_idx >= 0:
            flag[furthest_idx] = False
            cnt -= 1
            print("REMOVED")
 
        # --- Write MM atom coordinates for molecules within cutoff ---
        for mol_idx in range(setup.n_mm_mols):
            if not flag[mol_idx]:
                continue
            for ii, atom in enumerate(setup.mm_atoms):
                jj = mol_idx * n_mm_per_mol + ii
                new_pos = _centred_position(mm_coords[jj], qm_centre, box)
                qs = [atom.q1, atom.q2, atom.q3]
                for k in range(3):
                    head = atom.head(qs[k], layer_flag=" -1")
                    f_main[k].write(_format_coord_line(head, *new_pos, layer_tag=" L"))
                    f_chg[k].write(_format_chg_line(*new_pos, qs[k]))
 
        # Blank line separator after coordinates
        for f in f_main + f_chg:
            f.write("\n")
 
        # --- Connectivity: QM atoms (non-dummies) ---
        for i, atom in enumerate(setup.qm_atoms):
            if atom.is_dummy:
                continue
            for f in f_main:
                f.write(f"{i + 1:5d}\n")
 
        # --- Connectivity: MM atoms (offset by dummy count) ---
        # Fortran: for i in 2..cnt: write (i-2)*Sv_Atoms + So_Atoms + ii - Ndum
        # i in 2..cnt means cnt-1 included molecules, indexed by (i-2) = 0..cnt-2
        ndum = setup.n_dummies
        included_count = cnt - 1  # number of included MM molecules
        for m in range(included_count):  # m = 0, 1, ..., included_count - 1
            for ii in range(1, n_mm_per_mol + 1):  # 1-based atom index within molecule
                atom_index = m * n_mm_per_mol + n_qm + ii - ndum
                for f in f_main:
                    f.write(f"{atom_index:5d}\n")
 
        # Blank line separator after connectivity
        for f in f_main:
            f.write("\n")
 
        # --- VDW parameters section (only in main files) ---
        # The Fortran uses list-directed output which produces idiosyncratic
        # mixed decimal/scientific formatting. We use a cleaner fixed-decimal
        # format here — Gaussian parses both equivalently.
        for atom in setup.qm_atoms:
            for f in f_main:
                f.write(f" VDW {atom.type_:2s}  {atom.sigma:13.9f}  {atom.epsilon:13.9f}\n")
        for atom in setup.mm_atoms:
            for f in f_main:
                f.write(f" VDW {atom.type_:2s}  {atom.sigma:13.9f}  {atom.epsilon:13.9f}\n")
 
        # Trailing blank line in main files
        for f in f_main:
            f.write("\n")
 
        print(f"Total solvent molecules =        {cnt - 1}")
        print(f"Total atoms =        {(cnt - 1) * n_mm_per_mol + n_qm}")
 
    finally:
        for f in f_main + f_chg:
            f.close()
 
 
# ---------------------------------------------------------------------------
# Top-level
# ---------------------------------------------------------------------------
 
class ShellOniom:
    """Python port of Shell_Oniom.f90.
 
    Usage:
        shell = ShellOniom()
        shell.process('oniom.inp')
 
    Or equivalently from the command line:
        python Shell_Oniom.py
    """
 
    def process(self, oniom_inp_path: str = "oniom.inp") -> None:
        """Process all configurations specified in oniom.inp."""
        setup = parse_oniom_inp(oniom_inp_path)
        print(f"Processing {setup.n_conf} configurations")
        print(f"  prefix: {setup.outfile_prefix}")
        print(f"  cutoff: {setup.rcut} nm")
        print(f"  QM atoms: {len(setup.qm_atoms)} (central index: {setup.qm_central}, dummies: {setup.n_dummies})")
        print(f"  MM atoms per molecule: {len(setup.mm_atoms)} (central index: {setup.mm_central})")
        print(f"  MM molecules in box: {setup.n_mm_mols}")
 
        for cf in range(1, setup.n_conf + 1):
            datafile = f"conf_{cf}.gro"
            print(datafile)
            process_configuration(setup, cf)
 

