#!/usr/bin/python3
import yaml
import math
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
################################################################################
# Avogadro's constant (SI exact value, post-2019 redefinition), mol^-1
AVOGADRO = 6.02214076e23
 
 
@dataclass
class StateConditions:
    temperatures: list      # K, e.g. [298] or [298, 323, 348]
    pressures: list         # bar, e.g. [1] or [1, 10, 100]
    replicas: int           # number of replicas; runs as replica_1 .. replica_N
 
 
@dataclass
class Molecule:
    system_title: str
    gas_dipole: float                       # D (Debye)
    dielectric_constant: float              # experimental, dimensionless
    refractive_index: float                 # dimensionless
    sol_keyword: str                        # Gaussian SCRF keyword
    density: float                          # g/cm^3 (experimental)
    molar_mass: float                       # g/mol
    aa_surround: bool                       # True for AA bath, False for UA
    solvent_smiles: Optional[str]           # None if user provides solvent .gro
    mixture: bool                           # True if a solute is also present
    solute_smiles: Optional[str]            # None if pure OR if .gro provided
 
    # Derived at load time
    calculated_dielectric: float = field(init=False)
 
    def __post_init__(self):
        # Optical (high-frequency) dielectric correction:
        #     eps_calc = eps_exp - n^2 + 1
        # Used in PCM2 calculations.
        self.calculated_dielectric = round(
            self.dielectric_constant - self.refractive_index ** 2 + 1, 3
        )
 
 
@dataclass
class GromacsSoftware:
    ntmpi: int          # MPI threads
    ntomp: int          # OpenMP threads
    pin: str            # 'on' or 'off' — passed to gmx mdrun -pin
    pinoffset: int      # CPU offset for pinning
 
 
@dataclass
class GaussianSoftware:
    g09root_override: Optional[str]    # None = auto-detect
    scratch_dir: str                   # GAUSS_SCRDIR
    max_jobs: int                      # parallel job limit on this machine
    nproc: int                         # cores per Gaussian job (%nprocshared)
    mem: str                           # memory per job, e.g. "5GB"
 
 
@dataclass
class Software:
    gromacs: GromacsSoftware
    gaussian: GaussianSoftware
 
 
@dataclass
class Configuration:
    box_length_nm: float        # cubic box edge length, nm
    cluster_radius_nm: float    # cluster radius for dipole calculations, nm
 
 
@dataclass
class Sampling:
    n_configurations: int       # configurations extracted from MD trajectory
 
 
@dataclass
class ChargeScaling:
    qr1: float
    qr2: float
    qr3: float
 
 
@dataclass
class ElectronicStructureLevel:
    method: str
    basis: str
 
 
@dataclass
class ElectronicStructure:
    v0: ElectronicStructureLevel   # low-level optimisation
    v1: ElectronicStructureLevel   # high-level single point
 
 
@dataclass
class AdvancedSettings:
    configuration: Configuration
    sampling: Sampling
    charge_scaling: ChargeScaling
    electronic_structure: ElectronicStructure
 
 
@dataclass
class Settings:
    state: StateConditions
    molecule: Molecule
    software: Software
    advanced: AdvancedSettings
 
    # Derived flags
    is_mixture: bool = field(init=False)
    build_solvent_gro: bool = field(init=False)
    build_solute_gro: bool = field(init=False)
 
    # Derived molecule counts (cubic box, scientific units handled inside)
    N: float = field(init=False)
    initial_molecules: int = field(init=False)
 
    def __post_init__(self):
        self.is_mixture = self.molecule.mixture
 
        self.build_solvent_gro = self.molecule.solvent_smiles is not None
        self.build_solute_gro = (
            self.is_mixture and self.molecule.solute_smiles is not None
        )
 
        # Target molecule count for the cubic box.
        # density [g/cm^3] * volume [cm^3] / molar_mass [g/mol] * N_A [/mol]
        # L is in nm; 1 nm = 1e-7 cm, so volume = (L * 1e-7)^3 cm^3.
        L_nm = self.advanced.configuration.box_length_nm
        rho = self.molecule.density
        M = self.molecule.molar_mass
        volume_cm3 = (L_nm * 1e-7) ** 3
        self.N = (rho * volume_cm3 / M) * AVOGADRO
 
        # How many to insert: target + 5% buffer, minus the one already in
        # the box from the vacuum step.
        self.initial_molecules = int(math.ceil(self.N + 0.05 * self.N) - 1)
 
 
def load_settings(path: str = "Settings.yml") -> Settings:
    yaml_path = Path(path)
    with yaml_path.open("r") as f:
        data = yaml.safe_load(f)
 
    if not isinstance(data, dict):
        raise ValueError(
            f"{path} did not parse to a YAML mapping (dict). "
            f"Got {type(data).__name__}."
        )
 
    # --- State Conditions ---
    sc = data["State_Conditions"]
    state = StateConditions(
        temperatures=sc["Temperature"],
        pressures=sc["Pressure"],
        replicas=sc["Replicas"],
    )
 
    # --- Molecule Information ---
    mi = data["Molecule_Information"]
    molecule = Molecule(
        system_title=mi["System_Title"],
        gas_dipole=mi["Gas_Dip"],
        dielectric_constant=mi["Di_Const"],
        refractive_index=mi["Ref_Ind"],
        sol_keyword=mi["sol_keyword"],
        density=mi["density"],
        molar_mass=mi["mol_mass"],
        aa_surround=bool(mi["AA_Surround"]),
        solvent_smiles=mi["Solvent_Smiles"],
        mixture=bool(mi["Mixture"]),
        solute_smiles=mi["Solute_Smiles"],
    )
 
    # --- Software ---
    sw = data["Software"]
    gromacs_sw = GromacsSoftware(
        ntmpi=sw["GROMACS"]["ntmpi"],
        ntomp=sw["GROMACS"]["ntomp"],
        pin=sw["GROMACS"]["pin"],
        pinoffset=sw["GROMACS"]["pinoffset"],
    )
    gaussian_sw = GaussianSoftware(
        g09root_override=sw["Gaussian"]["g09root_override"],
        scratch_dir=sw["Gaussian"]["scratch_dir"],
        max_jobs=sw["Gaussian"]["max_jobs"],
        nproc=sw["Gaussian"]["nproc"],
        mem=sw["Gaussian"]["mem"],
    )
    software = Software(gromacs=gromacs_sw, gaussian=gaussian_sw)
 
    # --- Advanced Settings ---
    adv = data["Advanced_Settings"]
    advanced = AdvancedSettings(
        configuration=Configuration(
            box_length_nm=adv["configuration"]["box_length_nm"],
            cluster_radius_nm=adv["configuration"]["cluster_radius_nm"],
        ),
        sampling=Sampling(
            n_configurations=adv["sampling"]["n_configurations"],
        ),
        charge_scaling=ChargeScaling(
            qr1=adv["electrostatics"]["charge_scaling"]["qr1"],
            qr2=adv["electrostatics"]["charge_scaling"]["qr2"],
            qr3=adv["electrostatics"]["charge_scaling"]["qr3"],
        ),
        electronic_structure=ElectronicStructure(
            v0=ElectronicStructureLevel(
                method=adv["electronic_structure"]["v0"]["method"],
                basis=adv["electronic_structure"]["v0"]["basis"],
            ),
            v1=ElectronicStructureLevel(
                method=adv["electronic_structure"]["v1"]["method"],
                basis=adv["electronic_structure"]["v1"]["basis"],
            ),
        ),
    )
 
    return Settings(
        state=state,
        molecule=molecule,
        software=software,
        advanced=advanced,
    )
 

