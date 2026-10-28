#!/usr/bin/python3
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt


################################################################################
################################################################################    
################################################################################    
class Gro_Builder(object):
    
    def __init__(self):
        pass
################################################################################  
################################################################################
    def AA_Solvent(self, smiles, gro_name='Solvent.gro'):
        self._build_aa_gro(smiles, gro_name[:-4] if gro_name.endswith('.gro') else gro_name, 'Solv')
################################################################################
    def AA_Solute(self, smiles, gro_name='Solute.gro'):
        self._build_aa_gro(smiles, gro_name[:-4] if gro_name.endswith('.gro') else gro_name, 'Solu')
################################################################################      
    def UA_Solvent(self, smiles, gro_name='Solvent_UA.gro'):
        self._build_ua_gro(smiles, gro_name[:-4] if gro_name.endswith('.gro') else gro_name, 'Solv')
################################################################################
    @staticmethod
    def _write_xyz(path, symbols, coords_ang, comment):
        """Write a plain .xyz file (Angstrom) in the same atom order as the
        .gro, so the numbering can be checked in any molecular viewer."""
        with open(path, 'w') as f:
            f.write(f"{len(symbols)}\n{comment}\n")
            for sym, (x, y, z) in zip(symbols, coords_ang):
                f.write(f"{sym:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
################################################################################    
    def _build_aa_gro(self, smiles, molecule, resname):
        mol = Chem.MolFromSmiles(smiles)
        # Generate 3D coordinates
        mol = Chem.AddHs(mol)
        Chem.AllChem.EmbedMolecule(mol)
        Chem.AllChem.UFFOptimizeMolecule(mol)
        print(f'mol={mol}; molecule = {molecule}')
        
        img = Chem.Draw.MolToImage(mol)
        
        # LL: plotting does not work for me        
        fig, ax = plt.subplots()
        ax.imshow(img)
        ax.axis("off")
        fig.savefig(f'{molecule}.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        with open(f"{molecule}.gro", 'w') as f:
            # Write header
            f.write(f"{molecule}\n")
            f.write(f"{mol.GetNumAtoms()}\n")
    
            # Write atom section
            for i, atom in enumerate(mol.GetAtoms()):
                pos = (mol.GetConformer().GetAtomPosition(i))/10
                symbol=atom.GetSymbol()
                f.write(f"{1:5d}{resname:5s}{symbol:5s}{i+1:5d}{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}\n")
            box_size = 0.0
            f.write(f"{box_size:10.5f}{box_size:10.5f}{box_size:10.5f}")

        conf = mol.GetConformer()
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = [tuple(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        self._write_xyz(f"{molecule}.xyz", symbols, coords,
                        f"{molecule} from SMILES {smiles}; same atom order as {molecule}.gro")
        print(f"Gro_Builder built {molecule}.gro and {molecule}.xyz (AA)")


################################################################################  
    def _build_ua_gro(self, smiles, molecule, resname):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        if Chem.AllChem.EmbedMolecule(mol) == -1:
            raise RuntimeError(f"RDKit failed to embed 3D coordinates for {smiles}")
        Chem.AllChem.UFFOptimizeMolecule(mol)
    
        # Classify each atom
        heteroatoms = []      # (idx, element_symbol) — appear first
        carbons = []          # (idx, element_symbol) — appear after heteroatoms
        polar_hs = []         # (idx, 'H') — appear after all heavy atoms
        explicit_hs = []      # (idx, 'H') — aromatic/alkene Hs, appear last
    
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            sym = atom.GetSymbol()
            if sym == 'H':
                parent = atom.GetNeighbors()[0]
                if parent.GetSymbol() == 'C':
                    if parent.GetHybridization() == Chem.HybridizationType.SP3:
                        # Aliphatic H — drop entirely (absorbed into C's mass)
                        continue
                    else:
                        # Aromatic/alkene CH — keep explicit
                        explicit_hs.append((idx, 'H'))
                else:
                    # Polar H (bonded to O, N, etc.)
                    polar_hs.append((idx, 'H'))
            elif sym == 'C':
                carbons.append((idx, 'C'))
            else:
                heteroatoms.append((idx, sym))
        
        # Assemble in UA convention: heteroatoms, carbons, polar Hs, explicit Hs
        ordered = heteroatoms + carbons + polar_hs + explicit_hs
        
        # Write the .gro
        conformer = mol.GetConformer()
        with open(f"{molecule}.gro", 'w') as f:
            f.write(f"{molecule}\n")
            f.write(f"{len(ordered)}\n")
            for i, (atom_idx, name) in enumerate(ordered):
                pos = conformer.GetAtomPosition(atom_idx) / 10  # Å → nm
                f.write(
                    f"{1:5d}{resname:5s}{name:5s}{i+1:5d}"
                    f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}\n"
                )
            f.write(f"{0.0:10.5f}{0.0:10.5f}{0.0:10.5f}")
        symbols = [name for _, name in ordered]
        coords = [tuple(conformer.GetAtomPosition(idx)) for idx, _ in ordered]
        self._write_xyz(f"{molecule}.xyz", symbols, coords,
                        f"{molecule} (united-atom: aliphatic H merged into carbons) from SMILES {smiles}; same order as {molecule}.gro")
        print(f"Gro_Builder built {molecule}.gro and {molecule}.xyz (UA)")
################################################################################    
################################################################################    
