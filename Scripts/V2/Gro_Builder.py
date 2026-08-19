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
    def AA_Solvent(self, smiles):
        self._build_aa_gro(smiles, 'Solvent', 'Solv')
################################################################################
    def AA_Solute(self, smiles):
        self._build_aa_gro(smiles, 'Solute', 'Solu')
################################################################################      
    def UA_Solvent(self, smiles):
        self._build_ua_gro(smiles, 'Solvent_UA', 'Solv')    
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
        print(f"Gro_Builder built {molecule}.gro (UA)")
################################################################################    
################################################################################    
