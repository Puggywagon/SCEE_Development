#!/usr/bin/python3

import gromacs
import socket
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem



################################################################################
################################################################################
################################################################################
class SCEE_GROMACS(object):

    def __init__(self):
        hostname = socket.gethostname()
        self.param = {}

        print(hostname)


################################################################################
################################################################################
    def __setitem__(self, key, value):
        self.param[key] = value

        
    def __getitem__(self, key):
        return self.param[key]

    
################################################################################
################################################################################    
    def AA_Structure(self,mol,molecule,resname): # Will need different structure for if we include the option to gen UA solvents but I dunno how we do that for mixtures
        mol = Chem.MolFromSmiles(mol)
        # Generate 3D coordinates
        mol = Chem.AddHs(mol)
        Chem.AllChem.EmbedMolecule(mol)
        Chem.AllChem.UFFOptimizeMolecule(mol)
        print('here')
        print(f'mol={mol}; molecule = {molecule}')
        
        img = Chem.Draw.MolToImage(mol)
        print(img)
# LL: plotting does not work for me        
#        plt.imshow(img)
#        plt.axis("off")
#        plt.savefig(f'{molecule}.png', dpi=300)  
#        plt.clf()
        print('here 2')
    
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
        Gro_File=f"{molecule}.gro"

        print('Gro_Builder.AA_Structure() complete')
        
        return Gro_File


################################################################################
################################################################################
################################################################################
