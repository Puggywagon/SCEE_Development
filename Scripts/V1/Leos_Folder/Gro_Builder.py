#!/usr/bin/python3
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt

##################################################################################################
class Gro_Builder(object):
    def __init__(self):
        pass
###################################################################################################
    def Functional_Group(self,solmol):
       
        ketone_smarts = '[C](=O)[C]'
        carboxylic_acid_smarts = '[C](=O)[O]'
        ester_smarts = '[C](=O)[O][C]'
        aldehyde_smarts = '[C](=O)[H]'
        secondary_amine_smarts = 'N([C])([C])'
        tertiary_amine_smarts = 'N([C])([C])([C])'
        primary_amine_smarts = '[NH2]'
        thiol_smarts = '[S][H]'
        ether_smarts = '[O]([C])([C])'
        secondary_alcohol_smarts = '[C](O)(C)'
        tertiary_alcohol_smarts = '[C](O)(C)(C)'
        primary_alcohol_smarts = '[CH1](O)'
        diol_smarts = '[OH][CH][CH][OH]'
        amino_acid_smarts = '[C](=O)[O][C](N)'
        phenol_smarts = '[OH]c1ccccc1'
        aniline_smarts = '[N]c1ccccc1'
        benzene_smarts = 'c1ccccc1'

    
        if solmol.HasSubstructMatch(Chem.MolFromSmarts(ester_smarts)):
            functional_groups="Ester"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(aldehyde_smarts)):
            functional_groups="Aldehyde"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(secondary_amine_smarts)) or \
             solmol.HasSubstructMatch(Chem.MolFromSmarts(tertiary_amine_smarts)) or \
             solmol.HasSubstructMatch(Chem.MolFromSmarts(primary_amine_smarts)):
            if solmol.HasSubstructMatch(Chem.MolFromSmarts(secondary_amine_smarts)):
                functional_groups="Secondary Amine"
            elif solmol.HasSubstructMatch(Chem.MolFromSmarts(tertiary_amine_smarts)):
                functional_groups="Tertiary Amine"
            else:
                functional_groups="Primary Amine"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(thiol_smarts)):
            functional_groups="Thiol"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(ether_smarts)):
            functional_groups="Ether"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(secondary_alcohol_smarts)) or \
             solmol.HasSubstructMatch(Chem.MolFromSmarts(tertiary_alcohol_smarts)) or \
             solmol.HasSubstructMatch(Chem.MolFromSmarts(primary_alcohol_smarts)):
            if solmol.HasSubstructMatch(Chem.MolFromSmarts(secondary_alcohol_smarts)):
                functional_groups="Secondary Alcohol"
            elif solmol.HasSubstructMatch(Chem.MolFromSmarts(tertiary_alcohol_smarts)):
                functional_groups="Tertiary Alcohol"
            else:
                functional_groups="Primary Alcohol"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(diol_smarts)):
            functional_groups="Diol"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(amino_acid_smarts)):
            functional_groups="Amino Acid"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(phenol_smarts)):
            functional_groups="Aromatic"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(aniline_smarts)):
            functional_groups="Aromatic"
        elif solmol.HasSubstructMatch(Chem.MolFromSmarts(benzene_smarts)):
            functional_groups="Aromatic"
        else:
            functional_groups="Carbons only"
    
        return functional_groups
###################################################################################################
    def AA_Structure(self,mol,molecule,resname): # Will need different structure for if we include the option to gen UA solvents but I dunno how we do that for mixtures
        mol = Chem.MolFromSmiles(mol)
        # Generate 3D coordinates
        mol = Chem.AddHs(mol)
        Chem.AllChem.EmbedMolecule(mol)
        Chem.AllChem.UFFOptimizeMolecule(mol)
        
        img = Chem.Draw.MolToImage(mol)
        plt.imshow(img)
        plt.axis("off")
        plt.savefig(f'{molecule}.png', dpi=300)  
        plt.clf()
    
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
        return Gro_File
#####################################################################################################################################
