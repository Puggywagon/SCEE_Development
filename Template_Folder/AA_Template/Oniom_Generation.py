#!/usr/bin/python3
import pandas as pd
import re

from pathlib import Path
import Atoms

Atoms=Atoms.Atoms()

class Oniom_Generation(object):
    def __init__(self):
        pass
################################################################################
    def build_atomtype_table(self, topology_file):

        table = {}  # name -> record dict
        seen = set()
    
        def walk(path):
            path = Path(path).resolve()
            if path in seen or not path.exists():
                return
            seen.add(path)
            text = path.read_text()
            for record in self.get_atomtypes(text):
                table[record['name']] = record   # later overrides earlier
            for inc in self.get_included_files(text):
                walk((path.parent / inc).resolve())
        walk(Path(topology_file))
        return table
################################################################################
    def get_atoms_from_file(self, itp_path):
        with open(itp_path) as f:
            return self.get_atoms(f.read())

################################################################################
    def QM_Inputs(self,settings, role, Oniom):
        # Will need to add something that identifies what molecule is what?
        if role == 'pure':
            atoms_itp = 'Solvent.itp'
            topology = 'Pure.top'
        elif role == 'mixture':
            atoms_itp = 'Solute.itp'
            topology = 'Mixture.top'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        atom_dict = self.get_atoms_from_file(atoms_itp)
        atoms = pd.DataFrame(atom_dict)
        atoms.columns = ['id', 'at_type', 'res num', 'res_name', 'at_name', 'cg nr', 'charge', 'mass']
    
        
        Gros = []
        Masses = []
        spicy = 0
        for index, row in atoms.iterrows():
            Gro = atoms.iloc[index]['at_name']
            Mass = atoms.iloc[index]['mass']
            Gros.append(Gro)
            Masses.append(Mass)
            if Gro == 'OW':
                spicy = row['id']
        Gro_List, Gaus_List, Dummy_List, Total_Atoms = Atoms.Atom_Types(Gros, Masses)
    
        qilist = []
        for index, row in atoms.iterrows():
            qilist.append(atoms.iloc[index]['charge'])
        qmax_signed = max(qilist, key=abs)
        qmax = float(abs(qmax_signed))     
        
    
        if spicy == 0:
            Central_Atom = qilist.index(qmax_signed) + 1
        else:
            Central_Atom = spicy + 1
    
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms:.0f} {Central_Atom:.0f}\n')
        
        # Atomtype lookup from the full include tree
        atomtype_table = self.build_atomtype_table(topology)
    
        Sigma_List = []
        Epsilon_List = []
        for index, row in atoms.iterrows():
            at_type = row['at_type']
            record = atomtype_table.get(at_type)
            if record is None:
                raise ValueError(f"Atomtype {at_type!r} not found in {topology}")
            Gro_Atom_Types = record['bond_type'] if record['bond_type'] else at_type
            # Sigma/epsilon special-casing (water + polar Hs)
            if Gro_Atom_Types == 'OW':
                Sigma_List.append(record['sigma'] * 10)
                Epsilon_List.append(record['epsilon'])
            elif Gro_Atom_Types in ('HO', 'HN', 'HW'):
                Sigma_List.append(0.2673)
                Epsilon_List.append(0.0418)
            else:
                Sigma_List.append(record['sigma'] * 10)
                Epsilon_List.append(record['epsilon'])
                
        qr1 = settings.advanced.charge_scaling.qr1
        qr2 = settings.advanced.charge_scaling.qr2
        qr3 = settings.advanced.charge_scaling.qr3
    
        for G, Gro, Gaus, Dummy, qi, Sigma, Epsilon in zip(
            Gros, Gro_List, Gaus_List, Dummy_List, qilist, Sigma_List, Epsilon_List
        ):
            wi = qi / qmax
            if G == 'OW':
                q1 = -qmax * qr1
                q2 = -qmax * qr2
                q3 = -qmax * qr3
            elif G == 'MW':
                q1 = q2 = q3 = 0
            elif wi == 1:
                q1 = qmax * qr1
                q2 = qmax * qr2
                q3 = qmax * qr3
            else:
                q1 = wi * qmax * qr1
                q2 = wi * qmax * qr2
                q3 = wi * qmax * qr3
            
            with open(Oniom, 'a') as file:
                file.write(
                    f'{Gaus:2s} {Gro:0s} {Sigma:6.4f} {Epsilon:5.4f} '
                    f'{q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n'
                )
        
        # Trailing blank line separates QM block from MM block (or MM from molecule count)
        with open(Oniom, 'a') as file:
            file.write("\n")
        
        return Total_Atoms, qmax
################################################################################
    def MM_Inputs(self, settings, Oniom):
        # MM region — bath solvent. Source itp depends on aa_surround.
        if settings.molecule.aa_surround:
            atoms_itp = 'Solvent.itp'
        else:
            atoms_itp = 'Solvent_UA.itp'
    
        # Topology for atomtype lookup — same one used for the run
        if settings.is_mixture:
            topology = 'Mixture.top'
        else:
            topology = 'Pure.top'
    
        atom_dict = self.get_atoms_from_file(atoms_itp)
        atoms = pd.DataFrame(atom_dict)
        atoms.columns = ['id', 'at_type', 'res num', 'res_name', 'at_name', 'cg nr', 'charge', 'mass']
    
        Gros = []
        Masses = []
        spicy = 0
        for index, row in atoms.iterrows():
            Gro = atoms.iloc[index]['at_name']
            Mass = atoms.iloc[index]['mass']
            Gros.append(Gro)
            Masses.append(Mass)
            if Gro == 'OW':
                spicy = row['id']
        
        Gro_List, Gaus_List, Dummy_List, Total_Atoms = Atoms.Atom_Types(Gros, Masses)
        
        qilist = []
        for index, row in atoms.iterrows():
            qilist.append(atoms.iloc[index]['charge'])
        qmax_signed = max(qilist,key=abs)
        qmax = float(abs(qmax_signed))

        if spicy == 0:
            Central_Atom = qilist.index(qmax_signed) + 1
        else:
            Central_Atom = spicy + 1
    
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms:.0f} {Central_Atom:.0f}\n')
    
        # Atomtype lookup from the full include tree
        atomtype_table = self.build_atomtype_table(topology)
    
        Sigma_List = []
        Epsilon_List = []
        for index, row in atoms.iterrows():
            at_type = row['at_type']
            record = atomtype_table.get(at_type)
            if record is None:
                raise ValueError(f"Atomtype {at_type!r} not found in {topology}")
            Gro_Atom_Types = record['bond_type'] if record['bond_type'] else at_type
            # Sigma/epsilon special-casing (water + polar Hs)
            if Gro_Atom_Types == 'OW':
                Sigma_List.append(record['sigma'] * 10)
                Epsilon_List.append(record['epsilon'])
            elif Gro_Atom_Types in ('HO', 'HN', 'HW'):
                Sigma_List.append(0.2673)
                Epsilon_List.append(0.0418)
            else:
                Sigma_List.append(record['sigma'] * 10)
                Epsilon_List.append(record['epsilon'])

        qr1 = settings.advanced.charge_scaling.qr1
        qr2 = settings.advanced.charge_scaling.qr2
        qr3 = settings.advanced.charge_scaling.qr3
        
        for G, Gro, Gaus, qi, Sigma, Epsilon in zip(
            Gros, Gro_List, Gaus_List, qilist, Sigma_List, Epsilon_List
        ):
            wi = qi / qmax
            if G == 'OW':
                q1 = q2 = q3 = 0
            elif G == 'MW':
                q1 = -qmax * qr1
                q2 = -qmax * qr2
                q3 = -qmax * qr3
            elif wi == 1:
                q1 = qmax * qr1
                q2 = qmax * qr2
                q3 = qmax * qr3
            else:
                q1 = wi * qmax * qr1
                q2 = wi * qmax * qr2
                q3 = wi * qmax * qr3
        
            with open(Oniom, 'a') as file:
                file.write(
                    f'{Gaus:2s} {Gro:0s} {Sigma:6.4f} {Epsilon:5.4f} '
                    f'{q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n'
                )
    
#################################################################################            
    def get_atomtypes(self,txt):
        atomtype_dict = []                    
        tmp = re.findall(r'\[ *atomtypes *\] *\n+(.*?)^\s*$', txt, flags= re.MULTILINE | re.DOTALL)    
        if (len(tmp) > 0):
            lines = tmp[0].split('\n')
            for line in tmp[0].split('\n'):
                if (len(line)>0 and line[0] != ';'):
                    if (len(line)>0 and line[0] != '#'):
                        data = line.split()
                        if data[4] == 'A' or data[4] == 'D':
                            # 7-column form
                            if data[1].isdigit():
                                # atomic_number-form
                                bond_type = None
                                atomic_number = data[1]
                            else:
                                # bond_type-form
                                bond_type = data[1]
                                atomic_number = None
                            record = {
                                'name': data[0],
                                'bond_type': bond_type,
                                'atomic_number': atomic_number,
                                'mass': float(data[2]),
                                'charge': float(data[3]),
                                'ptype': data[4],
                                'sigma': float(data[5]),
                                'epsilon': float(data[6])
                            }
                        else:
                            # 8-column form
                            record = {
                                'name': data[0],
                                'bond_type': data[1],
                                'atomic_number': data[2],
                                'mass': float(data[3]),
                                'charge': float(data[4]),
                                'ptype': data[5],
                                'sigma': float(data[6]),
                                'epsilon': float(data[7])
                            }
                        atomtype_dict.append(record)
        return atomtype_dict
#################################################################################
    def get_atoms(self,txt):
        tmp = re.findall(r'\[ *atoms *\] *\n+(.*?)^\s*$', txt, flags= re.MULTILINE | re.DOTALL)
        atom_dict = []
        if (len(tmp) > 0):
            lines = tmp[0].split('\n')
            for line in tmp[0].split('\n'):
                if (len(line)>0 and line[0] != ';'):
                    data = line.split()
                    atom_dicts = {'id': float(data[0]),
                             'at_type': data[1],
                             'res nr': float(data[2]),
                             'residu name': data[3],
                             'at name': data[4],
                             'cg nr': data[5],
                             'charge': float(data[6]),
                             'mass': float(data[7])}
                    atom_dict.append(atom_dicts)
        return atom_dict 
#################################################################################
    def get_included_files(self,txt):
        itp_list = re.findall(r'^#include\s+"([./\w]+)" *', txt, re.MULTILINE)
        return itp_list        
#################################################################################        
    def Counting_Molecules(self,Oniom,Solvent_Molecules):
        with open(Oniom, 'a') as file:
                file.write(f'{Solvent_Molecules:.0f}\n')
################################################################################
    def Gen_File(self,Oniom,Configurations, Solvent, Cut_Off_Radius):
        with open(Oniom, 'w') as file:
            file.write(f'{Configurations:.0f}\n')
            file.write(f'{Solvent}\n')
            file.write(f'{Cut_Off_Radius:.1f}\n\n')
################################################################################
