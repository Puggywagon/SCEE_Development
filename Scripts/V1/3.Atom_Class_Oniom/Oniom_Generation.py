#!/usr/bin/python3
import pandas as pd
import re
import Atoms

Atoms=Atoms.Atoms()

class Oniom_Generation(object):
    def __init__(self):
        pass
################################################################################
    def Calc_Heavys(self,Solute):
        txt=Solute
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge','mass']
        
        Gros=[]
        Masses=[]
        for index, row in atoms.iterrows():
            Gro=atoms.iloc[index]['at_name']
            Mass=atoms.iloc[index]['mass']
            Gros.append(Gro)
            Masses.append(Mass)
        Heavy_Atoms, Total_Atoms=Atoms.Atom_Types(Gros,Masses)
        # We need qmax both for in the oniom generation and in the dipole moment scaling after the oniom generation step in the Run_SCEE script
        qilist=[]
        for index, row in atoms.iterrows():
            qlist=atoms.iloc[index]['charge']
            qilist.append(abs(qlist))
        qmax=max(qilist)
            
        return Heavy_Atoms,Total_Atoms,qmax
################################################################################
    def QM_Inputs(self,Solute,Oniom,dummy,split,qr1,qr2,qr3):
        txt=Solute
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge','mass']
        
        Gros=[]
        Masses=[]
        for index, row in atoms.iterrows():
            Gro=atoms.iloc[index]['at_name']
            Mass=atoms.iloc[index]['mass']
            Gros.append(Gro)
            Masses.append(Mass)
            if Gro_Atom_Types == 'OW':
                spicy=row['id']
            else:
                spicy=0
        Total_Atoms, Gro_List, Gaus_List, Dummy_List=Atoms.Atom_Types(Gros,Masses)
        
        qilist=[]
        for index, row in atoms.iterrows():
            qlist=atoms.iloc[index]['charge']
            qilist.append(abs(qlist))
        qmax=max(qilist)
        
        if spicy == 0:
            Central_Atom=qilist.index(qmax)+1
        else:
            Central_Atom=spicy+1
    
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms:.0f} {Central_Atom:.0f}\n')
        
        
        # This bit I will need a little help with. It is to do with the inconsistency of what is in the [ atomtypes ], I have seen versions with the dummy column and others without it. Variations causes errors in the atom_types dataframe. We need the [ atomtypes ] section as we need to extract sigma and epsilon from this list. So we will need some way of telling users to provide us access to those numbers ether through having an atom types section in the topology or giving us the forcefield file but this needs to be consistent... We have preassigned some of these however.
        if dummy == 'yes' and split == 'yes':
            f2 = open(f'oplsaaff.itp')
            forcefield = f2.read()        
            f2.close()
            atomtype_dict=self.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
        elif dummy == 'yes' and split == 'no':
            atomtype_dict=sel        f.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
        elif dummy == 'no' and split == 'yes':
            f2 = open(f'oplsaaff.itp')
            txt = f2.read()        
            f2.close()
            atomtype_dict=self.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','mass','q','Dummy_Bead','Sigma','Epsilon']
        else:
            atomtype_dict=self.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','mass','q','Dummy_Bead','Sigma','Epsilon']

        Sigma_List=[]
        Epsilon_List=[]
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            if Gro_Atom_Types == 'OW':
                Sigma=matching_type.iloc[0]['Sigma']*10
                Epsilon=matching_type.iloc[0]['Epsilon']
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            elif Gro_Atom_Types == 'HO':
                Sigma=0.2673
                Epsilon=0.0418
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            elif Gro_Atom_Types == 'HN':
                Sigma=0.2673
                Epsilon=0.0418
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            elif Gro_Atom_Types == 'HW':
                Sigma=0.2673
                Epsilon=0.0418
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            else:
                Sigma=matching_type.iloc[0]['Sigma']
                Epsilon=matching_type.iloc[0]['Epsilon']
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            
            
        for G,Gro,Gaus,Dummy,qi,Sigma,Epsilon in zip(Gros,Gro_List,Gaus_List,Dummy_List,qilist,Sigma_List,Epsilon_List):
            wi=qi/qmax            
            if G == 'OW':
                q1=-qmax*qr1  
                q2=-qmax*qr2      
                q3=-qmax*qr3
            elif G == 'MW':
                q1=0      
                q2=0  
                q3=0
            elif not G == 'MW' and not G == 'OW' and wi == 1:
                q1=qmax*qr1       
                q2=qmax*qr2        
                q3=qmax*qr3
            else:
                q1=wi*qmax*qr1         
                q2=wi*qmax*qr2       
                q3=wi*qmax*qr3

            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus:2s} {Gro:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus:2s} {Gro:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus:2s} {Gro:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')   
                           
        #Note we are going to have issues for when we start adding Cl, Na etc due to the formating above...                
        with open(Oniom, 'a') as file:
            file.write("\n")
################################################################################
    def MM_Inputs(self,solvent,Oniom,split,qr1,qr2,qr3,qmax):
        txt=solvent
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge','mass'] 
        
        Gros=[]
        Masses=[]
        for index, row in atoms.iterrows():
            Gro=atoms.iloc[index]['at_name']
            Mass=atoms.iloc[index]['mass']
            Gros.append(Gro)
            Masses.append(Mass)
            if Gro_Atom_Types == 'OW':
                spicy=row['id']
            else:
                spicy=0
        Total_Atoms, Gro_List, Gaus_List=Atoms.Atom_Types(Gros,Masses)
        
        qilist=[]
        for index, row in atoms.iterrows():
            qlist=atoms.iloc[index]['charge']
            qilist.append(abs(qlist))
        qmax=max(qilist)
        
        if spicy == 0:
            Central_Atom=qilist.index(qmax)+1
        else:
            Central_Atom=spicy+1
    
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms:.0f} {Central_Atom:.0f}\n')
            
            
         # Unless we tell users to include sigma and epsilon in the [ atoms ] we still need a way of extracting these for the molecules? I dunno if that is something that works for gromacs. I also dunno how we even going about generalising this for other software.
        if split == 'yes':
            f2 = open(f'oplsaaff.itp')
            txt = f2.read()        
            f2.close()
            atomtype_dict=self.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','mass','q','ptype','Sigma','Epsilon']
        else:
            atomtype_dict=self.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','mass','q','ptype','Sigma','Epsilon']
                                  
        Sigma_List=[]
        Epsilon_List=[]
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            if Gro_Atom_Types == 'OW':
                Sigma=matching_type.iloc[0]['Sigma']*10
                Epsilon=matching_type.iloc[0]['Epsilon']
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            elif Gro_Atom_Types == 'HO':
                Sigma=0.2673
                Epsilon=0.0418
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            elif Gro_Atom_Types == 'HN':
                Sigma=0.2673
                Epsilon=0.0418
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            elif Gro_Atom_Types == 'HW':
                Sigma=0.2673
                Epsilon=0.0418
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            else:
                Sigma=matching_type.iloc[0]['Sigma']
                Epsilon=matching_type.iloc[0]['Epsilon']
                Sigma_List.append(Sigma)
                Epsilon_List.append(Epsilon)
            
            
        for G,Gro,Gaus,qi,Sigma,Epsilon in zip(Gros,Gro_List,Gaus_List,qilist,Sigma_List,Epsilon_List):
            wi=qi/qmax            
            if G == 'OW':
                q1=-qmax*qr1  
                q2=-qmax*qr2      
                q3=-qmax*qr3
            elif G == 'MW':
                q1=0      
                q2=0  
                q3=0
            elif not G == 'MW' and not G == 'OW' and wi == 1:
                q1=qmax*qr1       
                q2=qmax*qr2        
                q3=qmax*qr3
            else:
                q1=wi*qmax*qr1         
                q2=wi*qmax*qr2       
                q3=wi*qmax*qr3

            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus:2s} {Gro:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus:2s} {Gro:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus:2s} {Gro:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')   
        #Note we are going to have issues for when we start adding Cl, Na etc due to the formating above...                
        with open(Oniom, 'a') as file:
            file.write("\n")
#################################################################################            
    def get_atomtypes(self,txt):
        atomtype_dict = []                    
        tmp = re.findall(r'\[ *atomtypes *\] *\n+(.*?)^\s*$', txt, flags= re.MULTILINE | re.DOTALL)    
        if (len(tmp) > 0):
            lines = tmp[0].split('\n')
            print(lines[0])
            for line in tmp[0].split('\n'):
                if (len(line)>0 and line[0] != ';'):
                    data = line.split()
                    atomtype_dicts = {'name': data[0],
                                 'bond_type': data[1],
                                 'atomic_number': int(data[2]), #Dunno why our script doesn't pick this up or flag it as missing?
                                 'mass': float(data[3]),
                                 'charge': float(data[4]),
                                 'ptype': data[5],
                                 'sigma': float(data[6]),
                                 'epsilon': float(data[7])}
                    atomtype_dict.append(atomtype_dicts)
        return atomtype_dict
          
#################################################################################
    def get_atoms(self,txt):
        tmp = re.findall(r'\[ *atoms *\] *\n+(.*?)^\s*$', txt, flags= re.MULTILINE | re.DOTALL)
        atom_dict = []
        if (len(tmp) > 0):
            lines = tmp[0].split('\n')
            print(lines[0])
            for line in tmp[0].split('\n'):
                if (len(line)>0 and line[0] != ';'):
                    data = line.split()
                    atom_dicts = {'id': float(data[0]),
                             'at_type': data[1],
                             'res nr': float(data[2]),
                             'residu name': data[3],
                             'at name': data[4],
                             'cg nr': data[5],
                             'mass': float(data[6]),
                             'charge': float(data[7])}
                    atom_dict.append(atom_dicts)
        return atom_dict 
#################################################################################
    def get_included_files(self,txt):
        itp_list = re.findall(r'^#include\s+"([./\w]+)" *', txt, re.MULTILINE)
        return itp_list        
#################################################################################        
    def Counting_Molecules(self,Topology,Oniom,Solvent_Molecules):
        with open(Oniom, 'a') as file:
                file.write(f'{Solvent_Molecules:.0f}\n')
################################################################################
    def Gen_File(self,Oniom,Configurations, system_title, Cut_Off_Radius):
        with open(Oniom, 'w') as file:
            file.write(f'{Configurations:.0f}\n')
            file.write(f'{system_title}\n')
            file.write(f'{Cut_Off_Radius:.1f}\n\n')
################################################################################
