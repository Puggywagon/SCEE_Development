#!/usr/bin/python3
import pandas as pd
import re

class Oniom_Generation(object):
    def __init__(self):
        pass
################################################################################
    def QM_Inputs_ff(self,Solute,Oniom):
        txt=Solute
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge']
    
        f2 = open(f'TrappeFF_optimum.itp')
        ForceField = f2.read()        
        f2.close()
        atomtype_dict=self.get_atomtypes_FF(ForceField)
        atom_types=pd.DataFrame(atomtype_dict)
        atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
                           
        qilist=[]
        atom_counts=0
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            
            qilist.append(abs(qlist))
            atom_counts += 1
        qmax=max(qilist)
        Central_Atom=qilist.index(qmax)+1
        Total_Atoms_str=atom_counts
        
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms_str:.0f} {Central_Atom:.0f}\n')
        
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            Gro_Atom_Types=matching_type.iloc[0]['type']
            Gaus_Atom_Types=self.gro_to_gaus(Gro_Atom_Types)
            if Gro_Atom_Types == 'HO':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HN':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HW':
                Sigma=0.2673
                Epsilon=0.0418
            else:
                Sigma=matching_type.iloc[0]['Sigma']
                Epsilon=matching_type.iloc[0]['Epsilon']
            qi=matching_type.iloc[0]['q']
            wi=qi/qmax
            if qi == qmax:
                q1=qmax*1.0        
                q2=qmax*1.4    
                q3=qmax*1.8
            else:
                q1=wi*qmax*1.0        
                q2=wi*qmax*1.4       
                q3=wi*qmax*1.8
            dummy='No'
            if dummy == 'Yes':
                Dummy=1
            elif dummy == 'No':
                Dummy=0
            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f}{Epsilon:7.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                 
        with open(Oniom, 'a') as file:
            file.write("\n")

################################################################################
    def MM_Inputs_ff(self,Solvent,Oniom):
        #Read in the data
        txt=Solvent
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge']      
    
        f2 = open(f'TrappeFF_optimum.itp')
        ForceField = f2.read()    
        f2.close()
        atomtype_dict=self.get_atomtypes_FF(ForceField)
        atom_types=pd.DataFrame(atomtype_dict)
        atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
               
        qilist=[]
        atom_counts=0
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            qilist.append(abs(qlist))
            atom_counts += 1
        qmax=max(qilist)
        Central_Atom=qilist.index(qmax)+1
        Total_Atoms_str=atom_counts
        
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms_str:.0f} {Central_Atom:.0f}\n')
   
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            Gro_Atom_Types=matching_type.iloc[0]['type']
            Gaus_Atom_Types=self.gro_to_gaus(Gro_Atom_Types)
            if Gro_Atom_Types == 'HO':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HN':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HW':
                Sigma=0.2673
                Epsilon=0.0418
            else:
                Sigma=matching_type.iloc[0]['Sigma']
                Epsilon=matching_type.iloc[0]['Epsilon']
            qi=matching_type.iloc[0]['q']
            wi=qi/qmax
            if qi == qmax:
                q1=qmax*1.0     
                q2=qmax*1.4     
                q3=qmax*1.8
            else:
                q1=wi*qmax*1.0      
                q2=wi*qmax*1.4       
                q3=wi*qmax*1.8
            if matching_type.iloc[0]['Dummy_Bead'] == str('D'):
                Dummy=1
            elif matching_type.iloc[0]['Dummy_Bead'] == str('A'):
                Dummy=0
            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f}{Epsilon:7.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
            
################################################################################
    def QM_Inputs_Water(self,Solute,Oniom):
        txt=Solute
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge']
        
        atomtype_dict=self.get_atomtypes(txt)
        atom_types=pd.DataFrame(atomtype_dict)
        atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
               
        qilist=[]
        atom_counts=0
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            Gro_Atom_Types=matching_type.iloc[0]['type']
            if Gro_Atom_Types == 'OW':
                Central_Atom=row['id']
            qilist.append(abs(qlist))
            atom_counts += 1
        qmax=max(qilist)
        Total_Atoms_str=atom_counts
    
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms_str:.0f} {Central_Atom:.0f}\n')
    
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            Gro_Atom_Types=matching_type.iloc[0]['type']
            Gaus_Atom_Types=self.gro_to_gaus(Gro_Atom_Types)
            if Gro_Atom_Types == 'OW':
                Sigma=matching_type.iloc[0]['Sigma']*10
                Epsilon=matching_type.iloc[0]['Epsilon']
            elif Gro_Atom_Types == 'HO':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HN':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HW':
                Sigma=0.2673
                Epsilon=0.0418
            else:
                Sigma=matching_type.iloc[0]['Sigma']
                Epsilon=matching_type.iloc[0]['Epsilon']
            qi=matching_type.iloc[0]['q']
            wi=qi/qmax
            if Gro_Atom_Types == 'OW':
                q1=-qmax*1.0   
                q2=-qmax*1.4      
                q3=-qmax*1.8
            elif Gro_Atom_Types == 'MW':
                q1=0      
                q2=0       
                q3=0
            else:
                q1=wi*qmax*1.0      
                q2=wi*qmax*1.4        
                q3=wi*qmax*1.8
            if matching_type.iloc[0]['Dummy_Bead'] == str('D'):
                Dummy=1
            elif matching_type.iloc[0]['Dummy_Bead'] == str('A'):
                Dummy=0
            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
            
            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f}{Epsilon:7.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                 
             
        with open(Oniom, 'a') as file:
            file.write("\n")

################################################################################
    def MM_Inputs_Water(self,Solvent,Oniom):
        #Read in the data
        txt=Solvent
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge']      
    
        atomtype_dict=self.get_atomtypes(txt)
        atom_types=pd.DataFrame(atomtype_dict)
        atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
                   
        qilist=[]
        atom_counts=0
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            Gro_Atom_Types=matching_type.iloc[0]['type']
            if Gro_Atom_Types == 'OW':
                Central_Atom=row['id']
            qilist.append(abs(qlist))
            atom_counts += 1
        qmax=max(qilist)
        Total_Atoms_str=atom_counts
        
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms_str:.0f} {Central_Atom:.0f}\n')
   
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            Gro_Atom_Types=matching_type.iloc[0]['type']
            Gaus_Atom_Types=self.gro_to_gaus(Gro_Atom_Types)
            if Gro_Atom_Types == 'OW':
                Sigma=matching_type.iloc[0]['Sigma']*10
                Epsilon=matching_type.iloc[0]['Epsilon']
            elif Gro_Atom_Types == 'HO':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HN':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HW':
                Sigma=0.2673
                Epsilon=0.0418
            else:
                Sigma=matching_type.iloc[0]['Sigma']
                Epsilon=matching_type.iloc[0]['Epsilon']
            qi=matching_type.iloc[0]['q']
            wi=qi/qmax
            if qi == qmax:
                q1=qmax*1.0        
                q2=qmax*1.4        
                q3=qmax*1.8
            else:
                q1=wi*qmax*1.0      
                q2=wi*qmax*1.4            
                q3=wi*qmax*1.8
            dummy='No'
            if dummy == 'Yes':
                Dummy=1
            elif dummy == 'No':
                Dummy=0
            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f}{Epsilon:7.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
 ################################################################################
    def QM_Inputs(self,Solute,Oniom):
        txt=Solute
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge']
        
        atomtype_dict=self.get_atomtypes(txt)
        atom_types=pd.DataFrame(atomtype_dict)
        atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
        
        qilist=[]
        atom_counts=0
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            qilist.append(abs(qlist))
            atom_counts += 1
        qmax=max(qilist)
        Central_Atom=qilist.index(qmax)+1
        Total_Atoms_str=atom_counts
    
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms_str:.0f} {Central_Atom:.0f}\n')
   
   
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            Gro_Atom_Types=matching_type.iloc[0]['type']
            Gaus_Atom_Types=self.gro_to_gaus(Gro_Atom_Types)
            if Gro_Atom_Types == 'HO':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HN':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HW':
                Sigma=0.2673
                Epsilon=0.0418
            else:
                Sigma=matching_type.iloc[0]['Sigma']
                Epsilon=matching_type.iloc[0]['Epsilon']
            qi=matching_type.iloc[0]['q']
            wi=qi/qmax
            if qi == qmax:
                q1=qmax*1.0        
                q2=qmax*1.4        
                q3=qmax*1.8
            else:
                q1=wi*qmax*1.0          
                q2=wi*qmax*1.4       
                q3=wi*qmax*1.8
            dummy='No'
            if dummy == 'Yes':
                Dummy=1
            elif dummy == 'No':
                Dummy=0
            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f}{Epsilon:7.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                 
             
        with open(Oniom, 'a') as file:
            file.write("\n")

################################################################################
    def MM_Inputs(self,Solvent,Oniom):
        #Read in the data
        txt=Solvent
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge']          

        atomtype_dict=self.get_atomtypes(txt)
        atom_types=pd.DataFrame(atomtype_dict)
        atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
                           
        qilist=[]
        atom_counts=0
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            qilist.append(abs(qlist))
            atom_counts += 1
        qmax=max(qilist)
        Central_Atom=qilist.index(qmax)+1
        Total_Atoms_str=atom_counts
    
        with open(Oniom, 'a') as file:
            file.write(f'{Total_Atoms_str:.0f} {Central_Atom:.0f}\n')
   
    
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            Gro_Atom_Types=matching_type.iloc[0]['type']
            Gaus_Atom_Types=self.gro_to_gaus(Gro_Atom_Types)
            if Gro_Atom_Types == 'HO':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HN':
                Sigma=0.2673
                Epsilon=0.0418
            elif Gro_Atom_Types == 'HW':
                Sigma=0.2673
                Epsilon=0.0418
            else:
                Sigma=matching_type.iloc[0]['Sigma']
                Epsilon=matching_type.iloc[0]['Epsilon']
            qi=matching_type.iloc[0]['q']
            wi=qi/qmax
            if qi == qmax:
                q1=qmax*1.0        
                q2=qmax*1.4        
                q3=qmax*1.8
            else:
                q1=wi*qmax*1.0        
                q2=wi*qmax*1.4        
                q3=wi*qmax*1.8
            if matching_type.iloc[0]['Dummy_Bead'] == str('D'):
                Dummy=1
            elif matching_type.iloc[0]['Dummy_Bead'] == str('A'):
                Dummy=0
            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f}{Epsilon:7.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
             
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
                                 'mass': float(data[2]),
                                 'charge': float(data[3]),
                                 'ptype': data[4],
                                 'sigma': float(data[5]),
                                 'epsilon': float(data[6])}
                    atomtype_dict.append(atomtype_dicts)
        return atomtype_dict
  #################################################################################            
    def get_atomtypes_FF(self,ForceField):
        txt=ForceField
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
                                 'mass': float(data[2]),
                                 'charge': float(data[3]),
                                 'ptype': data[4],
                                 'sigma': float(data[5]),
                                 'epsilon': float(data[6])}
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
                             'cg nr': float(data[5]),
                             'charge': float(data[6])}
                    atom_dict.append(atom_dicts)
        return atom_dict 
#################################################################################
    def get_Molecules(self,txt):
        tmp = re.findall(r'\[ *molecules *\] *\n+(.*?)^\s*$', txt, flags= re.MULTILINE | re.DOTALL)
        Molecules_dict=[]    
        if (len(tmp) > 0):
            lines = tmp[0].split('\n')
            print(lines[0])
            for line in tmp[0].split('\n'):
                if (len(line)>0 and line[0] != ';'):
                    data = line.split()
                    Molecules_dicts = {'Residue': data[0],
                                          'Number Of Molecules': float(data[1])}
                    Molecules_dict.append(Molecules_dicts)
        return Molecules_dict    
#################################################################################
    def get_included_files(self,txt):
        itp_list = re.findall(r'^#include\s+"([./\w]+)" *', txt, re.MULTILINE)
        return itp_list        
#################################################################################
    def gro_to_gaus(self,Gro_Atom_Types):
        if Gro_Atom_Types == 'OW':
            Gaus_Atom_Types = 'O'
            return Gaus_Atom_Types 
        
        elif Gro_Atom_Types == 'OH':
            Gaus_Atom_Types = 'O'
            return Gaus_Atom_Types 
        
        elif Gro_Atom_Types == 'HO':
            Gaus_Atom_Types = 'H'
            return Gaus_Atom_Types 
        
        elif Gro_Atom_Types == 'HW':
            Gaus_Atom_Types = 'H'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'NH':
            Gaus_Atom_Types = 'O'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'HN':
            Gaus_Atom_Types = 'H'
            return Gaus_Atom_Types 
             
        elif Gro_Atom_Types == 'HC':
            Gaus_Atom_Types = 'H'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'CM':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'CA':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'CB':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
        
        elif Gro_Atom_Types == 'CC':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'CD':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
        
        elif Gro_Atom_Types == 'CE':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
        
        elif Gro_Atom_Types == 'CF':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
        
        elif Gro_Atom_Types == 'CG':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'MW':
            Gaus_Atom_Types = 'Bq'
            return Gaus_Atom_Types 
                
        elif Gro_Atom_Types == 'OK':
            Gaus_Atom_Types = 'O'
            return Gaus_Atom_Types 
                
        elif Gro_Atom_Types == 'CK':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 

#################################################################################        
    def Counting_Molecules(self,Topology,Oniom):
        txt=Topology
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','residue_name','at_name','cg nr','charge']  
    
        Molecules_dict=self.get_Molecules(Topology)
        Molecules=pd.DataFrame(Molecules_dict)
        Molecules.columns=['Res_name','No of Molecules']   

        MM_Mol_list=[]
        for index, row in Molecules.iterrows():
            matching_type = atoms.loc[atoms['residue_name'] == row['Res_name']]
            MM_Molecule = row['No of Molecules']
            MM_Mol_list.append(MM_Molecule)
        MM_Molecules=max(MM_Mol_list)
        with open(Oniom, 'a') as file:
                file.write(f'{MM_Molecules:.0f}\n')
#################################################################################        
    def Counting_Molecules_itp(self,Topology,Solvent,Oniom):
        txt=Solvent
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','residue_name','at_name','cg nr','charge']  
    
        Molecules_dict=self.get_Molecules(Topology)
        Molecules=pd.DataFrame(Molecules_dict)
        Molecules.columns=['Res_name','No of Molecules']  
    
        MM_Mol_list=[]
        for index, row in Molecules.iterrows():
            matching_type = atoms.loc[atoms['residue_name'] == row['Res_name']]
            MM_Molecule = row['No of Molecules']
            MM_Mol_list.append(MM_Molecule)
        MM_Molecules=max(MM_Mol_list)
        with open(Oniom, 'a') as file:
            file.write(f'{MM_Molecules:.0f}\n\n')

################################################################################
    def Gen_File(self,Oniom,Configurations, System_Title, Cut_Off_Radius):
        with open(Oniom, 'w') as file:
            file.write(f'{Configurations:.0f}\n')
            file.write(f'{System_Title}\n')
            file.write(f'{Cut_Off_Radius:.1f}\n\n')
################################################################################
