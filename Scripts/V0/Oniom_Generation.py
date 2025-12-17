#!/usr/bin/python3
import pandas as pd
import re

class Oniom_Generation(object):
    def __init__(self):
        pass
################################################################################
    def Calc_Heavys(self,Solute,dummy,split):
        txt=Solute
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge','mass']
        
        if dummy == 'yes' and split == 'yes':
            f2 = open(f'oplsaaff.itp')
            forcefield = f2.read()        
            f2.close()
            atomtype_dict=self.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
        
        elif dummy == 'yes' and split == 'no':
            atomtype_dict=self.get_atomtypes(txt)
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
        
        print(atoms,atom_types)
        
        
        heavy_atoms=0
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            print(f"Shape of matching_type: {matching_type.shape}")
            print(f"Contents of matching_type:\n{matching_type.to_string()}")
            Gro_Atom_Types=matching_type.iloc[0]['type']
            if Gro_Atom_Types == 'OW':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'OH':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'NH':
                heavy_atoms +=1            
            elif Gro_Atom_Types == 'CM':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CA':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CB':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CC':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CD':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CE':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CF':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CG':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'OK':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CK':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'C!':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'C=':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CT':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CT':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CT_4':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'F':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'Cl':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'OC':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CO':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'C_2':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'O_2':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'NT':
                heavy_atoms +=1
            elif Gro_Atom_Types == 'CN':
                heavy_atoms += 1
        return heavy_atoms
################################################################################
    def Calc_Qmax(self,Solute,Oniom,dummy,split):
        txt=Solute
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge','mass']
        
        if dummy == 'yes' and split == 'yes':
            f2 = open(f'oplsaaff.itp')
            forcefield = f2.read()        
            f2.close()
            atomtype_dict=self.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
        
        elif dummy == 'yes' and split == 'no':
            atomtype_dict=self.get_atomtypes(txt)
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
          
        qilist=[]
        atom_counts=0
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            qilist.append(abs(qlist))
            atom_counts += 1
        Total_Atoms=atom_counts
        qmax=max(qilist)
        return qmax,Total_Atoms
################################################################################
    def QM_Inputs(self,Solute,Oniom,dummy,split,qr1,qr2,qr3,qmax):
        txt=Solute
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge','mass']
        
        if dummy == 'yes' and split == 'yes':
            f2 = open(f'oplsaaff.itp')
            forcefield = f2.read()        
            f2.close()
            atomtype_dict=self.get_atomtypes(txt)
            atom_types=pd.DataFrame(atomtype_dict)
            atom_types.columns=['name','type','MW','q','Dummy_Bead','Sigma','Epsilon']
        
        elif dummy == 'yes' and split == 'no':
            atomtype_dict=self.get_atomtypes(txt)
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
          
        qilist=[]
        atom_counts=0

        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            Gro_Atom_Types=matching_type.iloc[0]['type']
            if Gro_Atom_Types == 'OW':
                spicy=row['id']
            else:
                spicy=0
            qilist.append(abs(qlist))
            atom_counts += 1
        
        if spicy == 0:
            Central_Atom=qilist.index(qmax)+1
        else:
            Central_Atom=spicy+1
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
                q1=-qmax*qr1  
                q2=-qmax*qr2      
                q3=-qmax*qr3
            elif Gro_Atom_Types == 'MW':
                q1=0      
                q2=0  
                q3=0
            elif not Gro_Atom_Types == 'MW' and not Gro_Atom_Types == 'OW' and wi == 1:
                q1=qmax*qr1       
                q2=qmax*qr2        
                q3=qmax*qr3
            else:
                q1=wi*qmax*qr1         
                q2=wi*qmax*qr2       
                q3=wi*qmax*qr3
            if dummy == 'Yes':
                if matching_type.iloc[0]['Dummy_Bead'] == str('D'):
                    Dummy=1
                elif matching_type.iloc[0]['Dummy_Bead'] == str('A'):
                    Dummy=0
            else:
                Dummy=0
            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    if Gro_Atom_Types=='C_2':
                        Gro_Atom_Types='CO'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                    elif Gro_Atom_Types=='O_2':
                        Gro_Atom_Types='OC'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                    else:
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    if Gro_Atom_Types=='C_2':
                        Gro_Atom_Types='CO'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                    elif Gro_Atom_Types=='O_2':
                        Gro_Atom_Types='OC'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                    else:
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    if Gro_Atom_Types=='C_2':
                        Gro_Atom_Types='CO'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                    elif Gro_Atom_Types=='O_2':
                        Gro_Atom_Types='OC'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')
                    else:
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f} {Dummy:.0f}\n')   
        with open(Oniom, 'a') as file:
            file.write("\n")
################################################################################
    def MM_Inputs(self,solvent,Oniom,split,qr1,qr2,qr3,qmax):
        txt=solvent
        atom_dict=self.get_atoms(txt)
        atoms=pd.DataFrame(atom_dict)
        atoms.columns=['id','at_type','res num','res_name','at_name','cg nr','charge','mass'] 
        
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
                                  
        qilist=[]
        atom_counts=0
        
        for index, row in atoms.iterrows():
            matching_type = atom_types.loc[atom_types['name'] == row['at_type']]
            qlist=matching_type.iloc[0]['q']
            Gro_Atom_Types=matching_type.iloc[0]['type']
            if Gro_Atom_Types == 'OW':
                spicy=row['id']
            else:
                spicy=0
            qilist.append(abs(qlist))
            atom_counts += 1
        
        if spicy == 0:
            Central_Atom=qilist.index(qmax)+1
        else:
            Central_Atom=spicy+1
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
                q1=qmax*qr1        
                q2=qmax*qr2        
                q3=qmax*qr3
            else:
                q1=wi*qmax*qr1       
                q2=wi*qmax*qr2        
                q3=wi*qmax*qr3

            if q1>0 and q2>0 and q3>0:
                with open(Oniom, 'a') as file:
                    if Gro_Atom_Types=='C_2':
                        Gro_Atom_Types='CO'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
                    elif Gro_Atom_Types=='O_2':
                        Gro_Atom_Types='OC'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
                    else:
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
            elif q1==0 and q2==0 and q3==0:
                with open(Oniom, 'a') as file:
                    if Gro_Atom_Types=='C_2':
                        Gro_Atom_Types='CO'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
                    elif Gro_Atom_Types=='O_2':
                        Gro_Atom_Types='OC'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
                    else:
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')

            elif q1<0 and q2<0 and q3<0:
                with open(Oniom, 'a') as file:
                    if Gro_Atom_Types=='C_2':
                        Gro_Atom_Types='CO'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
                    elif Gro_Atom_Types=='O_2':
                        Gro_Atom_Types='OC'
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')
                    else:
                        file.write(f'{Gaus_Atom_Types:2s} {Gro_Atom_Types:0s} {Sigma:6.4f} {Epsilon:5.4f} {q1:>8.5f} {q2:>8.5f} {q3:>8.5f}\n')   
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
            Gaus_Atom_Types = 'N'
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
                
        elif Gro_Atom_Types == 'C!':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
                
        elif Gro_Atom_Types == 'C=':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
                
        elif Gro_Atom_Types == 'CT':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'HA':
            Gaus_Atom_Types = 'H'
            return Gaus_Atom_Types 
                
        elif Gro_Atom_Types == 'CT':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
                
        elif Gro_Atom_Types == 'CT_4':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
                
        elif Gro_Atom_Types == 'F':
            Gaus_Atom_Types = 'F'   
            return Gaus_Atom_Types  
                
        elif Gro_Atom_Types == 'Cl':
            Gaus_Atom_Types = 'Cl'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'OC':
            Gaus_Atom_Types = 'O'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'CO':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'C_2':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'O_2':
            Gaus_Atom_Types = 'O'
            return Gaus_Atom_Types 
            
        elif Gro_Atom_Types == 'NT':
            Gaus_Atom_Types = 'N'
            return Gaus_Atom_Types  
            
        elif Gro_Atom_Types == 'H':
            Gaus_Atom_Types = 'H'
            return Gaus_Atom_Types

        elif Gro_Atom_Types == 'CN':
            Gaus_Atom_Types = 'C'
            return Gaus_Atom_Types

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
