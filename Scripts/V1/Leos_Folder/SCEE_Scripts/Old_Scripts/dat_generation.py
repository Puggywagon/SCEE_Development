#!/usr/bin/python3
import pandas as pd
import re
import numpy as np
##############################################################################################################
class Dat_Generation(object):
    def __init__(self,filename):
        self.filename=filename
#############################################################################################################
    def gro_to_dat(self):
        f=open(self.filename)
        gro_file = f.readlines()        
        f.close()
        gro_file=gro_file[2:]
        gro_file=gro_file[:-1]
        gro_dict=self.get_gro(gro_file)  
        gro_df=pd.DataFrame(gro_dict)
        gro_df.columns=['Residue','gro_atom','atom_number','x','y','z','Vel x','Vel y','Vel z']
        dat_atoms=[]
        dat_xs=[]
        x_mean=gro_df['x'].mean()
        dat_ys=[]
        y_mean=gro_df['y'].mean()
        dat_zs=[]
        z_mean=gro_df['z'].mean()
        for index, row in gro_df.iterrows():
            gro_atoms=gro_df.iloc[index]['gro_atom']
            dat_atom=self.gro_to_dat_atoms(gro_atoms)
            dat_x=gro_df.iloc[index]['x']
            dat_x=(dat_x-x_mean)*10
            dat_y=gro_df.iloc[index]['y']
            dat_y=(dat_y-y_mean)*10
            dat_z=gro_df.iloc[index]['z']
            dat_z=(dat_z-z_mean)*10
            print(f'{dat_x},{dat_y},{dat_z}')
            dat_atoms.append(dat_atom)
            dat_xs.append(dat_x)
            dat_ys.append(dat_y)
            dat_zs.append(dat_z)
        return dat_atoms,dat_xs,dat_ys,dat_zs
#################################################################################            
    def get_gro(self,gro_file):
        gro_dict = []                
        for line in gro_file:
            data = line.split()
            gro_dicts = {'Residue': data[0],
                         'gro_atoms': data[1],
                         'atom_number': data[2],
                         'x': float(data[3]),
                         'y': float(data[4]),
                         'z': float(data[5]),
                         'Vel x': float(data[6]),
                         'Vel y': float(data[7]),
                         'Vel z': float(data[8])}
            gro_dict.append(gro_dicts)
        return gro_dict
####################################################################################################
    def gro_to_dat_atoms(self,gro_atoms):
    
        if gro_atoms == 'OW':
            dat_atom = 'O'
            return dat_atom 
        
        elif gro_atoms == 'OH':
            dat_atom = 'O'
            return dat_atom 
        
        elif gro_atoms == 'HO':
            dat_atom = 'H'
            return dat_atom 
        
        elif gro_atoms == 'HW':
            dat_atom = 'H'
            return dat_atom 
            
        elif gro_atoms == 'NH':
            dat_atom = 'N'
            return dat_atom 
            
        elif gro_atoms == 'HN':
            dat_atom = 'H'
            return dat_atom 
             
        elif gro_atoms == 'HC':
            dat_atom = 'H'
            return dat_atom 
            
        elif gro_atoms == 'CM':
            dat_atom = 'C'
            return dat_atom 
            
        elif gro_atoms == 'CA':
            dat_atom = 'C'
            return dat_atom 
            
        elif gro_atoms == 'CB':
            dat_atom = 'C'
            return dat_atom 
        
        elif gro_atoms == 'CC':
            dat_atom = 'C'
            return dat_atom 
            
        elif gro_atoms == 'CD':
            dat_atom = 'C'
            return dat_atom 
        
        elif gro_atoms == 'CE':
            dat_atom = 'C'
            return dat_atom 
        
        elif gro_atoms == 'CF':
            dat_atom = 'C'
            return dat_atom 
        
        elif gro_atoms == 'CG':
            dat_atom = 'C'
            return dat_atom 
            
        elif gro_atoms == 'MW':
            dat_atom = 'Bq'
            return dat_atom 
                
        elif gro_atoms == 'OK':
            dat_atom = 'O'
            return dat_atom 
            return dat_atom 
                
        elif gro_atoms == 'CK':
            dat_atom = 'C' 
            return dat_atom 
                
        elif gro_atoms == 'C!':
            dat_atom = 'C'
            return dat_atom 
                
        elif gro_atoms == 'C=':
            dat_atom = 'C'
            return dat_atom 
                
        elif gro_atoms == 'CT':
            dat_atom = 'C'
            return dat_atom 
            
        elif gro_atoms == 'HA':
            dat_atom = 'H'
            return dat_atom 
                
        elif gro_atoms == 'CT':
            dat_atom = 'C'
            return dat_atom 
                
        elif gro_atoms == 'CT_4':
            dat_atom = 'C'
            return dat_atom 
                
        elif gro_atoms == 'F':
            dat_atom = 'F'   
            return dat_atom  
                
        elif gro_atoms == 'Cl':
            dat_atom = 'Cl'
            return dat_atom 
            
        elif gro_atoms == 'OC':
            dat_atom = 'O'
            return dat_atom 
            
        elif gro_atoms == 'CO':
            dat_atom = 'C'
            return dat_atom 
            
        elif gro_atoms == 'C_2':
            dat_atom = 'C'
            return dat_atom 
            
        elif gro_atoms == 'O_2':
            dat_atom = 'O'
            return dat_atom 
            
        elif gro_atoms == 'NT':
            dat_atom = 'N'
            return dat_atom  
            
        elif gro_atoms == 'H':
            dat_atom = 'H'
            return dat_atom

        elif gro_atoms == 'CN':
            dat_atom = 'C'
            return dat_atom

        elif gro_atoms == 'CH':
            dat_atom = 'C'
            return dat_atom
########################################################################################################################            




