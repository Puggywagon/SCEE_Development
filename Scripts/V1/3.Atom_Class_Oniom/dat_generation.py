#!/usr/bin/python3
import pandas as pd
import re
import numpy as np
import Atoms

Atoms=Atoms.Atoms()
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
            Gros=gro_df.iloc[index]['gro_atom']
            dat_atom=Atoms.Atom_Types(Gros,Masses=0)
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
########################################################################################################################            




