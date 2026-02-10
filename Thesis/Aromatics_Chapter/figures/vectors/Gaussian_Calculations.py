#!/usr/bin/python3

import glob
import gromacs
import pandas as pd
import re
from subprocess import check_call
import os
import numpy as np

import Simulations_Analysis
Analysis=Simulations_Analysis.Simulations_Analysis() 
################################################################################
class Gaussian_Calculations(object):
    def __init__(self):
        self.runname = 'TIP4P2005'       

        #self.basis_v0 = 'aug-cc-pvtz'
        #self.method_v0 = 'b3lyp'
        #self.basis_v1 = 'aug-cc-pvtz'
        #self.method_v1 = 'b3lyp'
        self.basis_v0 = '6-31+G(d,p)'
        self.method_v0 = 'b3lyp'
        self.basis_v1 = 'aug-cc-pvtz'
        self.method_v1 = 'blyp'
################################################################################
    def read_multipoles(self, filename):
        f = open(filename, 'r')
        text_file = f.read()
        f.close()

        multipole = {}
            
        re_dipole = ' Dipole moment \(field-independent basis, Debye\):(?:.*\n){2}'
        raw_data = re.findall(re_dipole, text_file, re.M)
        
        if raw_data:
            data = raw_data[-1].split('\n')
            dipole_data = data[1].split()
            multipole['x'] = float(dipole_data[1])
            multipole['y'] = float(dipole_data[3])
            multipole['z'] = float(dipole_data[5])
            multipole['total dipole'] = float(dipole_data[7])        
        else:
            multipole['x'] = np.nan
            multipole['y'] = np.nan
            multipole['z'] = np.nan
            multipole['total dipole'] = np.nan
        
        re_quad = ' Quadrupole moment \(field-independent basis, Debye-Ang\):(?:.*\n){3}'
        raw_data = re.findall(re_quad, text_file, re.M)
        if raw_data:        
            data = raw_data[-1].split('\n')
            quad_data = data[1].split()
            multipole['xx'] = float(quad_data[1])
            multipole['yy'] = float(quad_data[3])
            multipole['zz'] = float(quad_data[5])
            quad_data = data[2].split()
            multipole['xy'] = float(quad_data[1])
            multipole['xz'] = float(quad_data[3])
            multipole['yz'] = float(quad_data[5])
        else:
            multipole['xx'] = np.nan
            multipole['yy'] = np.nan
            multipole['zz'] = np.nan
            multipole['xy'] = np.nan
            multipole['xz'] = np.nan
            multipole['yz'] = np.nan

        re_oct = ' Octapole moment \(field-independent basis, Debye-Ang\*\*2\):(?:.*\n){4}'
        raw_data = re.findall(re_oct, text_file, re.M)
        if raw_data:        
            data = raw_data[-1].split('\n')
            quad_data = data[1].split()
            multipole['xxx'] = float(quad_data[1])
            multipole['yyy'] = float(quad_data[3])
            multipole['zzz'] = float(quad_data[5])
            multipole['xyy'] = float(quad_data[7])
            quad_data = data[2].split()
            multipole['xxy'] = float(quad_data[1])
            multipole['xxz'] = float(quad_data[3])
            multipole['xzz'] = float(quad_data[5])
            multipole['yzz'] = float(quad_data[7])
            quad_data = data[3].split()
            multipole['yyz'] = float(quad_data[1])
            multipole['xyz'] = float(quad_data[3])
        else:
            multipole['xxx'] = np.nan
            multipole['yyy'] = np.nan
            multipole['zzz'] = np.nan
            multipole['xyy'] = np.nan
            multipole['xxy'] = np.nan
            multipole['xxz'] = np.nan
            multipole['xzz'] = np.nan
            multipole['yzz'] = np.nan
            multipole['yyz'] = np.nan
            multipole['xyz'] = np.nan
        
        return multipole
    
################################################################################
    def get_multipole_statistics(self,filelist):
        raw_dict = {}
        case_dict = {1: 'dipole_l', 2: 'dipole_m', 3: 'dipole_h'}     
        for Q in range(1,4): 
            for filename in filelist:
                multipole = self.read_multipoles(filename)
                config = filename.split('_')[-3].replace('c', '')
                if config not in raw_dict:
                    raw_dict[config] = {}
                raw_dict[config][case_dict[Q]] = multipole['total dipole']
                
        data_dict = {'config': [], 'dipole_l': [], 'dipole_m': [], 'dipole_h': []}    
        for config, value in raw_dict.items():
            dipole_l=value.get('dipole_l')
            dipole_m=value.get('dipole_m')
            dipole_h=value.get('dipole_h')
            if not any(np.isnan(x) for x in [dipole_l,dipole_m,dipole_h]):
                data_dict['config'].append(config)
                data_dict['dipole_l'].append(dipole_l)
                data_dict['dipole_m'].append(dipole_m)
                data_dict['dipole_h'].append(dipole_h)

        df = pd.DataFrame.from_dict(data_dict)
        return df           
################################################################################
################################################################################
    def init1(self,Gaus='Vacuum', workdir='./'):
        rundir = workdir+'opt-test_Vacuum/'    
        
        filelist = glob.glob(rundir + f'nvt_vacuum.out')
    
        df4 = self.get_multipole_statistics(filelist)
        df4['muG'] = df4.mean()
        df4.to_csv('Dipole_Vacuum.csv', index=False)        
        mu_Vacuum=df4['muG']
        return mu_Vacuum
################################################################################
    def init2(self,exp_diconst,Gaus='PCM1', workdir='./'):
        rundir = workdir+'opt-test_PCM1/'
        
        filelist = glob.glob(rundir + f'PCM1.out')
          
        df2 = self.get_multipole_statistics(filelist)
        df2['muL_PCM1'] = df2.mean()
        df2.to_csv('Dipole_PCM1.csv', index=False)
        PCM1= df2['muL_PCM1']
        return PCM1
################################################################################
    def init3(self,cal_diconst,Gaus='PCM2', workdir='./'):
        rundir = workdir+'opt-test_PCM2/'
        
        filelist = glob.glob(rundir + f'PCM2.out')
        df3 = self.get_multipole_statistics(filelist)
   
        df3['muL_PCM2'] = df3.mean()
        df3.to_csv('Dipole_PCM2.csv', index=False)
        PCM2=df3['muL_PCM2']
        
        return PCM2       
################################################################################
    def init4(self,Molecule,Gaus='SCEE_V1', workdir='./'):  
        filelist=[]
        for Q in range(1,4):
            filelist.append(glob.glob(workdir+rundir + f'*_c*_q{Q}_v1.out'))
        df1 = Molecule.get_multipole_statistics(filelist)
        print(df1.describe())
        print(df1.head())
        muL_list = []
        for index, row in df1.iterrows():
             x = np.array([self.num1, self.num2, self.num3])
             y = np.array([row['dipole_l'], row['dipole_m'], row['dipole_h']])
             coeff = np.polyfit(x, y, 2)
             b = coeff[1] - 1
             fact = b*b-4*coeff[0]*coeff[2]
             if (coeff[0] >= 1.0e-6) and (fact > 0):
                  muL = (-b - np.sqrt(fact)) / (2*coeff[0])
             else:
                 coeff = np.polyfit(x, y, 1)
                 b = coeff[0] - 1
                 muL = - coeff[1]/b
             muL_list.append(muL)
        df1['muL_SCEE'] = muL_list
        df1.to_csv('Dipole_SCEE.csv', index=False)
        SCEE=df1['muL_SCEE']
        return SCEE
################################################################################
