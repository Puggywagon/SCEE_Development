#~/usr/bin/python3


import numpy as np
from statistics import mean
import math 
import matplotlib.pyplot as plt
import pandas as pd
import glob

import Molecule
import Particle
import SimBox


################################################################################
################################################################################
################################################################################
def read_conf(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    
    lines.pop(0)
    nmol = lines.pop(0)
    natoms = 4
#    print(nmol)
#    print(lines[-1])
    tmp = lines.pop(-1)
    data = tmp.split()
    x = float(data[0])
    y = float(data[1])
    z = float(data[2])
    box = SimBox.SimBox(x, y, z)
#    print(tmp)

    molecule_dict = {}
    for line in lines:
        data = line.split()
        m_name = data[0]
        if (m_name not in molecule_dict):
            m = Molecule.Molecule(m_name)
            m.box = box
            molecule_dict[m_name] = m
            
        p_name = data[1]
        p = Particle.Particle(p_name)
        x = float(data[3])
        y = float(data[4])
        z = float(data[5])
        p.position = np.array([x, y, z])
        
        molecule_dict[m_name].add_particle(p)

    for name, m in molecule_dict.items():
        m.add_acceptor(m.particle_list[0])
        m.add_donor(m.particle_list[1], m.particle_list[0])
        m.add_donor(m.particle_list[2], m.particle_list[0])

    return molecule_dict


################################################################################
################################################################################
################################################################################
#df = pd.read_csv('./data/replica_1/5ns/298K/1.0Bar/junk_dipoles.csv')
df = pd.read_csv('Simulations/replica_8/298K/1.0Bar/Dipole.csv')
df.sort_values(by='config', inplace=True)
print(df.head())


file_list = glob.glob('Simulations/replica_8/298K/1.0Bar/*_*.gro')

#filename = 'conf_1.gro'
data_dict = {'config': [], 'donors': [], 'acceptors': [], 'nHB': []}
for filename in file_list:
    tmp = filename.split('_')[-1]
    config = tmp.replace('.gro', '')
    molecule_dict = read_conf(filename)
    
    molecule_list = list(molecule_dict.values())
    m0 = molecule_list.pop(0)
    
    donors = 0
    acceptors = 0
    for m in molecule_list:
        if (m0.is_acceptor(m)):
            acceptors += 1
        if (m0.is_donor(m)):
            #print(m.ID)
            donors += 1
    
    nHB = acceptors + donors        
    print(f'{config}:  nHB = {nHB}, donors={donors}, acceptors={acceptors}')
    data_dict['config'].append(config)
    data_dict['donors'].append(donors)
    data_dict['acceptors'].append(acceptors)
    data_dict['nHB'].append(nHB)

tmp = pd.DataFrame.from_dict(data_dict)
tmp.sort_values(by='config', inplace=True)

df['donors'] = tmp['donors']
df['acceptors'] = tmp['acceptors']
df['nHB'] = tmp['nHB']
    
df.to_csv('hbanalysis.csv', index=False)


################################################################################
################################################################################
################################################################################

