#~/usr/bin/python3


import numpy as np
from statistics import mean
import math 
import matplotlib.pyplot as plt
import pandas as pd
import glob
from pathlib import Path

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
Replicas=[0,1,2,3,4,5,6]
Temps=[298,400,500,600,700,800,900,1000]
Press=[1,2,5,10,20,50,100,200,500,1000,2000,5000]

missing=[]
read_fail=[]
base=Path("Replica_Dipoles")

for i in Replicas:
    if i!=6:
        for T in Temps:

            tdir = base / f"Replica_{i+1}" / f"{T}K"
            if not tdir.is_dir():
                # whole temperature folder missing in that replica
                missing.append((i, T, "Tdir missing"))
                continue
            for P in Press:
                print(f'T:{T},P:{P}')
                candidates = sorted(tdir.glob(f"{P}*Bar"))

                if not candidates:
                    missing.append((i, T, P))
                    continue
                pdir = candidates[0]
                csv_path = pdir / "Dipole.csv"

                if not csv_path.is_file():
                    missing.append((i, T, P))
                    continue

                try:
                    df = pd.read_csv(csv_path)
                except Exception as e:
                   read_fail.append((i, T, P, type(e).__name__))
                   continue

                if "config" not in df.columns:
                    read_fail.append((i, T, P, "no config col"))
                    continue

                df.sort_values(by='config', inplace=True)
                print(df.head())


                file_list = glob.glob(f'Data/replica_{i}/{T}K/{P}.0Bar/*_*.gro')

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
                print(tmp)
                df['donors'] = tmp['donors']
                df['acceptors'] = tmp['acceptors']
                df['nHB'] = tmp['nHB']
    
                df.to_csv(pdir / "hbanalysis.csv", index=False)
################################################################################
################################################################################
################################################################################

