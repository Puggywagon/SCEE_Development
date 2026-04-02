import subprocess
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
import glob
import csv
################################################################################
class Dipole_distributions(object):
    def __init__(self):
        pass
################################################################################
    def get_dipole(self,i):
        df = pd.read_csv(f'./Dipole.csv')
        df.columns = ['config', 'dipole_l', 'dipole_m', 'dipole_h','muL']
        dipole_moments=df['muL']
        dipole=np.mean(dipole_moments)
        std_dipole=np.std(dipole_moments)
    
        data = {
        'Replica': i,
        'Dipole Moment': dipole,
        'stdev dipole': std_dipole
        },
        df=pd.DataFrame(data)

        with open(f"./Results/Dipoles.csv", "a", newline="") as file:  # Open in append mode
            df.to_csv(file, header=False, index=False)
            return df
###############################################################################################
