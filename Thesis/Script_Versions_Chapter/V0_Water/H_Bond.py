import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
################################################################################
class H_Bond(object):
    def __init__(self):
        pass
########################################################################################################################
    def get_Hbonds(self,system_title,i,T,p):
        gromacs.environment.flags['capture_output'] = True
        input_str = '0 0\n'
        tmp = gromacs.tools.Hbond(f=f'{system_title}_QMMM_md3',
                                  s=f'{system_title}_QMMM_md3',
                                  input=input_str)
        chk, hbond_output, stderr = tmp.run()
        f = open(f'HBond.txt', 'w')
        f.write(hbond_output)
        f.close()

        lines = hbond_output.splitlines()
        average_hbonds_line = lines[-1]  # Assuming average is on the second-to-last line
        hbond = float(average_hbonds_line.split()[-5])  # Extract the numerical value
        total_hbond = 1800
        percent_hbond = (hbond / total_hbond) * 100
        Nearest_neighbour=(hbond / total_hbond)*4

        data = {
            'Replica': i,
            'Number of H-Bonds': hbond,
            'Percent of H-Bonds': percent_hbond,
            'Nearest Neighbour': Nearest_neighbour
        },
        df=pd.DataFrame(data)
    
        with open(f"./Results/Neighbouring.csv", "a", newline="") as file:  # Open in append mode
            df.to_csv(file, header=False, index=False)
            return df
###############################################################################################
