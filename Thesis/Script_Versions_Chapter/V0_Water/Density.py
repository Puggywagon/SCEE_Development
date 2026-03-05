import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
################################################################################
class Density(object):
    def __init__(self):
        pass
################################################################################
    def get_density(self,system_title,i,T,p):
        gromacs.environment.flags['capture_output'] = True
        input_str = 'Density 0\n'
        tmp = gromacs.tools.Energy(f=f'{system_title}_QMMM_md3',
                               s=f'{system_title}_QMMM_md3',
                               input=input_str)
        chk, density_output, stderr = tmp.run()

        f = open(f'Density.txt', 'w')
        f.write(density_output)
        f.close()

        density_lines = density_output.splitlines()
        density_line = density_lines[9]  # Assuming density is on the 5th line
        density = float(density_line.split()[1])

        data = {
            'Replica': i,
            'Density': density
        },
        df=pd.DataFrame(data)
        with open(f"./Results/Densities.csv", "a", newline="") as file:  # Open in append mode
            df.to_csv(file, header=False, index=False)
            return df
###############################################################################################
