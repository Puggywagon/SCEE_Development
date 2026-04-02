import gromacs
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
################################################################################
class Di_const(object):
    def __init__(self):
        pass
################################################################################
    def get_epsilon(self,system_title,i,T,p):
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        tmp = gromacs.tools.Dipoles(f=f'{system_title}_QMMM_md3',
                                s=f'{system_title}_QMMM_md3',
                                input=input_str)
        chk, dipole_output, stderr = tmp.run()

        f = open(f'Dipole.txt', 'w')
        f.write(dipole_output)
        f.close()    

        dipole_lines = dipole_output.splitlines()  # Split text into lines
        epsilon_line = dipole_lines[-1]  # Assuming epsilon is on the last line
        epsilon = float(epsilon_line.split()[2])  # Extract the numerical value

        data = {
            'Replica': i,
            'Epsilon': epsilon
        },
        df=pd.DataFrame(data)
    
        with open(f"./Results/Di_Const.csv", "a", newline="") as file:  # Open in append mode
            df.to_csv(file, header=False, index=False)
            return df
######################################################
