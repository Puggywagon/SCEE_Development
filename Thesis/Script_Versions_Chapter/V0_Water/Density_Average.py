import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
################################################################################
class Density_Average(object):
    def __init__(self):
        pass
######################################################
    def Density_Average(self):
        df=pd.read_csv("./Results/Densities.csv")
        df.columns=['Replica','Density']
    
        means=df['Density'].mean()
        stds=df['Density'].std()
        SE=(((df['Density'].std())/(math.sqrt(6)))*2)
    
        counting=df['Replica'].count()
          
        data={
        'Mean': means,
        'SE':SE,
        'counts': counting
        }
        df3=pd.DataFrame(data)
        with open(f"./Results/Density.csv", "a", newline="") as file:  # Open in append mode
            df3.to_csv(file, header=False, index=False)
            return df3
######################################################
