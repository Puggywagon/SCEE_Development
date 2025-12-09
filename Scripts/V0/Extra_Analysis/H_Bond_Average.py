import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
################################################################################
class H_Bond_Average(object):
    def __init__(self):
        pass
######################################################
    def H_bond_Average(self):
        df=pd.read_csv("./Results/Neighbouring.csv")
        df.columns=['Replica','Number of H-Bonds','Percent of H-Bonds','Nearest Neighbour']
        
        means=df['Nearest Neighbour'].mean()
        stds=df['Nearest Neighbour'].std()
        SE=(((df['Nearest Neighbour'].std())/(math.sqrt(6)))*2)
        
        counting=df['Replica'].count()
           
        data={
        'Mean Neighbours': means,
        'SE Neighbours':SE,
        'counts': counting
        }
        df3=pd.DataFrame(data)
        with open(f"./Results/Neighbours.csv", "a", newline="") as file:  # Open in append mode
            df3.to_csv(file, header=False, index=False)
            return df3
######################################################
