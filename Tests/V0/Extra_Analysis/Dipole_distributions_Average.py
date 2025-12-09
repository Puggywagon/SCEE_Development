import subprocess
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
import glob
import csv
################################################################################
class Dipole_distributions_Average(object):
    def __init__(self):
        pass
######################################################
    def Dipole_Average(self):
        df=pd.read_csv("./Results/Dipoles.csv")
        df.columns=['Replica','Dipole','std dipole']
        means=df['Dipole'].mean()
        stds=df['Dipole'].std()
        SE=(((df['Dipole'].std())/(math.sqrt(6)))*2)
        std_means=df['std dipole'].mean()
        counting=df['Replica'].count()

        data={
        'Mean': means,
        'SE':SE,
        'std means': std_means,
        'counts': counting
        }
        df=pd.DataFrame(data)
        with open(f"./Results/Average_Dipoles.csv", "a", newline="") as file:  # Open in append mode
            df.to_csv(file, header=False, index=False)
            return df
######################################################
