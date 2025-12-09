import gromacs
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
################################################################################
class Di_const_Average(object):
    def __init__(self):
        pass
######################################################
    def Di_Const_Average(self):
        df=pd.read_csv('./Results/Di_Const.csv')
        df.columns=['Replica','Diconstant']
        means=df['Diconstant'].mean()
        stds=df['Diconstant'].std()
        SE=(df['Diconstant'].std())/(math.sqrt(10))
        counting=df['Replica'].count
    
        data={
        'Mean': means,
        'SE':SE,
        'counts': counting
        }
        df2=pd.DataFrame(data)
        with open(f"./Results/Average_Di_Constants.csv", "a", newline="") as file:  # Open in append mode
            df2.to_csv(file, header=False, index=False)
            return df2
######################################################
