import pandas as pd
class Neighbours_Dipoles(object):
    def __init__(self):
        pass
#################################################################################
    def read_Average_Dipoles(self):
        df2=pd.read_csv('./Results/Average_Dipoles.csv')
        df2.columns = ['Temperature', 'Pressure', 'Dipoles','SE','std means']
        df2_dipoles = df3[['My Dipole']]

        df3=pd.read_csv('./Results/Neighbours.csv')
        df3.columns = ['Temperature', 'Pressure','Number of H-Bonds','Percent of H-Bonds','Nearest']
    
        df4=pd.concat([df3, df2_dipoles], axis=1)
    
        with open(f"./Results/Neighbours_Dipoles.csv", "a", newline="") as file:  # Open in append mode
            df4.to_csv(file, header=False, index=False)
            return df4    
#################################################################################
