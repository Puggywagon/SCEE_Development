import pandas as pd
class Density_Dipoles(object):
    def __init__(self):
        pass
#################################################################################
    def read_Density(self):
        df2=pd.read_csv('./Results/Average_Dipoles.csv')
        df2.columns = ['Temperature', 'Pressure', 'Dipoles','SE','std means']
        df2_dipoles = df3[['My Dipole']]

        df3=pd.read_csv('./Results/Density.csv')
        df3.columns = ['Temperature', 'Pressure', 'Density']
    
        df4=pd.concat([df3, df2_dipoles], axis=1)
    
        with open(f"./Results/Density_Dipoles.csv", "a", newline="") as file:  # Open in append mode
            df4.to_csv(file, header=False, index=False)
            return df4    
#################################################################################
