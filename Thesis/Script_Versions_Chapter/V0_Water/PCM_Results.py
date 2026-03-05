import pandas as pd
class PCM_Results(object):
    def __init__(self):
        pass
#################################################################################
    def read_PCM(self):
        df = pd.read_csv('./Results/pcm_dipoles.csv')
        df.columns = ['Temperature', 'Pressure', 'Density', 'Volume','Internal Energy','Enthalpy', 'Entropy','Cv','Cp','Sound','Joule','Viscosity','Therm','Phase','epsilon','g','PMC Dipole']
        
        df2=pd.read_csv('./Results/Average_Dipoles.csv')
        df2.columns = ['Temperature', 'Pressure', 'Dipoles','SE','std means']
        df2_dipoles = df3[['My Dipole']]

        df3=pd.read_csv('./Results/Neighbours.csv')
        df3.columns = ['Temperature', 'Pressure', 'Nearest']
        df3_neighbours = df4[['Nearest']]    
    
        df4=pd.concat([df, df2_dipoles, df3_neighbours], axis=1)
    
        with open(f"./Results/PCM_Dipoles.csv", "a", newline="") as file:  # Open in append mode
            df4.to_csv(file, header=False, index=False)
            return df4    
#################################################################################

