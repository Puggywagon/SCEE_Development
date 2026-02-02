import subprocess
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
import glob
import csv
################################################################################
class Hbond_Analysis_Average(object):
    def __init__(self):
        pass
######################################################
    def Hbond_Analysis(self,T,P):
        csvs=glob.glob(f'Replica_*/{T}K/{P}.0Bar/hbanalysis.csv')
        data_1={'donors': [], 'acceptors': [], 'nHB': []}
        for csv in csvs:
            df=pd.read_csv(csv)
            df.columns=['config', 'dipole_l', 'dipole_m', 'dipole_h','muL','donors','acceptors','nHB']
            data_1['donors']=df['donors']
            data_1['acceptors']=df['acceptors']
            data_1['nHB']=df['nHB']
        
        data_1=pd.DataFrame(data_1)
        
        mdonors=data_1['donors'].mean()
        macceptors=data_1['acceptors'].mean()
        mnHB=data_1['nHB'].mean()
        
        return mdonors, macceptors, mnHB      
######################################################
    def Cluster_Analysis(self,data):
        columns=['Time', 'Temperature', 'Pressure', 'Density','Epsilon','Number of HBonds','Percentage of H-bonds','Nearest Neighbour','Dipole Moment','State']
        df = pd.read_csv("../TIP4P_Results/5ns/Analysis/data.csv",header=None, names=columns)
        print(df)
        df['donors']=data['donors']
        df['acceptors']=data['acceptors']
        df['nHB']=data['nHB']
        clusters={'clusters': []}
        for i in range(0,len(df)):
            nhb=df['nHB']
            ner_neh=df['Nearest Neighbour']
            cluster=nhb[i]/ner_neh[i]
            clusters['clusters'].append(cluster)
        df['Clusters']=pd.DataFrame(clusters)    
        
        print(df)
        
        df.to_csv("../TIP4P_Results/5ns/Analysis/data_2.csv", index=False,header=True)
######################################################

Hbond_go=Hbond_Analysis_Average()

Temps=[298,400,500,600,700,800,900,1000]
Press=[1,2,5,10,20,50,100,200,500,1000,2000,5000]

data={'donors':[], 'acceptors':[], 'nHB':[]}       
for T in Temps:
    for P in Press:
        mdonors, macceptors, mnHB=Hbond_go.Hbond_Analysis(T,P)
        data['donors'].append(mdonors)
        data['acceptors'].append(macceptors)
        data['nHB'].append(mnHB)
data=pd.DataFrame(data)        
Hbond_go.Cluster_Analysis(data)
