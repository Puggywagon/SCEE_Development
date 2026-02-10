import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
######################################################
def test(T_list1,T_list3,P_list):
    df=pd.read_csv('Di_Const.csv')
    df.columns=['Replica','Temperature','Pressure','Diconstant']
    means=df.groupby(['Temperature','Pressure'])['Diconstant'].mean()
    stds=df.groupby(['Temperature','Pressure'])['Diconstant'].std()
    SE=((df.groupby(['Temperature','Pressure'])['Diconstant'].std())/(math.sqrt(7))*2)
    
    print(means)
    print(stds)
    print(SE)
    
    Temperature_list=[]
    Pressure_list=[]
    for T_list in [T_list1,T_list3]:    
        for T in T_list:
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            Temperature_list.append(T)
            for P in P_list:
    	        Pressure_list.append(P)
    Temperature_list.sort()
    data={
    'Temperature':Temperature_list,
    'Pressure':Pressure_list,
    'Mean': means,
    'SE':SE
    }
    df2=pd.DataFrame(data)
    with open(f"Diconsts2.csv", "a", newline="") as file:  # Open in append mode
        df2.to_csv(file, header=False, index=False)
        return df2
######################################################
def read_data():
    df3 = pd.read_csv('Diconsts2.csv')
    df3.columns = ['Temperature', 'Pressure', 'Epsilon','SE']
    return df3
##########################################################################################
def plt_epsilon(df2,T_list):    
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(14, 6))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        minval=tmp['SE']
        maxval=tmp['SE']
        yerr=[minval,maxval]
        plt.errorbar(tmp['Pressure'], tmp['Epsilon'], yerr=yerr, capsize=3, ecolor = "black", ls='None', color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'Dielectric Constant',fontsize=12)
    plt.xlabel(r'Pressure / bar',fontsize=12)
    plt.xscale('log')
    plt.xlim([0.5, 1.0e4])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Deielectric Constant vs Pressure',fontsize=14,weight='bold')
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'Epsilon_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'Epsilon_{Min}_{Max}.png', dpi=300)
    plt.close()
##########################################################################################
T_list1=[298]
for temp in range(400,1001,100):
    T_list1.append(temp)
T_list2=[298]
for temp in range(273,374,10):
    T_list2.append(temp)
T_list3=[]
for temp in range(273,374,10):
    T_list3.append(temp)
p_ref = [1,2,5]
P_list = []
for scale in np.logspace(0, 3, 4):
    P_list += [p*scale for p in p_ref]
test(T_list1,T_list3,P_list)
for T_list in [T_list1,T_list2]:    
    df3=read_data()
    plt_epsilon(df3,T_list)
