import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
################################################################################
def plot_experiment():
    df = pd.read_csv('experimental_results.csv')
    df.columns = ['Temperature (K)', 'Pressure (bar)', 'Phase']
    Tc = 647.096      # Critical temperature (Tc) 647.096 K
    pc = 220.640      # Critical pressure (Pc)    220.640 bar
    rhoc = 17.873728  # Critical density (Dc)     17.873728 mol/l

    tmp = df[df['Phase']=='liquid']
    x_list = list(tmp['Temperature (K)']) + [Tc]
    y_list = list(tmp['Pressure (bar)']) + [pc]
    plt.plot(x_list, y_list, ls='solid', color='black')
    plt.plot([Tc], [pc], ls='-',marker='o', color='black',label='Experiment')

    plt.xlim([250, 1101])
    plt.ylim([0.5, 1.0e4])
    plt.yscale('log')
    plt.ylabel('pressure / bar')
    plt.xlabel('temperature / K')
################################################################################
def plot_vega():
    df1 = pd.read_csv('vega_results.csv')
    df1.columns = ['T / K', 'p / bar', 'rhol','rhog']
      
    Tc = 640
    pc = 146
    rhoc = 0.31   # g cm^{-3} 
    
    x_list =  [Tc] + list(df1['T / K'])
    y_list =  [pc] + list(df1['p / bar'])
    plt.plot(x_list, y_list, color='black', ls='dashed')
    plt.plot([Tc], [pc], ls='--',marker='o', color='black',label='TIP4P/2005')
    
    plt.xlim([250, 1101])
    plt.ylim([0.5, 1.0e4])
    plt.yscale('log')
################################################################################
def liquid_boundary(T,p):
    if T==298 and p==1.0:
        return "liquid"
    elif T==298 and p==2.0:
        return "liquid"
    elif T==298 and p==5.0:
        return "liquid"
    elif T==298 and p==10.0:
        return "liquid"
    elif T==298 and p==20.0:
        return "liquid"
    elif T==298 and p==50.0:
        return "liquid"
    elif T==298 and p==100.0:
        return "liquid"
    elif T==298 and p==200.0:
        return "liquid"
    elif T==298 and p==500.0:
        return "liquid"
    elif T==298 and p==1000.0:
        return "liquid"
    elif T==298 and p==2000.0:
        return "liquid"
    elif T==298 and p==5000.0:
        return "liquid"
    elif T==400 and p==1.0:
        return "liquid"
    elif T==400 and p==2.0:
        return "liquid"
    elif T==400 and p==5.0:
        return "liquid"
    elif T==400 and p==10.0:
        return "liquid"
    elif T==400 and p==20.0:
        return "liquid"
    elif T==400 and p==50.0:
        return "liquid"
    elif T==400 and p==100.0:
        return "liquid"
    elif T==400 and p==200.0:
        return "liquid"
    elif T==400 and p==500.0:
        return "liquid"
    elif T==400 and p==1000.0:
        return "liquid"
    elif T==400 and p==2000.0:
        return "liquid"
    elif T==400 and p==5000.0:
        return "liquid"
    elif T==500 and p==20.0:
        return "liquid"
    elif T==500 and p==50.0:
        return "liquid"
    elif T==500 and p==100.0:
        return "liquid"
    elif T==500 and p==200.0:
        return "liquid"
    elif T==500 and p==500.0:
        return "liquid"
    elif T==500 and p==1000.0:
        return "liquid"
    elif T==500 and p==2000.0:
        return "liquid"
    elif T==500 and p==5000.0:
        return "liquid"
    elif T==600 and p==100.0:
        return "liquid"
    elif T==600 and p==200.0:
        return "liquid"
    elif T==600 and p==500.0:
        return "liquid"
    elif T==600 and p==1000.0:
        return "liquid"
    elif T==600 and p==2000.0:
        return "liquid"
    elif T==600 and p==5000.0:
        return "liquid"
    else:
        return 0
################################################################################
def gas_boundary(T,p):
    if T==500 and p==1.0:
        return "gas"
    elif T==500 and p==2.0:
        return "gas"
    elif T==500 and p==5.0:
        return "gas"
    elif T==500 and p==10.0:
        return "gas"
    elif T==600 and p==1.0:
        return "gas"
    elif T==600 and p==2.0:
        return "gas"
    elif T==600 and p==5.0:
        return "gas"
    elif T==600 and p==10.0:
        return "gas"
    elif T==600 and p==20.0:
        return "gas"
    elif T==600 and p==50.0:
        return "gas"
    elif T==700 and p==1.0:
        return "gas"
    elif T==700 and p==2.0:
        return "gas"
    elif T==700 and p==5.0:
        return "gas"
    elif T==700 and p==10.0:
        return "gas"
    elif T==700 and p==20.0:
        return "gas"
    elif T==700 and p==50.0:
        return "gas"
    elif T==700 and p==100.0:
        return "gas"
    elif T==800 and p==1.0:
        return "gas"
    elif T==800 and p==2.0:
        return "gas"
    elif T==800 and p==5.0:
        return "gas"
    elif T==800 and p==10.0:
        return "gas"
    elif T==800 and p==20.0:
        return "gas"
    elif T==800 and p==50.0:
        return "gas"
    elif T==800 and p==100.0:
        return "gas"
    elif T==900 and p==1.0:
        return "gas"
    elif T==900 and p==2.0:
        return "gas"
    elif T==900 and p==5.0:
        return "gas"
    elif T==900 and p==10.0:
        return "gas"
    elif T==900 and p==20.0:
        return "gas"
    elif T==900 and p==50.0:
        return "gas"
    elif T==900 and p==100.0:
        return "gas"
    elif T==1000 and p==1.0:
        return "gas"
    elif T==1000 and p==2.0:
        return "gas"
    elif T==1000 and p==5.0:
        return "gas"
    elif T==1000 and p==10.0:
        return "gas"
    elif T==1000 and p==20.0:
        return "gas"
    elif T==1000 and p==50.0:
        return "gas"
    elif T==1000 and p==100.0:
        return "gas"
    else:
        return 0
################################################################################
def get_state(T,p):
    liquid_condition="liquid"
    gas_condition="gas"
    cmap = plt.cm.rainbow
    T_list=[298] #Temp in K
    for temp in range(400,1001,100):
        T_list.append(temp)
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    for T, color in zip(T_list, color_list):
        T==T
        if liquid_condition==liquid_boundary(T,p):
            plt.plot(T, p, ls='None',color='blue', marker='o')
        elif gas_condition==gas_boundary(T,p):
            plt.plot(T, p, ls='None',color='green', marker='^')
        else:
            plt.plot(T, p, ls='None', color='red', marker='*')
#################################################################################
def plot_update():
    plt.plot([],[], ls='', color='blue', marker='o',label='Liquid')
    plt.plot([],[], ls='', color='green', marker='^',label='Gas')
    plt.plot([],[], ls='', color='red', marker='*',label='Supercritical')
                 
    plt.ylabel("Pressure / Bar",fontsize=24)
    plt.xlabel("Temperature / K",fontsize=24)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)  
    plt.xlim([250, 1101])
    plt.ylim([0.5, 1.0e4])
    plt.yscale('log')   
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': 16})
################################################################################
def plot_close():
    plt.savefig(f'./5ns/Analysis/Phase_Diagram.pdf', dpi=300)
    plt.savefig(f'./5ns/Analysis/Phase_Diagram.png', dpi=300)
    plt.close()
################################################################################
def pressure(T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    for p in P_list:
        state=get_state(T,p)
##############################################################################################
def temperature():
    T_list=[298] #Temp in K
    for temp in range(400,1001,100):
        T_list.append(temp)
    for T in T_list:
        pressure(T)
################################################################
def time():
    time_length =[5]

    for Time in time_length:
        plot_experiment()
        plot_vega()
        temperature()
        plot_update()
####################################################################################
plt.figure(figsize=(24, 20))
time()
plot_close()
