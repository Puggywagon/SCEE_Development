import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
################################################################################
def enter_dirP(p):
    dire='{}Bar'.format(p)
    os.chdir(dire)
################################################################################
def enter_dirT(T):
    dire='{}K'.format(T)
    os.chdir(dire)
################################################################################
def enter_dirtime(Time):
    dire='{}ns'.format(Time)
    os.chdir(dire)
################################################################################
def exit_dir():
    dire=f'../'
    os.chdir(dire)
################################################################################
def read_rdf():
    data = np.loadtxt("OO_RDF.xvg", skiprows=26) 
    r =data[:, 0]
    gr =data[:, 1]
    return r,gr
################################################################################
def cood_noleg(T,p,df2,color):
    tmp=df2[df2['Pressure']==p]
    density=900
    r=[]
    coord_numb=[]
    for i in range(0,len(tmp['r']),1):
        rs=tmp['r'].iloc[i]
        r.append(rs)
        r2=rs**2
        minr3,maxr3=quad(integrand2, 0.2,0.4, args=r2)
        r3=maxr3-minr3
        gr=tmp['gr'].iloc[i]
        mingr2,maxgr2=quad(integrand, 0.2,0.4, args=gr)
        gr2=maxgr2-mingr2
        numb=((gr*r3)-(r2*gr2))
        cood_num=4*np.pi*density*numb
        coord_numb.append(cood_num)
    
       
    plt.plot(r,coord_numb, color=color, label=f'p={p:.0f} bar')
    #plt.xlim(0,1.4)
    #plt.ylim(0,5)
    plt.xlabel('r / nm',fontsize=12)
    plt.ylabel('$N_{c}$',fontsize=12)
    #plt.yscale('function')
    plt.title('{}K'.format(T),fontsize=14,weight='bold')
################################################################################    
def integrand(r,gr_r2):
    return gr_r2
################################################################################    
def integrand2(r,r2):
    return r2
################################################################################    
def integrand3(r,r3):
    return r3
################################################################################
def plt_close_noleg(T):
    plt.savefig(f'../Analysis/RDFs/Coord_Numbs/Figs/{T}_Coords_noleg.pdf', dpi=300)
    plt.savefig(f'../Analysis/RDFs/Coord_Numbs/Figs/{T}_Coords_noleg.png', dpi=300)
    plt.close()
################################################################################
def get_density():
    gromacs.environment.flags['capture_output'] = True
    input_str = 'Density 0\n'
    tmp = gromacs.tools.Energy(f='Water_QMMM_md3',
                               s='Water_QMMM_md3',
                               input=input_str)
    chk, density_output, stderr = tmp.run()

    f = open(f'Density.txt', 'w')
    f.write(density_output)
    f.close()

    density_lines = density_output.splitlines()
    density_line = density_lines[9]  # Assuming density is on the 5th line
    density = float(density_line.split()[1])

    return density
################################################################################
def cood_leg(T,p,df2,color):
    tmp=df2[df2['Pressure']==p]
    density=900
    r=[]
    coord_numb=[]
    for i in range(0,371,1):
        rs=tmp['r'].iloc[i]
        r.append(rs)
        r2=rs**2
        r3=(r**3)/3
        #minr3,maxr3=quad(integrand2, 0.2,0.4, args=r2)
        #minr4,maxr4=quad(integrand3, 0.2,0.4, args=r3)
        gr=tmp['gr'].iloc[i]
        #mingr2,maxgr2=quad(integrand, 0,20, args=gr)
        
        gr_r2=gr*r2
        mingr2_r3,maxgr2_r3=quad(integrand, 0.2,0.4, args=gr_r2)
        
        numb=(maxgr2_r3-mingr2_r3)
        cood_num=4*np.pi*density*numb
        coord_numb.append(cood_num)
    
    plt.plot(r,coord_numb, color=color, label=f'p={p:.0f} bar')
    #plt.xlim(0,1.4)
    #plt.ylim(0,5)
    plt.xlabel('r / nm',fontsize=12)
    plt.ylabel('$N_{c}$',fontsize=12)
    #plt.yscale('function')
    plt.title('{}K'.format(T),fontsize=14,weight='bold')
################################################################################
def plt_close_leg(T):
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'../Analysis/RDFs/Coord_Numbs/Figs/{T}_Coords_leg.pdf', dpi=300)
    plt.savefig(f'../Analysis/RDFs/Coord_Numbs/Figs/{T}_Coords_leg.png', dpi=300)
    plt.close()
###############################################################################################
def df_csv(T,p,density,i,j):
    data = {
        'Pressure': p,
        'Density': density,
        'r': i,
        'gr': j
    },
    df=pd.DataFrame(data)
    with open(f"../../Analysis/RDFs/Coord_Numbs/CSV/RDF_Density{T}.csv", "a", newline="") as file:  # Open in append mode
        df.to_csv(file, header=False, index=False)
        return df
#################################################################################
def read_csv(T):
    df2 = pd.read_csv(f'../../Analysis/RDFs/Coord_Numbs/CSV/RDF_Density{T}.csv')
    df2.columns = ['Pressure', 'Density','r','gr']
    return df2
##############################################################################################
def pressure_noleg(T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(P_list))]
    plt.figure(figsize=(14, 6))
    for p,color in zip(P_list,color_list):
        enter_dirP(p)
        #density=get_density()
        #r,gr=read_rdf()
        #for i,j in zip(r,gr):
            #df_csv(T,p,density,i,j)
        df2=read_csv(T)
        cood_noleg(T,p,df2,color)
        exit_dir()
##############################################################################################
def pressure_leg(T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(P_list))]
    plt.figure(figsize=(14, 6))
    for p, color in zip(P_list, color_list):
        enter_dirP(p)
        df2=read_csv(T)
        cood_noleg(T,p,df2,color)
        exit_dir()
##############################################################################################
def temperature():
    T_list=[298] #Temp in K
    for temp in range(400,1001,100):
        T_list.append(temp)
    for T in T_list:
        enter_dirT(T)        
        pressure_noleg(T)
        plt_close_noleg(T)        
        pressure_leg(T)
        plt_close_leg(T)
        exit_dir()
################################################################
def time():
    time_length = [5]
    for Time in time_length:
        enter_dirtime(Time)
        temperature()
        exit_dir()
################################################################            
time()
                                
