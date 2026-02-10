import gromacs
import numpy as np
import matplotlib.pylab as plt
###############################################################
def plt_coordnumb(T,p,color):
    with open(f'./5ns/{T}K/{p}Bar/volume.txt') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Volume"):
                 volume = float(line.split()[1])
    
    DIR = f'./5ns/{T}K/{p}Bar/'
    xvgfile = DIR + 'OO_RDF.xvg'

    data = np.loadtxt(xvgfile, skiprows=25)

    N = 0.0
    r_list = data[:,0]
    N_list = [N]
    for i in range(len(data)-1):
        r0 = data[i  ,0]
        r1 = data[i+1,0]
        g0 = data[i  ,1]
        g1 = data[i+1,1]
        dr = r1 - r0
        num_density=900/volume
        f0 = 4.0*np.pi*r0*r0*g0*num_density
        f1 = 4.0*np.pi*r1*r1*g1*num_density
        N += 0.5*dr*(f0+f1)
        N_list.append(N)

    plt.plot(r_list, N_list,color=color, label=f'p={p:.0f} bar')
    plt.xlim([0.2, 0.5])
    plt.ylim([0.0, 4])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)  
    plt.xlabel('r / nm',fontsize=22)
    plt.ylabel('$N_{c}$',fontsize=22)
    #plt.title('{}K'.format(T),fontsize=14,weight='bold')
################################################################################
def plt_close(T):
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': 16})
    plt.savefig(f'./5ns/Analysis/RDFs/Coord_Numbs/Figs/{T}_Coords_noleg.pdf', dpi=300)
    plt.savefig(f'./5ns/Analysis/RDFs/Coord_Numbs/Figs/{T}_Coords_noleg.png', dpi=300)
    #plt.show()
    plt.close()    
##############################################################################################
def pressure(T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(P_list))]
    plt.figure(figsize=(24, 20))
    for p, color in zip(P_list, color_list):    
        plt_coordnumb(T,p,color)
##############################################################################################
def temperature():
    T_list=[298] #Temp in K
    for temp in range(400,1001,100):
        T_list.append(temp)
    for T in T_list:
        pressure(T)
        plt_close(T)
##############################################################################################
temperature()    
