import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as klp


import numpy as np
from scipy.integrate import odeint
import math as mt
import scipy.optimize as opt
import pylab

#------------------------------------------------------------------ Parameters
# kinetic parameters
k = {}
k['A'] = 45
k['B'] = 50

#----------------------------------------------------------- Control variables
# Running time (min)
t_f = 250
# Initial substrate titer (mΜ)
S1_initial = 250
# Initial intermediate titer (mM)
S2_initial = 0
# Initial cofactor titer (mΜ)
S4_initial = 4
# Enzyme titers (μΜ)
E_A = 1
E_B = 1

#-----------------------------------------------------------Algebraic equations

def Reaction_rates(S,t):  
    v = {}
    v['I'] = k['A']*E_A*S[4]*S[1]/ (S[1]+15)/ (S[4]+25)
    v['II'] = k['B']*E_B*S[2]*S[5]/ (S[2]+30)/ (S[5]+10)
    return(v)

#-------------------------------------------------------Differential equations
def Material_balances(S,t):         
    v = Reaction_rates(S,t)
    dS1dt = -v['I']
    dS2dt = +v['I']-v['II']
    dS3dt = +v['II']
    dS4dt = -v['I']+v['II']
    dS5dt = +v['I']-v['II']
    dSdt = [0, dS1dt, dS2dt, dS3dt, dS4dt, dS5dt]
    return(dSdt)

#--------------------------------------------------------------------Solver
t = np.linspace(0,t_f,200)  
S0 = [0, S1_initial, S2_initial, 0, S4_initial, 0]
S = odeint(Material_balances,S0,t)

import matplotlib.patches as patches


import matplotlib.pyplot as plt
import matplotlib as mpl
fig = plt.figure(figsize=(9,9))

plt.subplot(211)

ax = fig.gca()
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.major.width'] = 3
ax.plot(t,S[:,1],'-', color = 'black', linewidth=2, label='$\it{S}$$_1$')
ax.plot(t,S[:,2],'-', color = 'green', linewidth=2.5, label='$\it{S}$$_2$')
ax.plot(t,S[:,3],'-', color = 'grey', linewidth=2, label='$\it{S}$$_3$')
ax.plot(t,S[:,4],'-', color = 'blue', linewidth=2, label='$\it{c}$$_1$')
ax.plot(t,S[:,5],'--', color = 'red', linewidth=2, label='$\it{c}$$_2$')
#ax.plot(t,c1,'-.', color = 'black', linewidth=2, label='$\it{c}$$_2$')
plt.ylabel('$\it{c}$$_1$ / (mM)', fontsize=25)


#rect = patches.Rectangle((-1, -3), 53.5, 15, linestyle= '-.', linewidth=1,  edgecolor='black', facecolor='none')
#ax.add_patch(rect)

plt.ylabel('$\it{S}$$_i$ or $\it{c}$$_i$ / (mM)', fontsize=25)
plt.xlabel('$\it{t}$ / (min)', fontsize=25)

ax.grid(False)
#ax.legend(fontsize=25)
#plt.legend(bbox_to_anchor=(1.05, 1), fontsize=25, loc='upper left')
#plt.ylim(0,8)


plt.xticks(fontsize=25)
#plt.xticks([0, 20, 40, 60, 80, 100])
plt.yticks(fontsize=25)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
plt.savefig('simu crn a.png', bbox_inches='tight')
plt.savefig('simu crn a.pdf', bbox_inches='tight')
#plt.show()


plt.axvline(x=5, linewidth=2.5, color = 'black', linestyle = '--')
plt.axvline(x=200, linewidth=2.5, color = 'black', linestyle = '--')


plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)

plt.subplot(212)
ax = fig.gca()
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.major.width'] = 3
#ax.plot(t,S[:,1],'-', color = 'black', linewidth=2, label='$\it{S}$$_1$')
ax.plot(t,S[:,2],'-', color = 'green', linewidth=2.5, label='$\it{S}$$_2$')
#ax.plot(t,S[:,3],'-', color = 'grey', linewidth=2, label='$\it{S}$$_3$')
ax.plot(t,S[:,4],'-', color = 'blue', linewidth=2, label='$\it{c}$$_1$')
ax.plot(t,S[:,5],'--', color = 'red', linewidth=2, label='$\it{c}$$_2$')
#ax.plot(t,c1,'-.', color = 'black', linewidth=2, label='$\it{c}$$_2$')
plt.ylabel('$\it{c}$$_1$ / (mM)', fontsize=25)

plt.axvline(x=5, linewidth=2.5, color = 'black', linestyle = '--')
plt.axvline(x=200, linewidth=2.5, color = 'black', linestyle = '--')

#rect = patches.Rectangle((-1, -3), 53.5, 15, linestyle= '-.', linewidth=1,  edgecolor='black', facecolor='none')
#ax.add_patch(rect)

plt.ylabel('$\it{S}$$_i$ or $\it{c}$$_i$ / (mM)', fontsize=25)
plt.xlabel('$\it{t}$ / (min)', fontsize=25)

ax.grid(False)
#ax.legend(fontsize=25)
#plt.legend(bbox_to_anchor=(1.05, 1), fontsize=25, loc='upper left')
#plt.ylim(0,8)


plt.xticks(fontsize=25)
plt.yticks([0, 1, 2, 3, 4])
plt.yticks(fontsize=25)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
plt.savefig('First simu.pdf', bbox_inches='tight')
#plt.show()

