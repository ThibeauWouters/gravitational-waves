import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':10})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

FFT = False
FFT = True
COWLING = True
COWLING = False

filepath = "../"

simname = [ \
        "output.log" ,\
        ]

plot_color = ["black", "r","g","b","fuchsia"]
ls = ["-","-.", "--", ":","-."]
lw = [1,1.5,1.5,1.5,1.5]

fig = plt.figure(1)
fig, (ax1,ax2) = plt.subplots(2, 1)
#plt.subplots_adjust(bottom=0.13, left=0.14, right=0.96, top=0.96, hspace=0.45)
#plt.subplots_adjust(right=0.96, top=0.96, hspace=0.40)

df_TD = pd.read_csv("GW_TD.csv", delimiter = '\t', usecols = ['t','I11_dot','I11_ddot','I22_dot','I22_ddot','hp'])
df_TD = df_TD[['t','I11_dot','I11_ddot','I22_dot','I22_ddot','hp']]
time=df_TD['t']
t_min = time[0]
t_max = time[len(time)-1]
ax1.set_ylabel('$\\ddot I_{11} - 2\\ddot I_{22}$')
ax1.set_xlabel('$t$ $\\rm{(s)}$')
ax1.plot(df_TD['t'] , df_TD['hp'], linestyle = '-' ,linewidth = 1.0 )
#ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlim(t_min,t_max)
#ax1.set_ylim(-3,0)

ax2.set_ylabel('$\\rm{PSD}$')
ax2.set_xlabel('$f$ $\\rm{(Hz)}$')
ax2.set_xlim(00,8500)
#ax2.set_ylim(-1E19,4E20)

if (FFT):
  #FFT part
  df_FD = pd.read_csv("GW_FD.csv", delimiter = '\t', usecols = ['f','fft','psd'])
  df_FD = df_FD[['f','fft','psd']]
  ax2.plot(df_FD['f'], df_FD['psd'] ,linestyle = '-' ,linewidth = 1.0, label = '$\\rm{FFT}$' )

  if (COWLING):
    ax2.axvline(x=2328, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=4300, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=6320, color='g',label = '$H_2$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=8153, color='b',label = '$H_3$', linestyle='--',linewidth = 0.6)
  else:
    ax2.axvline(x=1168, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=4030, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=1685, color='g',label = '$^2f$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=1004, color='b',label = '$i_{2}$', linestyle='--',linewidth = 0.6)
  
  ax2.legend(loc='best')
  ax2.grid(True)


fig.savefig("./GW_combine.png", bbox_inches="tight")
