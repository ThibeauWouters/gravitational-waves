import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

from matplotlib import rc_context
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

FFT = False
FFT = True
COWLING = True
COWLING = False

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
ax1.set_ylabel('$\\ddot I_{zz} - \\ddot I_{RR}$')
ax1.set_xlabel('$t$ $\\rm{(s)}$')
ax1.plot(df_TD['t'] , df_TD['hp'], linestyle = '-' ,linewidth = 1.0 )
#ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlim(t_min,t_max)

if (FFT):
  #FFT part
  df_FD = pd.read_csv("GW_FD.csv", delimiter = '\t', usecols = ['f','fft','psd'])
  df_FD = df_FD[['f','fft','psd']]
  ax2.plot(df_FD['f']/1000, df_FD['psd'] ,linestyle = '-' ,linewidth = 1.0, label = '$\\rm{FFT}$' )

  if (COWLING):
    ax2.axvline(x=2.328, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=4.300, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=6.320, color='g',label = '$H_2$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=8.153, color='b',label = '$H_3$', linestyle='--',linewidth = 0.6)
  else:
    ax2.axvline(x=1.168, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=4.030, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=1.685, color='g',label = '$^2f$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=1.004, color='b',label = '$i_{2}$', linestyle='--',linewidth = 0.6)
  
  ax2.legend(loc='center right')

ax2.grid(True)
ax2.set_ylabel('$\\rm{PSD}$')
ax2.set_xlabel('$f$ $\\rm{[kHz]}$')
#ax2.set_yscale('log')
ax2.set_xlim(0.0,8.5)


fig.tight_layout()
fig.subplots_adjust(hspace=0.35)
fig.savefig("./GW_combine.png", bbox_inches="tight")
