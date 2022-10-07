import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

from matplotlib import rc_context
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':15})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


FFT = False
FFT = True
COWLING = True
COWLING = False

filepath = "../"
savepath = './'

simname = [ \
        "output.log" ,\
        ]

plot_color = ["black", "r","g","b","fuchsia"]
ls = ["-","-.", "--", ":","-."]
lw = [1,1.5,1.5,1.5,1.5]

fig = plt.figure(1)
fig, (ax1,ax2) = plt.subplots(2, 1)
plt.subplots_adjust(bottom=0.13, left=0.14, right=0.96, top=0.96, hspace=0.45)
plt.subplots_adjust(right=0.96, top=0.96, hspace=0.40)

df_TD = pd.read_csv("GW_TD_xy.csv", delimiter = '\t', usecols = ['t','I11_dot','I11_ddot','I22_dot','I22_ddot','I12_ddot','hp','hc'])
df_TD = df_TD[['t','I11_dot','I11_ddot','I22_dot','I22_ddot','I12_ddot','hp','hc']]
time=1000*df_TD['t']
t_min = time[0]
t_max = time[len(time)-1]
ax1.set_ylabel('$h$ ($d$ = 100 Mpc)')
ax1.set_xlabel('$t$ $\\rm{[ms]}$')
ax1.plot(1000*df_TD['t'] , df_TD['hp'], linestyle = '-' ,linewidth = 1.0 ,label='$h_+$',color='g')

ax1.plot(1000*df_TD['t'] , df_TD['hc'], linestyle = '-' ,linewidth = 1.0 ,label='$h_\\times$',color='r')
ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlim(t_min,t_max)
#ax1.set_ylim(-3,0)

ax2.set_ylabel('$\\rm{FFT}$')
ax2.set_xlabel('$f$ $\\rm{[kHz]}$')
ax2.set_xlim(0.0,6.0)
#ax2.set_yscale('log')
#ax2.set_ylim(1E-31,1E-27)

if (FFT):
  #FFT part
  df_FD = pd.read_csv("GW_FD_xy.csv", delimiter = '\t', usecols = ['f','fft_hp','psd_hp','fft_hc','psd_hc'])
  df_FD = df_FD[['f','fft_hp','psd_hp','fft_hc','psd_hc']]
  ax2.plot(0.001*df_FD['f'], df_FD['fft_hp'], linestyle = '-' ,linewidth = 1.0, label = '$\\rm{FFT}(h_+)$',color='g' )
  ax2.plot(0.001*df_FD['f'], df_FD['fft_hc'], linestyle = '-' ,linewidth = 1.0, label = '$\\rm{FFT}(h_\\times)$',color='r' )
  #ax2.plot(0.001*df_FD['f'], df_FD['psd_hp'], linestyle = '-' ,linewidth = 1.0, label = '$\\rm{PSD}(h_+)$',color='g' )
  #ax2.plot(0.001*df_FD['f'], df_FD['psd_hc'], linestyle = '-' ,linewidth = 1.0, label = '$\\rm{PSD}(h_\\times)$',color='r' )
  if (COWLING):
    #ax2.axvline(x=2.328, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=4.300, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=6.320, color='g',label = '$H_2$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=8.153, color='b',label = '$H_3$', linestyle='--',linewidth = 0.6)
  else:
    #ax2.axvline(x=1.168, color='black',label = '$F$', linestyle='--',linewidth = 1.0)
    #ax2.axvline(x=4.030, color='r',label = '$H_1$', linestyle='--',linewidth = 1.0)
    #ax2.axvline(x=1.685, color='black',label = '$^2f$', linestyle='--',linewidth = 1.0)
    #ax2.axvline(x=2.512, label = '$^2p$', linestyle='--',linewidth = 1.0)
    #ax2.axvline(x=4.030-2.512, label = '$H_1-^2p$', linestyle='--',linewidth = 1.0)
    #ax2.axvline(x=1.004, color='b',label = '$i_{2}$', linestyle='--',linewidth = 1.0)

  #ax2.set_yscale('log')
  ax2.legend(loc='best')
  ax2.grid(True)


fig.savefig(savepath + 'fig_GW_xy.png', bbox_inches="tight")
