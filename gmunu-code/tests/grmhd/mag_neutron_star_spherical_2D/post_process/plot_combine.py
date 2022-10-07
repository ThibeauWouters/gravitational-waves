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

FFT = True
FFT = False
COWLING = True
COWLING = False

filepath = "../"

fig = plt.figure(1)
fig, (ax1,ax2) = plt.subplots(2, 1)
plt.subplots_adjust(bottom=0.13, left=0.14, right=0.96, top=0.96, hspace=0.45)
plt.subplots_adjust(right=0.96, top=0.96, hspace=0.40)


df_TD = pd.read_csv(filepath+'output.log', index_col = False, delim_whitespace=True, usecols=['global_time','rho_max','alp_min','M_rest'])
df_TD = df_TD[['global_time','rho_max','alp_min','M_rest']]
data_time = df_TD['global_time'].to_numpy()
data_rho_max = df_TD['rho_max'].to_numpy()
y = data_rho_max / data_rho_max[0] - 1.0
y = y * 1000.
x = data_time / 2.03001708e05 * 1000 # s to ms

ax1.set_ylabel('$[\\rho_{c}(t)/\\rho_{c}(0) - 1] \\times 10^{3}$')
#ax1.set_ylabel('$[M(t)/M(0) - 1]$')
ax1.set_xlabel('$t$ $\\rm{(ms)}$')
ax1.set_xlim(x[0],x[-1])
#ax1.set_ylim(-3,0)

ax1.plot(x , y, linestyle = '-' ,linewidth = 1.0 )
#ax1.legend(loc='best')
ax1.grid(True)

if (FFT):
  #FFT part
  df_FD = pd.read_csv("vel_FD.csv", delimiter = '\t', usecols = ['f','fft1','psd1','fft2','psd2'])
  df_FD = df_FD[['f','fft1','psd1','fft2','psd2']]
  ax2.plot(0.001*df_FD['f'], df_FD['fft1'] ,linestyle = '-' ,linewidth = 1.0, label = '$\\rm{FFT}$$(v^{r})$' )
  ax2.plot(0.001*df_FD['f'], df_FD['fft2'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\theta})$' )

  if (COWLING):
    ax2.axvline(x=2.328, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=4.300, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=6.320, color='g',label = '$H_2$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=8.153, color='b',label = '$H_3$', linestyle='--',linewidth = 0.6)
  else:
    ax2.axvline(x=1.168, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=4.030, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=1.685, color='g',label = '$^2f$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=1.004, color='b',label = '$i_{2}$', linestyle='--',linewidth = 0.6)
  
  ax2.legend(loc='center right')


ax2.set_ylabel('$\\rm{FFT\ of\ }$ $v^{r}(t)$ $\\rm{and\ }$ $v^{\\theta}(t)$')
ax2.set_xlabel('$f$ $\\rm{[kHz]}$')
ax2.grid(True)
ax2.set_yscale('log')
ax2.set_xlim(0.0,8.5)
#ax2.set_ylim(1E-8,1E-3)

fig.savefig("./central_rho_combine.png", bbox_inches="tight")
