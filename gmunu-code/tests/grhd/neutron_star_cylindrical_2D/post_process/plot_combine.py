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

#it,global_time,dt,rho,veloc1,eps,press,alp,beta1,psi = np.loadtxt(filepath+"output.log", unpack=True ,dtype='float')
df_TD = pd.read_csv(filepath+'output.log', index_col = False, delim_whitespace=True, usecols = ['global_time','rho_max'])
df_TD = df_TD[['global_time','rho_max']]

t = df_TD['global_time'].to_numpy() / 2.03001708e05 *1000 # s to ms
delta_rho = df_TD['rho_max'].to_numpy() 
delta_rho = delta_rho / delta_rho[0] - 1.0

ax1.plot(t , delta_rho, linestyle = '-' ,linewidth = 1.0 )
#ax1.legend(loc='best')
ax1.grid(True)
ax1.set_ylabel('$\\rho_{\\max}(t)/\\rho_{\\max}(0) - 1$')
#ax1.set_ylabel('$[M(t)/M(0) - 1]$')
ax1.set_xlabel('$t$ $\\rm{(ms)}$')
ax1.set_xlim(t[0],t[-1])

if (FFT):
  #FFT part
  df_FD = pd.read_csv("W_vel_FD.csv", delimiter = '\t', usecols = ['f','fft1','psd1','fft2','psd2'])
  df_FD = df_FD[['f','fft1','psd1','fft2','psd2']]
  ax2.plot(df_FD['f'], df_FD['fft1'] ,linestyle = '-' ,linewidth = 1.0, label = '$\\rm{FFT}$$(v^{r})$' )
  ax2.plot(df_FD['f'], df_FD['fft2'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(v^{\\theta})$' )

  if (COWLING):
    ax2.axvline(x=2328, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=4300, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=6320, color='g',label = '$H_2$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=8153, color='b',label = '$H_3$', linestyle='--',linewidth = 0.6)
  else:
    ax2.axvline(x=1168, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=4030, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=1685, color='g',label = '$^2f$', linestyle='--',linewidth = 0.6)
    ax2.axvline(x=1004, color='b',label = '$i_{2}$', linestyle='--',linewidth = 0.6)
  
  ax2.legend(loc='center right')

ax2.grid(True)
ax2.set_ylabel('$\\rm{FFT\ of\ }$ $v^{r}(t)$ $\\rm{and\ }$ $v^{\\theta}(t)$')
ax2.set_xlabel('$f$ $\\rm{(Hz)}$')
ax2.set_xlim(100,6500)

fig.tight_layout()
#plt.subplots_adjust(hspace=0.05, wspace=0.1)
fig.subplots_adjust(hspace=0.35)
fig.savefig("./central_rho_combine.png", bbox_inches="tight")
