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
COWLING = False
COWLING = True

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

#ax1.set_ylabel('$[\\rho_{c}(t)/\\rho_{c}(0) - 1] \\times 10^{3}$')
ax1.set_ylabel('$[\\rho_{c}(t)/\\rho_{c}(0) - 1]$')
ax1.set_xlabel('$t$ $\\rm{(ms)}$')
ax1.set_xlim(0,10)
#ax1.set_ylim(-3,0)

ax2.set_ylabel('$\\rm{FFT\ of\ }$ $v^{r}(t)$ $\\rm{and\ }$ $v^{\\theta}(t)$')
ax2.set_xlabel('$f$ $\\rm{(Hz)}$')
ax2.set_xlim(100,8500)
ax2.set_ylim(0,0.04)

#it,global_time,dt,rho,veloc1,eps,press,alp,beta1,psi = np.loadtxt(filepath+"output.log", unpack=True ,dtype='float')
x, y = np.loadtxt(filepath+"output.log", unpack=True ,dtype='float',usecols=(1,3),skiprows=1)
y = y / y[0] - 1.0
x = x / 2.03001708e05 * 1000 # s to ms
#y = y * 1E3 # 

ax1.plot(x , y, linestyle = '-' ,linewidth = 1.0 )
#ax1.legend(loc='best')
ax1.grid(True)

if (FFT):
  #FFT part
  df_FD = pd.read_csv("W_vel_FD.csv", delimiter = '\t', usecols = ['f','fft1','psd1','fft2','psd2'])
  df_FD = df_FD[['f','fft1','psd1','fft2','psd2']]
  ax2.plot(df_FD['f'], df_FD['fft1'] ,linestyle = '-' ,linewidth = 1.0, label = '$\\rm{FFT}$$(Wv^{r})$' )
  ax2.plot(df_FD['f'], df_FD['fft2'] ,linestyle = '-' ,linewidth = 1.0, label='$\\rm{FFT}$$(Wv^{\\theta})$' )

  #if (COWLING):
    #ax2.axvline(x=2328, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=4300, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=6320, color='g',label = '$H_2$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=8153, color='b',label = '$H_3$', linestyle='--',linewidth = 0.6)
  #else:
    #ax2.axvline(x=1168, color='black',label = '$F$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=4030, color='r',label = '$H_1$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=1685, color='g',label = '$^2f$', linestyle='--',linewidth = 0.6)
    #ax2.axvline(x=1004, color='b',label = '$i_{2}$', linestyle='--',linewidth = 0.6)
  
  ax2.legend(loc='best')
  ax2.grid(True)


fig.savefig("./central_rho_combine.png", bbox_inches="tight")
