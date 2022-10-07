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

#sigma~1E2 -> tau ~ 2.3E2

plot_color = ["black", "r","g","b","fuchsia"]
ls = ["-","-.", "--", ":","-."]
lw = [1,1.5,1.5,1.5,1.5]

#plt.ylabel('$[M(t)/M(0) - 1]$')
#plt.xlabel('$t$ $\\rm{(ms)}$')

plt.xlabel('$t$')
#plt.xlim(0,0.002)
#plt.ylim(1e-6,1e-5)

df_TD = pd.read_csv('../output.log', index_col = False, delim_whitespace=True, usecols = ['global_time','B_avg','EB_tor','EB_pol','divB_avg', 'divE_avg'])
df_TD = df_TD[['global_time','B_avg','EB_tor','EB_pol','divB_avg', 'divE_avg']]

t=df_TD['global_time'].to_numpy() / 2.03001708e05 # s
B_avg=df_TD['B_avg'].to_numpy()
EB_tor=df_TD['EB_tor'].to_numpy()
EB_pol=df_TD['EB_pol'].to_numpy()
divB=df_TD['divB_avg'].to_numpy()
divE=df_TD['divE_avg'].to_numpy()

plt.plot(t , EB_tor / EB_tor[0], linestyle = '-' ,linewidth = 1.0)

# plot expected decay rate
sigma_cgs =  2.03E4 * 10.**1
tau_diss = sigma_cgs * 1E-8 * 0.13
ind_of_exp = 1.0
#prefactor = (EB_tor[-1] / EB_tor[0] / np.exp( - t[-1] / tau_diss )**ind_of_exp)
prefactor = 1.0
plt.plot(t , prefactor * np.exp( - t / tau_diss )**ind_of_exp, linestyle = '--' ,linewidth = 1.0 )

plt.xlim(0,t[-1])
plt.xlabel('$t$ $\\rm{(s)}$')
plt.ylabel('$\\mathcal{E}_{\\rm{tor}}(t) / \\mathcal{E}_{\\rm{tor}}(0)$')

plt.yscale('log')
#plt.xlim(0,0.5)
plt.ylim(1e-1,2e0)
#plt.legend(loc='best')
plt.grid(True)


plt.savefig("./EB_tor.png", bbox_inches="tight")


############
### plot two
plt.clf()

plt.plot(t , divE, linestyle = '-' ,linewidth = 1.0)
#plt.xlim(0,t[-1])
plt.yscale('log')
plt.ylim(1e-13,1e-10)
#plt.legend(loc='best')
plt.grid(True)

plt.savefig("./divE.png", bbox_inches="tight")
