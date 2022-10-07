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

plot_color = ["black", "r","g","b","fuchsia"]
ls = ["-","-.", "--", ":","-."]
lw = [1,1.5,1.5,1.5,1.5]

#read data
df_TD = pd.read_csv('../output.log', index_col = False, delim_whitespace=True, usecols = ['global_time','EB_tor','EB_pol'])
df_TD = df_TD[['global_time','EB_tor','EB_pol']]

t=df_TD['global_time'].to_numpy() / 2.03001708e05 # s
EB_tor=df_TD['EB_tor'].to_numpy()
EB_pol=df_TD['EB_pol'].to_numpy()

# plot
plt.plot(t , EB_tor / EB_tor[0], linestyle = '-' ,linewidth = 1.0)
#plt.plot(t , EB_tor[-1] / EB_tor[0] * np.exp( - t / 2.03E13 )**2 , linestyle = '--' ,linewidth = 1.1 )
#plt.plot(t , np.exp( - t / 1.E8 )**2. , linestyle = '--' ,linewidth = 1.1 )

#plt.plot(t , EB_tor, linestyle = '-' ,linewidth = 1.0)
#plt.plot(t , EB_tor[-1] * np.exp( - t[-1] / 1.0E8 )**2 * np.exp( - t / 1.0E8 )**2 , linestyle = '--' ,linewidth = 1.1 )

plt.xlim(0,t[-1])
plt.xlabel('$t$ $\\rm{(s)}$')
plt.ylabel('$\\mathcal{E}_{\\rm{tor}}(t) / \\mathcal{E}_{\\rm{tor}}(0)$')

plt.yscale('log')
#plt.xlim(0,0.5)
#plt.ylim(1e-3,2e0)
plt.legend(loc='best')
plt.grid(True)


plt.savefig("./EB_tor.png", bbox_inches="tight")


############
### plot two
plt.clf()

plt.plot(t , EB_pol, linestyle = '-' ,linewidth = 1.0)

plt.xlim(0,t[-1])

plt.yscale('log')
plt.legend(loc='best')
plt.grid(True)


plt.savefig("./EB_pol.png", bbox_inches="tight")

