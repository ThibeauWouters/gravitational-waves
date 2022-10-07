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

savepath = './'

filepath = "../"

simname = [ \
        "output.log" ,\
        ]

plot_color = ["black", "r","g","b","fuchsia"]
ls = ["-","-.", "--", ":","-."]
lw = [1,1.5,1.5,1.5,1.5]

# read data
df_TD = pd.read_csv(filepath+'output.log', index_col = False, delim_whitespace=True, usecols = ['global_time','EB_tor','EB_tot'])
df_TD = df_TD[['global_time','EB_tor','EB_tot']]
data_time = df_TD['global_time'].to_numpy() / 2.03001708e05 * 1000 # s to ms
EB_tor = df_TD['EB_tor'].to_numpy()
EB_tot = df_TD['EB_tot'].to_numpy()
EB_pol = EB_tot - EB_tor

fig = plt.figure(1)
fig, (ax1) = plt.subplots(1, 1)
#plt.subplots_adjust(bottom=0.13, left=0.14, right=0.96, top=0.96, hspace=0.45)
#plt.subplots_adjust(right=0.96, top=0.96, hspace=0.40)

#ax1.set_ylabel('$[\\rho_{c}(t)/\\rho_{c}(0) - 1]$')
ax1.set_xlabel('$t$ $\\rm{(ms)}$')
#ax1.set_xlim(0,10)
#ax1.set_ylim(-3,0)

ax1.plot(data_time , EB_tor / EB_tot[0], linestyle = '-' ,linewidth = 1.0, label='$\\mathcal{E}_{\\rm{tor}}(t) / \\mathcal{E}_{\\rm{mag}}(0)$' )
ax1.plot(data_time , EB_pol / EB_tot[0], linestyle = '-' ,linewidth = 1.0, label='$\\mathcal{E}_{\\rm{tor}}(t) / \\mathcal{E}_{\\rm{mag}}(0)$' )
ax1.legend(loc='best')
ax1.grid(True)


fig.savefig(savepath + "B_energy.png", bbox_inches="tight")
