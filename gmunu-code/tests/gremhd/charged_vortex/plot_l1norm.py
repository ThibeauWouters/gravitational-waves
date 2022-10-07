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

df_TD = pd.read_csv('./output.log', index_col = False, delim_whitespace=True, usecols = ['global_time','q_tot','L1_norm_q'])
df_TD = df_TD[['global_time','q_tot','L1_norm_q']]

t=df_TD['global_time'].to_numpy()
L1_norm_q=df_TD['L1_norm_q'].to_numpy()
q_tot=df_TD['q_tot'].to_numpy()


fig = plt.figure(1)
fig, (ax1,ax2) = plt.subplots(2, 1)

ax1.plot(t , np.abs(1.0 - q_tot / q_tot[0]) )
ax1.set_ylabel('$ | 1 - q / q (0) | $')
ax1.set_yscale('log')
ax1.grid(True)

ax2.plot(t , L1_norm_q)
ax2.set_xlabel('$t$')
ax2.set_ylabel('$ || \\rho_e - {\\rho_e}(0) ||_1$')
ax2.set_yscale('log')
ax2.grid(True)

fig.savefig("./charge_conservation.png", bbox_inches="tight")


