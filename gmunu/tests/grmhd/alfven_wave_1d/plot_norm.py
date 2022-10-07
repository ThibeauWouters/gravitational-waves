import yt
import numpy as np
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot
from matplotlib import rc_context

def one_over_x2_pt(x):
    return (1.0E0 / x**2 )
def one_over_x3_pt(x):
    return (1.0E0 / x**3 )
one_over_x2 = np.vectorize(one_over_x2_pt)
one_over_x3 = np.vectorize(one_over_x3_pt)

df = pd.read_csv("norms_wenoyc3.csv", delimiter = '\t', usecols = ['nx','L1_norm'])
#df['nx'] = df['nx'].astype('int')
df = df[['nx','L1_norm']]

log_nx = np.log10(df['nx'])
log_err = np.log10(df['L1_norm'])

inx_max = df['nx'].idxmax()
ref_x = [log_nx[0], log_nx[inx_max]]

order = 2.0
shift_value = log_err[0] + order * log_nx[0]
ref_2nd_y = [ log_err[0], - order * log_nx[inx_max] + shift_value]
order = 3.0
shift_value = log_err[0] + order * log_nx[0]
ref_3rd_y = [ log_err[0], - order * log_nx[inx_max] + shift_value]

conv_rate = np.array([2.0])
for i in range(1,inx_max):
    conv_rate = np.append(conv_rate,[ np.abs(np.log10(0.5*df['L1_norm'][i]/df['L1_norm'][i-1])) ] )


# now, plot the date
plt.plot(log_nx, log_err, color='red', linestyle='dashed', marker='o',markerfacecolor='black', markersize=8)
plt.plot(ref_x,ref_2nd_y,color='grey', linestyle='-.',label='$ \\propto N^{-2}$')
plt.plot(ref_x,ref_3rd_y,color='blue', linestyle='-.',label='$ \\propto N^{-3}$')
plt.xlabel('$\\log_{10} (N)$')
plt.ylabel('$\\log_{10} (L1 norm)$')
#plt.xlim(0,1)
plt.legend(loc='best')
plt.grid(True)
#plt.show()
plt.savefig('l1_norm.png')
plt.clf()

plt.plot(df['nx'][1:inx_max],conv_rate[1:inx_max], color='black', linestyle='dashed', marker='o',markerfacecolor='black', markersize=8)
plt.axhline(y=2.0, color='grey', linestyle='-')
plt.xlabel('$N$')
plt.ylabel('L1 rate')
plt.xscale('log')
plt.grid(True)
plt.savefig('l1_order.png')

