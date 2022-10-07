import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':15})
from matplotlib import rc_context
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

import yt
import numpy as np

from scipy.special import erf

eta = 1.0E-2
T_final = 1.0E1

def B_sol(x,t):
    return erf( 0.5 * x / np.sqrt(eta * t) )
def E_sol(x,t):
    return np.sqrt(eta/np.pi/t) * np.exp( - 0.25 * x**2 / (eta * t) )

def gamma(field, data):
    return  np.sqrt( 1.0 + (data['W_vel1']**2 + data['W_vel2']**2) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']
def veloc2(field, data):
    return data['W_vel2'] / data['gamma']


print ( erf(1.) )
print ( B_sol(0.,1.) )


# Load the dataset
ds = yt.load('output0001.dat')
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')

# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['x'])

#print (ds.field_list)
#print (ds.derived_field_list)

# plot the data
fig, axs = plt.subplots(3,1, sharex='col')

axs[0].plot(np.array(plotdata['x'][srt]), np.array(B_sol(np.array(plotdata['x'][srt]),T_final)), color="black",linestyle = '-',linewidth = 3 ) 
axs[0].plot(np.array(plotdata['x'][srt]), np.array(plotdata['Bvec2'][srt]), color="red", linestyle = ':', linewidth=2.5, label='By')
axs[0].plot(np.array(plotdata['x'][srt]), np.array(plotdata['Bvec1'][srt]), linestyle = ':', linewidth=2.5, label='Bx')
axs[0].plot(np.array(plotdata['x'][srt]), np.array(plotdata['Bvec3'][srt]), linestyle = ':', linewidth=2.5, label='Bz')
axs[0].set(xlabel='$x$', ylabel='$B^y$')
axs[0].set_xlim(-1.5,1.5)
axs[0].grid(True)

axs[1].plot(np.array(plotdata['x'][srt]), np.array(E_sol(np.array(plotdata['x'][srt]),T_final)), color="black",linestyle = '-',linewidth = 3 ) 
axs[1].plot(np.array(plotdata['x'][srt]), np.array(plotdata['Evec1'][srt]), linestyle = ':', linewidth=2.5, label='Ex')
axs[1].plot(np.array(plotdata['x'][srt]), np.array(plotdata['Evec2'][srt]), linestyle = ':', linewidth=2.5, label='Ey')
axs[1].plot(np.array(plotdata['x'][srt]), np.array(plotdata['Evec3'][srt]), color="red", linestyle = ':', linewidth=2.5, label='Ez')
#axs[1,0].set_ylim(0.0, 1.2)
axs[1].set(xlabel='$x$', ylabel='$E^z$')
axs[1].grid(True)
axs[1].legend(loc='upper right')

axs[2].plot(np.array(plotdata['x'][srt]), np.array(plotdata['W_vel1'][srt]), linestyle = ':', linewidth=2.5, label='Wvx')
axs[2].plot(np.array(plotdata['x'][srt]), np.array(plotdata['W_vel2'][srt]), linestyle = ':', linewidth=2.5, label='Wvy')
axs[2].plot(np.array(plotdata['x'][srt]), np.array(plotdata['W_vel3'][srt]), color="red", linestyle = ':', linewidth=2.5, label='Wvz')
#axs[1,0].set_ylim(0.0, 1.2)
axs[2].grid(True)
axs[2].legend(loc='upper right')

fig.tight_layout()

#plt.legend(loc='best')
with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('results.png')
