import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot
from matplotlib import rc_context

import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import AxesGrid

def gamma(field, data):
    return  np.sqrt( 1.0 + (data['W_vel1']**2 + data['W_vel2']**2) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']
def veloc2(field, data):
    return data['W_vel2'] / data['gamma']
def B2(field, data):
    return  (data['Bvec1']**2 + data['Bvec2']**2)
def abs_divB_dx_over_B(field, data):
    return ( np.abs( data['divB'] / np.sqrt( data['B2'] ) ) * 0.05 )

# Load the dataset
ds = yt.load('output0040.dat', unit_system='code')
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')
ds.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')
ds.add_field( ('amrvac','abs_divB_dx_over_B'), function=abs_divB_dx_over_B, sampling_type='cell')

#plot = yt.LinePlot( ds, [('amrvac', 'gamma')], (-6.0,0.0), (6.0,0.0), 256 )
x_projection= yt.LineBuffer( ds, (-6.0,0.0), (6.0,0.0), 256 )
y_projection= yt.LineBuffer( ds, (0.0,-6.0), (0.0,6.0), 256 )

fig, axs = plt.subplots(4,2, sharex='col',sharey='row')
# x projection
axs[0,0].plot(x_projection[('amrvac', 'x')], x_projection[('amrvac', 'rho')],'*',lw=0.1)
#axs[0,0].set_ylim(0.0, 1.2)
axs[0,0].set( ylabel='$\\rho$')

axs[1,0].plot(x_projection[('amrvac', 'x')], x_projection[('amrvac', 'press')],'*',lw=0.1)
axs[1,0].set( ylabel='$P$')

axs[2,0].plot(x_projection[('amrvac', 'x')], x_projection[('amrvac', 'B2')],'*',lw=0.1)
axs[2,0].set( ylabel='$B^2$')

axs[3,0].plot(x_projection[('amrvac', 'x')], x_projection[('amrvac', 'gamma')],'*',lw=0.1)
axs[3,0].set(xlabel='$x$', ylabel='$\\Gamma$')

# y projection
axs[0,1].plot(y_projection[('amrvac', 'y')], y_projection[('amrvac', 'rho')],'*',lw=0.1)

axs[1,1].plot(y_projection[('amrvac', 'y')], y_projection[('amrvac', 'press')],'*',lw=0.1)

axs[2,1].plot(y_projection[('amrvac', 'y')], y_projection[('amrvac', 'B2')],'*',lw=0.1)

axs[3,1].plot(y_projection[('amrvac', 'y')], y_projection[('amrvac', 'gamma')],'*',lw=0.1)
axs[3,1].set(xlabel='$y$')

fig.tight_layout()

plt.subplots_adjust(hspace=0.1, wspace=0.0)

with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('slices.png')

#plt.subplot(4, 1, 1)
#plt.plot(plot[('amrvac', 'x')],plot[('amrvac', 'gamma')])

