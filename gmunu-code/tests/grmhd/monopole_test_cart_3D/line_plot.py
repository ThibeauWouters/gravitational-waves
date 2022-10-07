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

fig, axs = plt.subplots(1,1, sharex='col',sharey='row')

# Load the dataset
ds = yt.load('output0007.dat', unit_system='code')
x_projection= yt.LineBuffer( ds, (-2.0,0.0,0.0), (2.0,0.0,0.0), 256 )

# x projection
axs.plot(x_projection[('amrvac', 'x')], x_projection[('amrvac', 'divB')])#,lw=1.0)
#axs.plot(x_projection[('amrvac', 'x')], x_projection[('amrvac', 'Bphi')])#,lw=1.0)
axs.set_xlim(-2.0, 2.0)
axs.set_ylim(-2.0, 2.0)
axs.set(xlabel='$x$', ylabel='$\\nabla \\cdot \\vec{B}$')

fig.tight_layout()

plt.grid(True)
plt.subplots_adjust(hspace=0.1, wspace=0.0)

with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('slices.png')

#plt.subplot(4, 1, 1)
#plt.plot(plot[('amrvac', 'x')],plot[('amrvac', 'gamma')])

