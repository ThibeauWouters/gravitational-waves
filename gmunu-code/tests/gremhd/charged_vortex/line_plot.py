import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot
from matplotlib import rc_context

import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import AxesGrid

from matplotlib import rc_context
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def charge(x):
    return 0.7 / (x**2 + 1.)**2

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
ds = yt.load('output0001.dat', unit_system='code')
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')
ds.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')
ds.add_field( ('amrvac','abs_divB_dx_over_B'), function=abs_divB_dx_over_B, sampling_type='cell')

nx = ds.domain_dimensions[0]
ny = ds.domain_dimensions[1]

x_min = ds.domain_left_edge[0]
x_max = ds.domain_right_edge[0]
y_min = ds.domain_left_edge[1]
y_max = ds.domain_right_edge[1]

x_projection= yt.LineBuffer( ds, (x_min,0.0), (x_max,0.0), nx )
y_projection= yt.LineBuffer( ds, (0.0,y_min), (0.0,y_max), ny )

fig, axs = plt.subplots(2,1,sharey='row')
# x projection
axs[0].plot(np.array(x_projection[('amrvac','x')]), np.array(charge(np.array(x_projection[('amrvac','x')]))), color="black",linestyle = '-',linewidth = 2 ) 
axs[0].plot(x_projection[('amrvac', 'x')], x_projection[('amrvac', 'divE')],'o',color='red',markersize=3.5)
axs[0].set(xlabel='$x$', ylabel='$q$')
axs[0].set(xlim=(x_min,x_max))
axs[0].grid(True)

# y projection
axs[1].plot(np.array(y_projection[('amrvac','y')]), np.array(charge(np.array(y_projection[('amrvac','y')]))), color="black",linestyle = '-',linewidth = 2 ) 
axs[1].plot(y_projection[('amrvac', 'y')], y_projection[('amrvac', 'divE')],'o',color='red',markersize=3.5)
axs[1].set(xlabel='$y$', ylabel='$q$')
axs[1].set(xlim=(y_min,y_max))
axs[1].grid(True)

fig.tight_layout()

#plt.subplots_adjust(hspace=0.2, wspace=0.0)

with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('slices.png')

#plt.subplot(4, 1, 1)
#plt.plot(plot[('amrvac', 'x')],plot[('amrvac', 'gamma')])

