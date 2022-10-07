import matplotlib
import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import rc_context
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

run_name = [
	'k1_', 'k5_', 'k10_'
]

data_num = [
	'0000.dat',
	'0001.dat',
	'0002.dat',
	'0003.dat',
	'0004.dat',
	'0005.dat',
	'0006.dat',
	'0007.dat',
	'0008.dat'
]

mylabel = ['$\\kappa = 1$','$\\kappa = 5$','$\\kappa = 10$']
ls = ["-","-.",":", "--","-."]

y_max = 1.0
x_max = 2.0

fig, axs = plt.subplots(3,3, sharex='col',sharey='row')

for j, rn in enumerate(run_name):
   for i, dn in enumerate(data_num):
      fn = rn + dn
      print( fn )
      # Load the dataset
      ds = yt.load(fn, unit_system='code')
      nx = ds.domain_dimensions[0]
      x_projection= yt.LineBuffer( ds, (ds.domain_left_edge[0],0.0,0.0), (ds.domain_right_edge[0],0.0,0.0), nx )
      axs[int(i/3),int(i%3)].text(0.3*x_max, 0.7*y_max, '$t = %1.1f $' % ds.current_time)
      axs[int(i/3),int(i%3)].plot(x_projection[('amrvac', 'x')], x_projection[('amrvac', 'divB')], label=mylabel[j], ls=ls[j])#,lw=1.0)
      axs[int(i/3),int(i%3)].set_xticks([-2,0,2])
      axs[int(i/3),int(i%3)].set_xticklabels([])
      axs[int(i/3),int(i%3)].set_yticks([-1,0,1])
      axs[int(i/3),int(i%3)].set_yticklabels([])
      axs[int(i/3),int(i%3)].set_xlim(-x_max, x_max)
      axs[int(i/3),int(i%3)].set_ylim(-y_max, y_max)
      axs[int(i/3),int(i%3)].grid(True)


axs[2,1].set(xlabel='$x$')
axs[2,1].set_xticklabels([-2,0,2])
axs[1,0].set(ylabel='$\\nabla \\cdot \\vec{B}$')
axs[1,0].set_yticklabels([-1,0,1])

axs[0,1].legend(loc='upper center', bbox_to_anchor=(0.5,1.5), ncol=3)
plt.subplots_adjust(hspace=0.0, wspace=0.0)
plt.savefig('slices.png')

