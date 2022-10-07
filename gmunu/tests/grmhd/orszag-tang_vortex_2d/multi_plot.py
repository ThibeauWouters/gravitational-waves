import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot
from matplotlib import rc_context

import yt
import numpy as np
import matplotlib.pyplot as plt
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
ds = yt.load('output0485.dat', unit_system='code')
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')
ds.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')
ds.add_field( ('amrvac','abs_divB_dx_over_B'), function=abs_divB_dx_over_B, sampling_type='cell')

fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 2),
                axes_pad = 1.0,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")

fields=['gamma','B2','abs_divB_dx_over_B','divB']
p = yt.SlicePlot(ds,2, fields)
p.set_log('gamma', False)
p.set_log('B2', False)
#p.set_log('divB', False)
p.set_log('press', False)
#p.set_cmap(field="rho", cmap='hot')
#p.annotate_timestamp(corner='upper_right',draw_inset_box=True)
#p.annotate_timestamp(corner='upper_left', draw_inset_box=True)

# For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
# axes.
for i, field in enumerate(fields):
    p.set_xlabel('x-pi')
    p.set_ylabel('y-pi')
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

p._setup_plots()

with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('result.png')

