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
def mag_press(field, data):
    return  0.5 * ( data['B2']/data['gamma']**2 + ( data['Bvec1']*data['veloc1'] + data['Bvec2']*data['veloc2'] )**2 )
def abs_divB_dx_over_B(field, data):
    return ( np.abs( data['divB'] / np.sqrt( data['B2'] ) ) * 0.05 )

# Load the dataset
ds = yt.load('output0010.dat', unit_system='code')
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')
ds.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')
ds.add_field( ('amrvac','mag_press'), function=mag_press, sampling_type='cell')
ds.add_field( ('amrvac','abs_divB_dx_over_B'), function=abs_divB_dx_over_B, sampling_type='cell')

fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 3),
                axes_pad = 1.0,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")

fields=['Bvec1','B2', 'divB']
p = yt.SlicePlot(ds,2, fields)
p.set_log('B2', False)
p.set_log('Bvec1', False)
p.set_log('mag_press', False)

#p.set_cmap(field="rho", cmap='hot')
#p.annotate_timestamp(corner='upper_right',draw_inset_box=True)
#p.annotate_timestamp(corner='upper_left', draw_inset_box=True)

# For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
# axes.
for i, field in enumerate(fields):
    p.set_xlabel('$x$')
    p.set_ylabel('$y$')
    #p.annotate_streamlines('Bvec1', 'Bvec2', density = 0.5, factor = 10)
    #p.annotate_contour('B2', label=True)
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

p._setup_plots()

with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('result.png')

