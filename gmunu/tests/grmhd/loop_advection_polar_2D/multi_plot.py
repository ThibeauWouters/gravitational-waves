import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot
from matplotlib import rc_context

import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

# Load the dataset
ds = yt.load('output0100.dat', unit_system='code')

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

fields=['B2','mag_press','W_vel1','W_vel2']
p = yt.SlicePlot(ds, 'z', fields)
p.set_log('Bvec1', False)
p.set_log('Bvec2', False)
p.set_log('W_vel1', False)
p.set_log('W_vel2', False)
p.set_log('divB', False)
p.set_log('B2', False)
#p.set_log('mag_press', False)
p.set_zlim("mag_press", 1E-12, 1E-7)

#p.set_cmap(field="rho", cmap='hot')
#p.annotate_timestamp(corner='upper_right',draw_inset_box=True)
#p.annotate_timestamp(corner='upper_left', draw_inset_box=True)

# For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
# axes.
for i, field in enumerate(fields):
	p.set_xlabel('$x$')
	p.set_ylabel('$y$')
	plot = p.plots[field]
	plot.figure = fig
	plot.axes = grid[i].axes
	plot.cax = grid.cbar_axes[i]

p._setup_plots()

with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('result.png')

