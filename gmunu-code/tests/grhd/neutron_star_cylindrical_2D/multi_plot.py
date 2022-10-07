import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot
from matplotlib import rc_context

import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

debug = True
debug = False


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

# Load the dataset
ds = yt.load('output0000.dat', geometry_override='cylindrical', unit_system='code')
#print (ds.field_list)
#print (ds.derived_field_list)

fields=['rho','eps', 'alp', 'grid_level']
Radius = ds.domain_width[0]
p = yt.SlicePlot(ds,2, fields)
p.set_log('W_vel1', False)
p.set_log('W_vel2', False)
p.set_log('W_vel3', False)
p.set_log('alp', False)
p.set_log('grid_level', False)
p.set_cmap(field="rho", cmap='hot')
p.set_cmap(field="eps", cmap='RED TEMPERATURE')
#p.set_cmap(field="Bvec3", cmap='Plasma')
#p.set_zlim("rho", 1e-13, 1.3E-3)
#Radius = 50
#p.set_center((Radius/2.0,0.))
#p.set_width(Radius,Radius)
p.annotate_timestamp(corner='upper_right',draw_inset_box=True)

# For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
# axes.
for i, field in enumerate(fields):
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

p._setup_plots()

with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('hydro.png')

