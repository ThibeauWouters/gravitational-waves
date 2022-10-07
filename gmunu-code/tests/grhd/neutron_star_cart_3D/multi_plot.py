import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

savepath = '/users/chi-kit.cheong/public_html/'

units_override = dict( 
length_unit=(1.0, 'cm'),
time_unit=(1.0, 's'),
mass_unit=(1.0, 'g'),
)

fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 2),
                axes_pad = 1.0,
                label_mode = "each",
                share_all = False,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="5%",
                cbar_pad="0%")

# Load the dataset
ds = yt.load('output0086.dat', unit_system='code')
cuts=['z','z','x','x']
#fields=['press','rho',('index', 'grid_level'),'rho']
fields=['press','rho','W_vel2','rho']
annotate_flag=[True,False,True,False]

for i, (direction, field, annotate) in enumerate(zip(cuts, fields, annotate_flag)):
    # Load the data and create a single plot
    p = yt.SlicePlot(ds, direction, field)
    p.zoom(3)
    p.set_cmap(field="rho", cmap='hot')
    p.set_cmap(field="press", cmap='afmhot')
    p.set_log('alp', False)
    #p.set_log('psi', False)
    #p.set_log('beta1', False)
    p.set_log('W_vel1', False)
    p.set_log('W_vel2', False)
    #p.set_log('beta3', False)
    p.set_log(('index', 'grid_level'), False)
    if (annotate and direction == 'z'):
        #p.annotate_streamlines('W_vel1', 'W_vel2', density = 0.5, factor = 5)
        p.annotate_streamlines('W_vel1', 'W_vel2')
        p.annotate_grids()
    elif (annotate and direction == 'x'):
        p.annotate_grids()

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Since there are only two colorbar axes, we need to make sure we don't try
    # to set the temperature colorbar to cbar_axes[4], which would if we used i
    # to index cbar_axes, yielding a plot without a temperature colorbar.
    # This unnecessarily redraws the Density colorbar three times, but that has
    # no effect on the final plot.
    #if field == "rho":
    #    plot.cax = grid.cbar_axes[0]
    #elif field == "temperature":
    #    plot.cax = grid.cbar_axes[1]

    # Finally, redraw the plot.
    p._setup_plots()

plt.savefig(savepath + "NS_multiplot.pdf")
