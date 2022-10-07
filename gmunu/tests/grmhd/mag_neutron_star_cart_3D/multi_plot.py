import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def gamma(field, data):
    return  np.sqrt( 1.0 + data['psi']**4 * (data['W_vel1']**2 + data['W_vel2']**2) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']
def veloc2(field, data):
    return data['W_vel2'] / data['gamma']
def B2_tor(field, data):
    return  data['psi']**4 * (data['Bvec1']**2 + data['Bvec2']**2)
def B2(field, data):
    return  data['psi']**4 * (data['Bvec1']**2 + data['Bvec2']**2 + data['Bvec3']**2)

savepath = './'

units_override = dict( 
length_unit=(1.0, 'cm'),
time_unit=(1.0, 's'),
mass_unit=(1.0, 'g'),
)

fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 2),
                axes_pad = 2.0,
                label_mode = "each",
                share_all = False,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="5%",
                cbar_pad="0%")

# Load the dataset
ds = yt.load('output0005.dat', unit_system='code')
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
#ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
#ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')
ds.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')
#ds.add_field( ('amrvac','mag_press'), function=mag_press, sampling_type='cell')

cuts=['z','z','x','x']
fields=['divB','rho','B2','rho']
annotate_flag=[False,True,True,False]

for i, (direction, field, annotate) in enumerate(zip(cuts, fields, annotate_flag)):
    # Load the data and create a single plot
    p = yt.SlicePlot(ds, direction, field)
    p.zoom(3)
    p.set_cmap(field="rho", cmap='hot')
    #p.set_log('Bvec1', False)
    if (annotate and direction == 'z'):
        p.annotate_streamlines('Bvec1', 'Bvec2')
    #if (annotate and direction == 'x'):
    #    p.annotate_grids()

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

plt.savefig(savepath + "hydro.pdf")
