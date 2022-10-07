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

vector_factor = 0.1

def gamma(field, data):
    return  np.sqrt( 1.0 + data['psi']**4 * (data['W_vel1']**2 + data['W_vel2']**2 + data['r']**2*data['W_vel3']**2) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']
def veloc2(field, data):
    return data['W_vel2'] / data['gamma']
def B2_tor(field, data):
    return  data['psi']**4 * (data['r']**2 * data['Bvec3']**2)
def B2_pol(field, data):
    return  data['psi']**4 * ( data['Bvec1']**2 + data['Bvec2']**2 )
def B_tor(field, data):
    return  data['psi']**2 * (data['r'] * data['Bvec3'])
def B_pol(field, data):
    return  data['psi']**2 * np.sqrt( data['Bvec1']**2 + data['Bvec2']**2 )
def B1(field, data):
    return  np.log( np.abs(data['Bvec1']) ) * data['Bvec1']
def B2(field, data):
    return  np.log( np.abs(data['Bvec2']) ) * data['Bvec2']

units_override = dict( 
length_unit=(1.0/6.77140812e-06, 'cm'),
time_unit=(1.0/2.03001708e05, 's'),
mass_unit=(1.0/5.02765209e-34, 'g'),
)

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
#ds = yt.load('output0000.dat', geometry_override='cylindrical', unit_system='code')
ds = yt.load('output0000.dat', geometry_override='cylindrical', units_override=units_override)
#ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','B_tor'), function=B_tor, sampling_type='cell', units='code_length')
ds.add_field( ('amrvac','B2_tor'), function=B2_tor, sampling_type='cell', units='code_length**2')
ds.add_field( ('amrvac','B2_pol'), function=B2_pol, sampling_type='cell')
ds.add_field( ('amrvac','B_pol'), function=B_pol, sampling_type='cell')
ds.add_field( ('amrvac','B1'), function=B1, sampling_type='cell')
ds.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')

fields=['rho','divB','B2_pol','B_pol']
plot_width = 30
plot_height = 30
p = yt.SlicePlot(ds,2, fields)
p.set_log('W_vel1', False)
p.set_log('W_vel2', False)
p.set_log('W_vel3', False)
p.set_log('Bvec1', False)
p.set_log('Bvec2', False)
p.set_log('Bvec3', False)
p.set_log('B_tor', False)
p.set_log('B2_tor', False)
#p.set_log('B2_pol', False)
p.set_cmap(field="rho", cmap='hot')
p.set_cmap(field="B_tor", cmap='jet')
p.set_cmap(field="B_pol", cmap='jet')
p.set_cmap(field="B2_tor", cmap='jet')
p.set_cmap(field="B2_pol", cmap='jet')
#p.set_zlim("B2_tor", 0, 1E-13)
#p.set_zlim("B_pol", 2E-6, 1E-9)
p.set_zlim("B2_pol", 2E-7, 1E-13)
#p.set_zlim("B_pol", 1E-6, 1E-9)
#p.set_zlim("B2_pol", 1E-12, 1E-18)
# full domain
p.set_center((plot_width/2.0, 0.0))
# half domain
#p.set_center((plot_width/2.0, plot_height/2.0))
p.set_width(plot_width, plot_height)
p.annotate_timestamp(corner='upper_right',draw_inset_box=True)
p.annotate_contour('rho',take_log=True,
    ncont=1,
    clim=(1e-5, 1e-4),
    label=False,
    plot_args={"colors": "red", "linewidths": 2},
)
#p.annotate_streamlines('Bvec1', 'Bvec2', factor=16, density=1, display_threshold=1e-9, field_color = 'rho', plot_args={ "linewidth": 0.9})
p.annotate_streamlines('Bvec1', 'Bvec2', factor=16, density=1, plot_args={ "linewidth": 0.9})
#p.annotate_quiver('Bvec1', 'Bvec2', factor=20, normalize=False, 
#    plot_args={"color": "black" 
    #, "linewidths": 0.2
#    }
#    )
#p.annotate_streamlines('Evec1', 'Evec2')
#p.annotate_contour('Bvec3',take_log=False, clim=(0,2e-3), ncont=10, label=True)

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

