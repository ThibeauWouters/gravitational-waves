import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.colors import LogNorm

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



level=1

# Load the dataset
ds_i = yt.load('None0000.dat', unit_system='code')
ds_i.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds_i.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds_i.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')
ds_i.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')
ds_i.add_field( ('amrvac','mag_press'), function=mag_press, sampling_type='cell')
mag_prs_i = ds_i.covering_grid(level, left_edge=ds_i.domain_left_edge, dims=ds_i.domain_dimensions * ds_i.refine_by**level)['mag_press']
#mag_prs_i = ds_i.covering_grid(level, left_edge=ds_i.domain_left_edge, dims=ds_i.domain_dimensions * ds_i.refine_by**level)['B2']

ds_f = yt.load('None0001.dat', unit_system='code')
ds_f.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds_f.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds_f.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')
ds_f.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')
ds_f.add_field( ('amrvac','mag_press'), function=mag_press, sampling_type='cell')
mag_prs_f = ds_f.covering_grid(level, left_edge=ds_f.domain_left_edge, dims=ds_f.domain_dimensions * ds_f.refine_by**level)['mag_press']
#mag_prs_f = ds_f.covering_grid(level, left_edge=ds_f.domain_left_edge, dims=ds_f.domain_dimensions * ds_f.refine_by**level)['B2']

fig, ax = plt.subplots(1)
plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
out_plot = np.abs( mag_prs_f[:,:,0].value - mag_prs_i[:,:,0].value ) #/ mag_prs_f[:,:,0].value
#image = ax.imshow(np.rot90(np.log10(out_plot)))
image = ax.imshow(np.rot90(out_plot), norm=LogNorm(), vmin=1E-20, vmax = 1E-6)
ax.set_axis_off()

fig.colorbar(image, cax=cax, orientation='vertical')

plt.savefig('compare.png')

