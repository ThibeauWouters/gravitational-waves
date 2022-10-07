import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':35})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

savepath = './'

def gamma(field, data):
    return  np.sqrt( 1.0 + data['psi']**4 * (data['W_vel1']**2 + data['W_vel2']**2 + data['W_vel3']**2) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']
def veloc2(field, data):
    return data['W_vel2'] / data['gamma']

def angle_phi(field, data):
    return np.arctan2(data['y'],data['x'])

def B2_tor(field, data):
    return  data['psi']**4 * (- np.sin(data['angle_phi']) * data['Bvec1']
                              + np.cos(data['angle_phi']) * data['Bvec2'] )**2
def B2(field, data):
    return  data['psi']**4 * (data['Bvec1']**2 + data['Bvec2']**2 + data['Bvec3']**2)
def B2_pol(field, data):
    return data['B2'] - data['B2_tor']

units_override = dict( 
length_unit=(1.0/6.77140812e-06, 'cm'),
time_unit=(1.0/2.03001708e05, 's'),
mass_unit=(1.0/5.02765209e-34, 'g'),
)


# Load the dataset
#ds = yt.load('output0000.dat', unit_system='code')
ds = yt.load('output0000.dat', units_override=units_override)
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')
ds.add_field( ('amrvac','angle_phi'), function=angle_phi, sampling_type='cell')
ds.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')
ds.add_field( ('amrvac','B2_tor'), function=B2_tor, sampling_type='cell')
ds.add_field( ('amrvac','B2_pol'), function=B2_pol, sampling_type='cell')

ds.print_stats()

proj = ds.slice('x',0)
p = proj.to_pw('B2_pol', origin='native')
p.set_cmap(field="B2_pol", cmap='jet')
p.set_colorbar_label('B2_pol', '$ | B_{\\rm{pol}} |^2 $')
#p.set_log('B2_pol', False)
p.set_zlim('B2_pol', 1e-21, 1e-13)
p.annotate_timestamp(corner='upper_right',draw_inset_box=True)
p.annotate_contour('rho',take_log=True,
    ncont=1,
    clim=(1e-5, 1e-4),
    label=False,
    plot_args={"colors": "red", "linewidths": 2},
)
#p.annotate_streamlines('Bvec2', 'Bvec3', factor=16, density=1, plot_args={ "linewidth": 0.9})
p.annotate_streamlines('Bvec2', 'Bvec3', factor=16, density=1, plot_args={ "linewidth": 1.0, "color":'green'})
#p.set_font({'family': 'sans-serif', 'style': 'italic','weight': 'bold', 'size': 15})
#p.set_xlabel('$y$')
#p.set_ylabel('$z$')


p.save(savepath + 'yz.png')

#proj = ds.slice('y',0)
#p = proj.to_pw('rho', origin='native')
#p.save(savepath + 'xz.png')
#
#proj = ds.slice('z',0)
#p = proj.to_pw('rho', origin='native')
#p.annotate_streamlines('W_vel1', 'W_vel2', density = 0.5, factor = 10)
##p.set_width(60,60)
#p.save(savepath + 'xy.png')
#
#proj = ds.slice('z',0)
##p.set_log(('index', 'grid_indices'), False)
#p = proj.to_pw(('index', 'grid_level'), origin='native')
#p.set_log(('index', 'grid_level'), False)
#
#p.save(savepath + 'gridlevel_xy.png')
