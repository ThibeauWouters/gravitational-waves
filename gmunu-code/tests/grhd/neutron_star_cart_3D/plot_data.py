import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np

#savepath = '/users/chi-kit.cheong/public_html/'
savepath = './'

units_override = dict( 
length_unit=(1.0, 'cm'),
time_unit=(1.0, 's'),
mass_unit=(1.0, 'g'),
)

# Load the dataset
ds = yt.load('output0004.dat', unit_system='code')

ds.print_stats()

proj = ds.slice('x',0)
p = proj.to_pw('rho', origin='native')
p.save(savepath + 'yz.png')

proj = ds.slice('y',0)
p = proj.to_pw('rho', origin='native')
p.save(savepath + 'xz.png')

proj = ds.slice('z',0)
p = proj.to_pw('rho', origin='native')
p.annotate_streamlines('W_vel1', 'W_vel2', density = 0.5, factor = 10)
p.set_width(60,60)
p.save(savepath + 'xy.png')

proj = ds.slice('z',0)
#p.set_log(('index', 'grid_indices'), False)
p = proj.to_pw(('index', 'grid_level'), origin='native')
p.set_log(('index', 'grid_level'), False)

p.save(savepath + 'gridlevel_xy.png')
