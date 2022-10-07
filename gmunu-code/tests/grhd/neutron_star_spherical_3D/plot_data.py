import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np


debug = True
debug = False

units_override = dict( 
length_unit=(1.0, 'cm'),
time_unit=(1.0, 's'),
mass_unit=(1.0, 'g'),
)

# Load the dataset
ds = yt.load('output0001.dat', geometry_override='spherical', unit_system='code')

Radius = ds.domain_width[0]

ds.print_stats()

r_slice = 1
theta_slice = 1.E-1
phi_slice = 1.E-1 #np.pi/4


proj = ds.slice('r', r_slice)
p = proj.to_pw('rho', origin='native')
p.save('r.png')

proj = ds.slice('phi', phi_slice)
p = proj.to_pw('W_vel3', origin='native')
p.set_center((Radius/2.0,Radius/2.0))
p.set_width(Radius,Radius)
p.save('phi.png')

proj = ds.slice('theta', theta_slice)
p = proj.to_pw('W_vel3', origin='native')
p.set_center((Radius/2,Radius/2))
p.set_width(Radius,Radius)

p.save('xy.png')
