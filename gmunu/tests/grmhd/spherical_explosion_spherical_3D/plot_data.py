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
ds = yt.load('output0000.dat', geometry_override='spherical', unit_system='code')

Radius = ds.domain_width[0]

ds.print_stats()
#print(ds.field_list)
#print(ds.derived_field_list)

#im, sc = yt.volume_render(ds, "density")
#im.write_png("original.png", sigma_clip=8.0)

#yt.ProjectionPlot(ds, ('phi', np.pi/2.), "rho").save('proj.png')

proj = ds.slice('r', 0)
p = proj.to_pw('rho', origin='native')
p.save('r.png')

proj = ds.slice('phi', 0.0 * np.pi)
p = proj.to_pw('rho', origin='native')
#p.set_center((Radius,Radius/4))
#p.set_width(2.0*Radius,Radius/2)
#p.set_width(12,12)
p.save('phi.png')

proj = ds.slice('theta', np.pi/4)
p = proj.to_pw('rho', origin='native')
#p.set_log('rho', False)
#p.set_log(('index', 'grid_level'), False)
#p.set_log(('index', 'grid_indices'), False)
#p = proj.to_pw('W_vel1', origin='native')
#p.set_center((25/2,25/2))
#p.set_width(25,25)

p.save('xy.png')
