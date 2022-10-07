import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np


debug = False
debug = True

units_override = dict( 
length_unit=(1.0, 'cm'),
time_unit=(1.0, 's'),
mass_unit=(1.0, 'g'),
)

# Load the dataset
ds = yt.load('output0000.dat', geometry_override='spherical', unit_system='code')

ds.print_stats()
print(ds.field_list)
print(ds.derived_field_list)

Radius = ds.domain_width[0]
p = yt.SlicePlot(ds,2, 'density')
p.set_cmap(field="density", cmap='hot')
p.set_center((Radius/2.0,Radius/2.0))
p.set_width(Radius,Radius)
#p.set_log('density', True)
#p.annotate_cell_edges()

# Save the line plot
#p.show()
p.save('hydro.png')

if (False):
#if (True):
  p = yt.SlicePlot(ds,2, 'W_vel3')
  p.set_log('W_vel3', False)
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('W_vel3.png')
  p = yt.SlicePlot(ds,2, 'W_vel2')
  p.set_log('W_vel2', False)
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('W_vel2.png')
  p = yt.SlicePlot(ds,2, 'W_vel1')
  p.set_log('W_vel1', False)
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('W_vel1.png')

if (debug):
  p = yt.SlicePlot(ds,2, 'psi')
  p.set_log('psi', False)
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('psi.png')
