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

ds.print_stats()
#print(ds.field_list)
#print(ds.derived_field_list)

Radius = ds.domain_width[0]
p = yt.SlicePlot(ds,2, 'Bvec3')
p.set_cmap(field="Bvec3", cmap='hot')
p.set_center((Radius/2.0,Radius/2.0))
p.set_width(Radius,Radius)
#p.set_log('density', True)
#p.annotate_cell_edges()

# Save the line plot
#p.show()
p.save('results.png')

if (False):
#if (True):
  p = yt.SlicePlot(ds,2, 'veloc3')
  p.set_log('veloc3', False)
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('veloc3.png')
  p = yt.SlicePlot(ds,2, 'veloc2')
  p.set_log('veloc2', False)
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('veloc2.png')
  p = yt.SlicePlot(ds,2, 'veloc1')
  p.set_log('veloc1', False)
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('veloc1.png')

if (debug):
  p = yt.SlicePlot(ds,2, 'rho')
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('rho.png')
  p = yt.SlicePlot(ds,2, 'psi')
  p.set_log('psi', False)
  p.set_center((Radius/2.0,Radius/2.0))
  p.set_width(Radius,Radius)
  p.save('psi.png')
