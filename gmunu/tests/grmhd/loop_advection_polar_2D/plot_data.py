import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np

units_override = dict( 
length_unit=(1.0, 'cm'),
time_unit=(1.0, 's'),
mass_unit=(1.0, 'g'),
)

# Load the dataset
ds = yt.load('output0005.dat', geometry_override='polar', unit_system='code')
#ds.add_field( ('amrvac','B2'), function=B2, sampling_type='cell')

p = yt.SlicePlot(ds,'z', 'B2')
#p.set_zlim("mag_press", 1E-12, 1E-7)
#p.set_log('mag_press', False)
#p = yt.SlicePlot(ds,'z', 'B2')
p.set_log('B2', False)

#p = yt.SlicePlot(ds,'z', 'Bvec1')
#p.set_log('Bvec1', False)
#p = yt.SlicePlot(ds,'z', 'Bvec2')
#p.set_log('Bvec2', False)

#p.set_center((Radius/2.0,Radius/2.0))
#p.set_width(Radius,Radius)

# Save the line plot
#p.show()
p.save('results.png')

