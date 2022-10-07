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
ds = yt.load('output_res_1024_limiter_woodward_0000.dat', unit_system='code')

ds.print_stats()
#print(ds.field_list)
#print(ds.derived_field_list)

# create a slice parallel to the yz-plane, in the middle of the x-range
p = yt.SlicePlot(ds,'z', 'rho')
#p.annotate_contour('rho')
#p.annotate_grids()
p.set_log('rho', False)

#p.annotate_mesh_lines(plot_args={'color':'black'})

# Save the line plot
#p.show()
p.save('results.png')
