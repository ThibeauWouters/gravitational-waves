import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np

# Load the dataset
ds = yt.load('output0082.dat', geometry_override='cylindrical', unit_system='code')

ds.print_stats()
#print(ds.field_list)
#print(ds.derived_field_list)

# create a slice parallel to the yz-plane, in the middle of the x-range
#p = yt.SlicePlot(ds, 2, 'W_vel2')
p = yt.SlicePlot(ds, 2, 'rho')
#p.set_cmap(field="rho", cmap='hot')
#p.annotate_contour('rho')
#p.annotate_grids()
#p.set_log('rho', False)
p.set_log('W_vel2', False)

#p.annotate_mesh_lines(plot_args={'color':'black'})

# Save the line plot
#p.show()
p.save('results.png')
