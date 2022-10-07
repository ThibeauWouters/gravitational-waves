import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np

def log_Energy(field, data):
    return np.log(data['nu_E']/1.0)

# Load the dataset
ds = yt.load('output0037.dat', unit_system='code')

ds.print_stats()
#print(ds.field_list)
#print(ds.derived_field_list)

ds.add_field( ('amrvac','output_E'), function=log_Energy, sampling_type='cell')
#p = yt.SlicePlot(ds, 2, 'output_E')
#p.set_colorbar_label('output_E', 'log(E/Emax)')
p = yt.SlicePlot(ds, 2, 'nu_E')
p.set_log("nu_E", True)
#p.set_cmap(field="nu_E", cmap='hot')
#p.annotate_quiver('nu_mom1', 'nu_mom2', factor=50, scale=20., normalize=False,  plot_args={'clim':(1.0E-10, 1E2), 'pivot':'middle'})
#p.annotate_quiver('nu_mom1', 'nu_mom2', factor=50, normalize=False,  plot_args={'clim':(1.0E-10, 1E2), 'pivot':'middle'})
#p.annotate_cell_edges()
#p.annotate_grids()
p.save('results.png')
