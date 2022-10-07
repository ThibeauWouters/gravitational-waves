import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

import yt
import numpy as np



# Load the dataset
ds = yt.load('output0900.dat', geometry_override='spherical', unit_system='code')

# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['r'])

#print (ds.field_list)
#print (ds.derived_field_list)


# plot the data
#plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['alp'][srt]))
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['rho'][srt]))
plt.xlim(0,50)
plt.yscale('log')


# Save the line plot
plt.legend(loc='best')
plt.grid(True)
#plt.show()
plt.savefig('results.png')
