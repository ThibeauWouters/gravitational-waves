import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

import yt
import numpy as np

def gamma(field, data):
    return  np.sqrt( 1.0 + (data['W_vel1']**2 ) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']


# Load the dataset
ds = yt.load('shocktube0040.dat')
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')

# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['x'])

#print (ds.field_list)
#print (ds.derived_field_list)


# plot the data
plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['rho'][srt]/10), label='rho/10')
plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['press'][srt]/20), label='P/20')
plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['veloc1'][srt]), label='v')
#plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['grid_level'][srt]), label='level')


# Save the line plot
plt.legend(loc='best')
plt.grid(True)
#plt.show()
plt.savefig('results.pdf')
