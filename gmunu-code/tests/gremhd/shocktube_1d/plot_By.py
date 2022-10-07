import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

from matplotlib import rc_context

import yt
import numpy as np

iprob = 1

def gamma(field, data):
    return  np.sqrt( 1.0 + (data['W_vel1']**2 + data['W_vel2']**2) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']
def veloc2(field, data):
    return data['W_vel2'] / data['gamma']


# Load the dataset
ds = yt.load('shocktube0040.dat')
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')

# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['x'])

#print (ds.field_list)
#print (ds.derived_field_list)
x, y = np.loadtxt("./solutions/balsara/By.txt", unpack=True ,dtype='float',delimiter=',')
plt.plot(x , y, color="black",linestyle = '-',linewidth = 1 )
plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['Bvec2'][srt]),'*',lw=0.1)
#plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['W_vel1'][srt]),'*',lw=0.1, label='v1')
plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['Evec2'][srt]),'*',lw=0.1, label='E2')
plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['Evec3'][srt]),'*',lw=0.1, label='E3')

plt.grid(True)
plt.legend(loc='best')
plt.savefig('results.png')
