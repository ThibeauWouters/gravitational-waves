import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

from matplotlib import rc_context

import yt
import numpy as np

A = 0.2
v0 = 0.2
theta = 0.0

def rho_sol_pt(x,y,t):
    return 1.0 + A * np.sin(2.0 * np.pi * ( ( x * np.cos(theta) + y * np.sin(theta) ) - ( v0 * np.cos(theta) ) * t ) )

rho_sol = np.vectorize(rho_sol_pt)

def gamma(field, data):
    return  np.sqrt( 1.0 + (data['W_vel1']**2) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']

# Load the dataset
#ds = yt.load('output0020.dat')
ds = yt.load('output_res_64_limiter_weno5_0001.dat')
T_final = float(ds.current_time)
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')

# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['x'])

#print (ds.field_list)
#print (ds.derived_field_list)
# plot the data
#plt.plot(np.array(plotdata['x'][srt]), np.array(rho_sol(np.array(plotdata['x'][srt]),0.0,T_final)), color="black",linestyle = '-',linewidth = 3 ) 
plt.plot(np.array(plotdata['x'][srt]), np.array(rho_sol(np.array(plotdata['x'][srt]),0.0,T_final)) - np.array(plotdata['rho'][srt]), 'r+')

plt.xlabel('$x$')
plt.ylabel('$\\rho$')
plt.xlim(0,1)
#plt.ylim(-1e-2,1e-2)
#plt.legend(loc='best')
plt.grid(True)
#plt.show()
plt.savefig('results.png')

