import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot
import glob

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

data_files = glob.glob('./output*.dat')

ds = yt.load('output0000.dat')
# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )
# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['x'])


time = np.array([0.0])
datay = np.array(plotdata['x'][srt])
dataz = [np.array(plotdata['Bvec2'][srt])]

for i, fn in enumerate(data_files):
   ds = yt.load(fn, unit_system='code')
   ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
   ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')
   ds.add_field( ('amrvac','veloc2'), function=veloc2, sampling_type='cell')

   # cutting the x-axis through the y=0,z=0 
   plotdata = ds.ortho_ray( 0, (0, 0) )
   # Sort the ray values by 'x' so there are no discontinuities in the line plot
   srt = np.argsort(plotdata['x'])
   time = np.append(time,[ds.current_time])
   dataz = np.append(dataz,[np.array(plotdata['Bvec2'][srt])],axis=0)

dataz = np.transpose(dataz)
X, Y = np.meshgrid(time, datay)

#plt.grid(True)
#plt.legend(loc='best')
#plt.savefig('results.png')

#print ( time )
print ( dataz.shape )
#print ( dataz )

plt.pcolor(X,Y,dataz, vmin=-1.0, vmax=1.0)
plt.savefig('results.png')
