import matplotlib
#matplotlib.use('Agg')
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

if (iprob != 0):
   # plot the data
   fig, axs = plt.subplots(2,3)
   
   x, y = np.loadtxt("./solutions/balsara/rho.txt", unpack=True ,dtype='float',delimiter=',')
   axs[0,0].plot(x , y, color="black",linestyle = '-',linewidth = 1 ) 
   axs[0,0].plot(np.array(plotdata['x'][srt]), np.array(plotdata['rho'][srt]),'*',lw=0.1)
   axs[0,0].set_title('$\\rho$')
   #axs[0,0].set_ylim(0.0, 1.2)
   axs[0,0].set(xlabel='$x$') #, ylabel='$\\rho$')
   
   x, y = np.loadtxt("./solutions/balsara/p.txt", unpack=True ,dtype='float',delimiter=',')
   axs[1,0].plot(x , y, color="black",linestyle = '-',linewidth = 1 ) 
   axs[1,0].plot(np.array(plotdata['x'][srt]), np.array(plotdata['press'][srt]),'*',lw=0.1)
   axs[1,0].set_title('$P$')
   #axs[1,0].set_ylim(0.0, 1.2)
   axs[1,0].set(xlabel='$x$') #, ylabel='$P$')
   
   x, y = np.loadtxt("./solutions/balsara/vx.txt", unpack=True ,dtype='float',delimiter=',')
   axs[0,1].plot(x , y, color="black",linestyle = '-',linewidth = 1 ) 
   axs[0,1].plot(np.array(plotdata['x'][srt]), np.array(plotdata['veloc1'][srt]),'*',lw=0.1)
   axs[0,1].set_title('$v^x$')
   axs[0,1].set(xlabel='$x$') #, ylabel='$v^x$')
   
   x, y = np.loadtxt("./solutions/balsara/vy.txt", unpack=True ,dtype='float',delimiter=',')
   axs[1,1].plot(x , y, color="black",linestyle = '-',linewidth = 1 ) 
   axs[1,1].plot(np.array(plotdata['x'][srt]), np.array(plotdata['veloc2'][srt]),'*',lw=0.1)
   axs[1,1].set_title('$v^y$')
   axs[1,1].set(xlabel='$x$') #, ylabel='$v^y$')
   
   x, y = np.loadtxt("./solutions/balsara/By.txt", unpack=True ,dtype='float',delimiter=',')
   axs[0,2].plot(x , y, color="black",linestyle = '-',linewidth = 1 ) 
   axs[0,2].plot(np.array(plotdata['x'][srt]), np.array(plotdata['Bvec2'][srt]),'*',lw=0.1)
   axs[0,2].set_title('$B^y$')
   axs[0,2].set(xlabel='$x$') #, ylabel='$B^y$')
   
   x, y = np.loadtxt("./solutions/balsara/gamma.txt", unpack=True ,dtype='float',delimiter=',')
   axs[1,2].plot(x , y, color="black",linestyle = '-',linewidth = 1 ) 
   axs[1,2].plot(np.array(plotdata['x'][srt]), np.array(plotdata['gamma'][srt]),'*',lw=0.1)
   axs[1,2].set_title('$\\Gamma$')
   axs[1,2].set_ylim(0.9, 1.5)
   axs[1,2].set(xlabel='$x$') #, ylabel='$B^y$')
   
   fig.tight_layout()
   
   # Save the line plot
   #plt.legend(loc='best')
   #plt.grid(True)
   #plt.show()
   #plt.savefig('results.png', bbox_inches='tight')
   with rc_context({'mathtext.fontset': 'stix'}):
     plt.savefig('results.png')

else:
   # plot the data
   plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['rho'][srt]/10), label='rho/10')
   plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['press'][srt]/20), label='P/20')
   plt.plot(np.array(plotdata['x'][srt]), np.array(plotdata['veloc1'][srt]), label='v')
   
   
   # Save the line plot
   plt.legend(loc='best')
   plt.grid(True)
   #plt.show()
   plt.savefig('results.png')
