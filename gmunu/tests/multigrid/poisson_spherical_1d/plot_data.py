import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

import yt
import numpy as np

rho0 = 1.

def phi_solution(r):
  if (r<=1.):
    return np.pi * rho0 * (r**4/5. - 2./3.*r**2 + 1.)
  else:
    return 8.*np.pi / (15.*r)



# Load the dataset
ds = yt.load('output0000.dat', geometry_override='spherical', unit_system='code')

# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['r'])

# code unit to cgs unit.
#plotdata['rho'] = np.array(plotdata['rho']) / rho_gf
#print (plotdata['r'])

#print (ds.field_list)
#print (ds.derived_field_list)

# plot the data
analytic_solution=  (np.array([phi_solution(r) for r in np.array(plotdata['r'][srt])]))

plt.subplot(2,1,1)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['rho'][srt]),label='rho')
plt.ylabel('rho')
#plt.ylim(1E-12,1E-2)
plt.grid(True)
#plt.yscale('log')
plt.subplot(2,1,2)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['psi'][srt]),'+')
plt.plot(np.array(plotdata['r'][srt]), analytic_solution,label='Phi')
#plt.plot(phi_solution(plotdata['r'][srt]),label='Phi')
#plt.ylim(0,4)
plt.ylabel('Phi')
plt.grid(True)

plt.savefig('results.png',papertype='a0')
plt.close()


