import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot
from matplotlib import rc_context

import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import AxesGrid

def gamma(field, data):
    return  np.sqrt( 1.0 + np.square( data['r']*data['W_vel2'] ) ) 
    #return  np.sqrt( 1.0 + data['psi']**4*(data['W_vel1']**2 + data['r']**2*data['W_vel2']**2 + data['r']**2*np.sin(data['theta'])**2*data['W_vel3']**2) ) 
def sqrt_vphi2(field, data):
    return data['psi']**2*data['r']*np.sin(data['theta'])*data['W_vel2'] / data['gamma']

# Load the dataset
ds_initial = yt.load('output0000.dat', geometry_override='spherical', unit_system='code')
#ds_initial.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')

#print (ds_initial.field_list)
#print (ds_initial.derived_field_list)

ds_final = yt.load('output0001.dat', geometry_override='spherical', unit_system='code')


# cutting the x-axis through the y=0,z=0 
#theta = 0.
#theta = np.pi/4.0
theta = 1.57

plotdata_initial = ds_initial.ortho_ray( 0, (theta, 0) )
plotdata_final = ds_final.ortho_ray( 0, (theta, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata_initial['r'])

gamma_initial = np.array(plotdata_initial['W_vel1']**2)
gamma_initial = gamma_initial + np.array((plotdata_initial['r'])**2 * plotdata_initial['W_vel2']**2 + (plotdata_initial['r'])**2*np.sin(plotdata_initial['theta'])**2*plotdata_initial['W_vel3']**2)
gamma_initial = gamma_initial * np.array(plotdata_initial['psi']**4)
gamma_initial = np.sqrt(1.0 + gamma_initial)

sqrt_vphi2_initial = np.array(plotdata_initial['psi']**2*plotdata_initial['r']*np.sin(plotdata_initial['theta'])*plotdata_initial['W_vel3'])
sqrt_vphi2_initial = sqrt_vphi2_initial / gamma_initial


gamma_final = np.array(plotdata_final['W_vel1']**2)
gamma_final = gamma_final + np.array((plotdata_final['r'])**2 * plotdata_final['W_vel2']**2 + (plotdata_final['r'])**2*np.sin(plotdata_final['theta'])**2*plotdata_final['W_vel3']**2)
gamma_final = gamma_final * np.array(plotdata_final['psi']**4)
gamma_final = np.sqrt(1.0 + gamma_final)

sqrt_vphi2_final = np.array(plotdata_final['psi']**2*plotdata_final['r']*np.sin(plotdata_final['theta'])*plotdata_final['W_vel3'])
sqrt_vphi2_final = sqrt_vphi2_final / gamma_final

# rescale vphi2
vphi2_max = max(sqrt_vphi2_initial)
sqrt_vphi2_initial = sqrt_vphi2_initial / vphi2_max
sqrt_vphi2_final = sqrt_vphi2_final / vphi2_max

#plt.plot(np.array(plotdata_final['r'][srt]), np.array((plotdata_final['rho'][srt]) / (plotdata_initial['rho'][0]) ), 'r.', label='$\\rho / \\rho_c$' )
#plt.plot(np.array(plotdata_initial['r'][srt]), np.array((plotdata_initial['rho'][srt]) / (plotdata_initial['rho'][0]) ), color='black', linestyle='-', linewidth = 1)
#plt.yscale("log")

#plt.plot(np.array(plotdata_final['r'][srt]), gamma_final, 'r.', label='$ \\widetilde{ \\sqrt{v^\\phi v_\\phi} } $' )
#plt.plot(np.array(plotdata_initial['r'][srt]), gamma_initial, color='green', linestyle='-', linewidth = 1)

plt.plot(np.array(plotdata_final['r'][srt]), np.array(plotdata_final['beta3'][srt]), 'r.', label='$ \\widetilde{ \\sqrt{v^\\phi v_\\phi} } $' )
plt.plot(np.array(plotdata_initial['r'][srt]), np.array(plotdata_initial['beta3'][srt]), color='green', linestyle='-', linewidth = 1)
plt.plot(np.array(plotdata_final['r'][srt]), np.array(plotdata_final['W_vel3'][srt]), 'b.', label='$ \\widetilde{ \\sqrt{v^\\phi v_\\phi} } $' )
plt.plot(np.array(plotdata_initial['r'][srt]), np.array(plotdata_initial['W_vel3'][srt]), color='black', linestyle='-', linewidth = 1)
#plt.plot(np.array(plotdata_final['r'][srt]), sqrt_vphi2_final, 'b.', label='$ \\widetilde{ \\sqrt{v^\\phi v_\\phi} } $' )
#plt.plot(np.array(plotdata_initial['r'][srt]), sqrt_vphi2_initial, color='black', linestyle='-', linewidth = 1)

#plt.plot(np.array(plotdata_final['r'][srt]), np.array((plotdata_final['Bvec3'][srt]) / max(plotdata_initial['Bvec3']) ), 'g.', label='$B^\\phi /B^\\phi_c $' )
#plt.plot(np.array(plotdata_initial['r'][srt]), np.array((plotdata_initial['Bvec3'][srt]) / max(plotdata_initial['Bvec3']) ), color='black', linestyle='-', linewidth = 1)

plt.xlim(0,15)


# Save the line plot
plt.legend(loc='best')
plt.grid(True)

with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('slices.png')

