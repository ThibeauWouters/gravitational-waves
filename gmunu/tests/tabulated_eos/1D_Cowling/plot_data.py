import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

import yt
import numpy as np

#define the parameters to transform the code unit to cgs unit.
ggrav = 6.673e-8
clite = 2.99792458e10
clite_g = 2.99792458e10

rho_gf = 1.61930347e-18
press_gf = 1.80171810e-39
eps_gf = 1.11265006e-21
time_gf = 2.03001708e05
mass_gf = 5.02765209e-34
length_gf = 6.77140812e-06
energy_gf = 5.59424238e-55
lum_gf = 2.7556091e-60

mev_to_erg = 1.60217733e-6
erg_to_mev = 6.24150636e5
amu_cgs = 1.66053873e-24
massn_cgs = 1.674927211e-24
amu_mev = 931.49432e0
kb_erg = 1.380658e-16
kb_mev = 8.61738568e-11
temp_mev_to_kelvin = 1.1604447522806e10
planck = 6.626176e-27
avo = 6.0221367e23
hbarc_mevcm = 1.97326966e-11


units_override = dict( 
length_unit=(1.0/length_gf, 'cm'),
time_unit=(1.0/time_gf, 's'),
mass_unit=(1.0/mass_gf, 'g'),
)

# Load the dataset
ds0 = yt.load('output0000.dat', geometry_override='spherical',units_override=units_override, unit_system='cgs')
ds = yt.load('output0400.dat', geometry_override='spherical',units_override=units_override, unit_system='cgs')

# cutting the x-axis through the y=0,z=0 
plotdata0 = ds0.ortho_ray( 0, (0, 0) )
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['r'])

# code unit to cgs unit.
#plotdata['rho'] = np.array(plotdata['rho']) / rho_gf
#print (plotdata['r'])

#print (ds.field_list)
#print (ds.derived_field_list)

# plot the data

plt.subplot(4,1,1)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['rho'][srt]),label='rho')
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata0['rho'][srt]),label='rho')
plt.ylabel('rho')
#plt.ylim(1E-12,1E-2)
plt.grid(True)
plt.yscale('log')

plt.subplot(4,1,2)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['press'][srt]),label='p')
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata0['press'][srt]),label='p')
plt.ylabel('p')
plt.grid(True)
plt.yscale('log')

plt.subplot(4,1,3)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['eps'][srt]),label='eps')
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata0['eps'][srt]),label='eps')
plt.ylabel('eps')
plt.grid(True)
plt.yscale('log')

plt.subplot(4,1,4)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['veloc1'][srt]), label='v')
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata0['veloc1'][srt]), label='v')
plt.ylabel('v1')
plt.grid(True)
#plt.yscale('log')

# Save the line plot
#plt.legend(loc='best')
#plt.show()
plt.savefig('hydro.png',papertype='a0')
plt.close()

plt.subplot(3,1,1)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['alp'][srt]),label='alp')
plt.ylabel('alp')
plt.grid(True)

plt.subplot(3,1,2)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['psi'][srt]),label='psi')
plt.ylabel('psi')
plt.grid(True)

plt.subplot(3,1,3)
plt.plot(np.array(plotdata['r'][srt]), np.array(plotdata['beta1'][srt]), label='beta1')
plt.ylabel('beta1')
plt.grid(True)

# Save the line plot
#plt.legend(loc='best')
#plt.show()
plt.savefig('metric.png',papertype='a0')
