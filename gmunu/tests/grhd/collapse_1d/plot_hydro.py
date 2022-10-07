import matplotlib.pyplot as plt
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
length_gf = 6.77140812e-06*1.0e5 # to km
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
length_unit=(1.0/length_gf, 'km'),
time_unit=(1.0/time_gf, 's'),
mass_unit=(1.0/mass_gf, 'g'),
)



# Load the dataset
ds = yt.load('output0002.dat', geometry_override='spherical', unit_system='code')

# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata['r'])

#print (ds.field_list)
#print (ds.derived_field_list)


fig, axs = plt.subplots(2,1,sharex=True, gridspec_kw={'height_ratios': [3, 1]})
fig.subplots_adjust(hspace=0)

# plot data
axs[0].plot(np.array(plotdata['r'][srt])/length_gf, np.array(plotdata['rho'][srt])/ rho_gf)

axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_ylabel('$\\rho\ \\rm{[g/cm^3]}$')
axs[0].grid(True)

# plot grid level
axs[1].plot(np.array(plotdata['r'][srt]), np.array(plotdata['grid_level'][srt]), label='level')
axs[1].set_xlabel('$r\ \\rm{[km]}$')
axs[1].set_ylabel('$\\rm{grid\ level}$')
axs[1].set_xscale('log')
axs[1].grid(True)


# Save the line plot
#plt.show()
plt.savefig('rho.pdf', bbox_inches="tight")

