import matplotlib.pyplot as plt
import yt
import numpy as np
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':15})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

def E_r(field, data):
    return  ( data['psi']**2*data['Evec1'] )
def E_z(field, data):
    return  ( data['psi']**2*data['Evec2'] )
def E_phi(field, data):
    return  ( data['psi']**2*data['r']*data['Evec3'] )
def B_r(field, data):
    return  ( data['psi']**2*data['Bvec1'] )
def B_z(field, data):
    return  ( data['psi']**2*data['Bvec2'] )
def B_phi(field, data):
    return  ( data['psi']**2 * data['r'] * data['Bvec3'] )

# Load the dataset
ds_initial = yt.load('output0000.dat', geometry_override='cylindrical', unit_system='code')
ds_initial.add_field( ('amrvac','E_r'), function=E_r, sampling_type='cell')
ds_initial.add_field( ('amrvac','E_z'), function=E_z, sampling_type='cell')
ds_initial.add_field( ('amrvac','E_phi'), function=E_phi, sampling_type='cell', units='code_length')
ds_initial.add_field( ('amrvac','B_r'), function=B_r, sampling_type='cell')
ds_initial.add_field( ('amrvac','B_z'), function=B_z, sampling_type='cell')
ds_initial.add_field( ('amrvac','B_phi'), function=B_phi, sampling_type='cell', units='code_length')

ds_final = yt.load('output0007.dat', geometry_override='cylindrical', unit_system='code')
ds_final.add_field( ('amrvac','E_r'), function=E_r, sampling_type='cell')
ds_final.add_field( ('amrvac','E_z'), function=E_z, sampling_type='cell')
ds_final.add_field( ('amrvac','E_phi'), function=E_phi, sampling_type='cell', units='code_length')
ds_final.add_field( ('amrvac','B_r'), function=B_r, sampling_type='cell')
ds_final.add_field( ('amrvac','B_z'), function=B_z, sampling_type='cell')
ds_final.add_field( ('amrvac','B_phi'), function=B_phi, sampling_type='cell', units='code_length')

fig, axs = plt.subplots(7,3, sharex='col',sharey='row', figsize=(20,15))

# cutting the x-axis through the y=0,z=0 
z_cut = 5.

# r, theta, z
plotdata_initial = ds_initial.ortho_ray('r', (0, z_cut) )
plotdata_final = ds_final.ortho_ray( 'r', (0, z_cut) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt_initial = np.argsort(plotdata_initial['r'])
srt_final = np.argsort(plotdata_final['r'])

# plot density
axs[0,0].set_title('$z = 5$')
axs[0,0].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['rho'][srt_final])),'r.',markersize=6)
axs[0,0].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['rho'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[0,0].set_xlim(0, 15 )
#axs[0,0].set_ylim(1E-10, 1E-2)
axs[0,0].set_yscale('log')
axs[0,0].set( ylabel='$\\rho $')
axs[0,0].grid(True)

# plot E_r
axs[1,0].set( ylabel='$ { \\sqrt{E^r E_r} } $')
axs[1,0].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_r'][srt_final])),'r.',markersize=6)
axs[1,0].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_r'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[1,0].grid(True)
# plot E_z
axs[2,0].set( ylabel='$ { \\sqrt{E^z E_z} } $')
axs[2,0].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_z'][srt_final])),'r.',markersize=6)
axs[2,0].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_z'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[2,0].grid(True)
# plot E_phi
axs[3,0].set( ylabel='$ { \\sqrt{E^\\varphi E_\\varphi} } $')
axs[3,0].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_phi'][srt_final])),'r.',markersize=6)
axs[3,0].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_phi'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[3,0].grid(True)

# plot B_r
axs[4,0].set( ylabel='$ { \\sqrt{B^r B_r} } $')
axs[4,0].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_r'][srt_final])),'r.',markersize=6)
axs[4,0].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_r'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[4,0].grid(True)
# plot B_z
axs[5,0].set( ylabel='$ { \\sqrt{B^z B_z} } $')
axs[5,0].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_z'][srt_final])),'r.',markersize=6)
axs[5,0].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_z'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[5,0].grid(True)
# plot B_phi
axs[6,0].set( ylabel='$ { \\sqrt{B^\\varphi B_\\varphi} } $')
axs[6,0].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_phi'][srt_final])),'r.',markersize=6)
axs[6,0].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_phi'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[6,0].grid(True)

axs[6,0].set( xlabel='$ R $')


z_cut = 3.

# r, theta, z
plotdata_initial = ds_initial.ortho_ray('r', (0, z_cut) )
plotdata_final = ds_final.ortho_ray( 'r', (0, z_cut) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt_initial = np.argsort(plotdata_initial['r'])
srt_final = np.argsort(plotdata_final['r'])

# plot density
axs[0,1].set_title('$z = 3$')
axs[0,1].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['rho'][srt_final])),'r.',markersize=6)
axs[0,1].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['rho'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[0,1].set_xlim(0, 15 )
#axs[0,1].set_ylim(1E-10, 1E-2)
axs[0,1].set_yscale('log')
axs[0,1].grid(True)

# plot E_r
axs[1,1].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_r'][srt_final])),'r.',markersize=6)
axs[1,1].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_r'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[1,1].grid(True)
# plot E_z
axs[2,1].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_z'][srt_final])),'r.',markersize=6)
axs[2,1].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_z'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[2,1].grid(True)
# plot E_phi
axs[3,1].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_phi'][srt_final])),'r.',markersize=6)
axs[3,1].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_phi'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[3,1].grid(True)

# plot B_r
axs[4,1].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_r'][srt_final])),'r.',markersize=6)
axs[4,1].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_r'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[4,1].grid(True)
# plot B_z
axs[5,1].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_z'][srt_final])),'r.',markersize=6)
axs[5,1].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_z'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[5,1].grid(True)
# plot B_phi
axs[6,1].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_phi'][srt_final])),'r.',markersize=6)
axs[6,1].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_phi'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[6,1].grid(True)

axs[6,1].set( xlabel='$ R $')

# cutting the x-axis through the y=0,z=0 
z_cut = 0.

# r, theta, z
plotdata_initial = ds_initial.ortho_ray('r', (0, z_cut) )
plotdata_final = ds_final.ortho_ray( 'r', (0, z_cut) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt_initial = np.argsort(plotdata_initial['r'])
srt_final = np.argsort(plotdata_final['r'])

# plot density
axs[0,2].set_title('$z = 0$')
axs[0,2].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['rho'][srt_final])),'r.',markersize=6)
axs[0,2].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['rho'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[0,2].set_xlim(0, 15 )
#axs[0,2].set_ylim(1E-10, 1E-2)
axs[0,2].set_yscale('log')
axs[0,2].grid(True)

# plot E_r
axs[1,2].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_r'][srt_final])),'r.',markersize=6)
axs[1,2].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_r'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[1,2].grid(True)
# plot E_z
axs[2,2].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_z'][srt_final])),'r.',markersize=6)
axs[2,2].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_z'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[2,2].grid(True)
# plot E_phi
axs[3,2].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['E_phi'][srt_final])),'r.',markersize=6)
axs[3,2].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['E_phi'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[3,2].grid(True)

# plot B_r
axs[4,2].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_r'][srt_final])),'r.',markersize=6)
axs[4,2].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_r'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[4,2].grid(True)
# plot B_z
axs[5,2].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_z'][srt_final])),'r.',markersize=6)
axs[5,2].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_z'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[5,2].grid(True)
# plot B_phi
axs[6,2].plot(np.array(plotdata_final['r'][srt_final]), np.array((plotdata_final['B_phi'][srt_final])),'r.',markersize=6)
axs[6,2].plot(np.array(plotdata_initial['r'][srt_initial]), np.array((plotdata_initial['B_phi'][srt_initial])), color='blue', linestyle='-', linewidth = 1)
axs[6,2].grid(True)
#axs[6,2].set_ylim(-1E-10, 1E-10)
#axs[6,2].set_yscale('log')

axs[6,2].set( xlabel='$ R $')


plt.subplots_adjust(hspace=0.1, wspace=0.1)
plt.savefig('./fig_MHD_NS_cylindrical_2D_slices.pdf', bbox_inches="tight")
