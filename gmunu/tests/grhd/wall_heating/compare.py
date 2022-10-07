import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

from matplotlib import rc_context

import yt
import numpy as np

# alpha = 0 for Cartesian
# alpha = 1 for cylindrical
# alpha = 2 for spherical
alpha = 1

if alpha == 0:
    x1 = 'x'
    coord = 'Cartesian'
elif alpha == 1:
    x1 = 'r'
    coord = 'cylindrical'
else:
    x1 = 'r'
    coord = 'spherical'

rho_0 = 1.0E0
polytropic_index = 4./3.

gamma_in = 1.0E1
v_in = np.sqrt(1.0 - 1.0 / gamma_in**2)
v_sh = ( polytropic_index - 1.0 ) * gamma_in * np.abs(v_in) / (gamma_in+1.0)

def rho_sol_pt(x,t):
    if (np.abs(x) >= v_sh * t):
        return (1.0E0 + np.abs(v_in / x )* t)**alpha * rho_0
    else:
        rho1 = (1.0E0 + np.abs(v_in) / v_sh)**alpha * rho_0
        return rho1 * (gamma_in * polytropic_index + 1.0) / (polytropic_index - 1.0)

def D_sol_pt(x,t):
    if (np.abs(x) >= v_sh * t):
        return (1.0E0 + np.abs(v_in / x )* t)**alpha * rho_0 * gamma_in
    else:
        rho1 = (1.0E0 + np.abs(v_in) / v_sh)**alpha * rho_0
        return rho1 * (gamma_in * polytropic_index + 1.0) / (polytropic_index - 1.0)
def cons_D(gamma,rho):
    return gamma * rho

def gamma(field, data):
    return  np.sqrt( 1.0 + (data['W_vel1']**2) )
def veloc1(field, data):
    return data['W_vel1'] / data['gamma']

rho_sol = np.vectorize(rho_sol_pt)
D_sol = np.vectorize(D_sol_pt)

# Load the dataset
#ds = yt.load('output_res_20480000.dat', geometry_override=coord)
#ds = yt.load('output0000.dat', geometry_override=coord)
#ds = yt.load('output0010.dat', geometry_override=coord)
ds = yt.load('output0020.dat', geometry_override=coord)
T_final = float(ds.current_time)
ds.add_field( ('amrvac','gamma'), function=gamma, sampling_type='cell')
ds.add_field( ('amrvac','veloc1'), function=veloc1, sampling_type='cell')

# cutting the x-axis through the y=0,z=0 
plotdata = ds.ortho_ray( 0, (0, 0) )

# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt = np.argsort(plotdata[x1])

#print (ds.field_list)
#print (ds.derived_field_list)

# plot the data
plt.plot(np.array(plotdata[x1][srt]), np.array(D_sol(np.array(plotdata[x1][srt]),T_final)), color="black",linestyle = '-',linewidth = 3 ) 
plt.plot(np.array(plotdata[x1][srt]), np.array(cons_D(plotdata['gamma'][srt],plotdata['rho'][srt])), 'r+')
#plt.plot(np.array(plotdata[x1][srt]), np.array(rho_sol(np.array(plotdata[x1][srt]),T_final)), color="black",linestyle = '-',linewidth = 3 ) 
#plt.plot(np.array(plotdata[x1][srt]), np.array(plotdata['rho'][srt]), 'r+',label='rho')
#plt.plot(np.array(plotdata[x1][srt]), np.array(plotdata['press'][srt]), color="red", linestyle = ':', linewidth=2.5, label='press')
#plt.plot(np.array(plotdata[x1][srt]), np.array(plotdata['eps'][srt]), color="red", linestyle = ':', linewidth=2.5, label='eps')
#plt.plot(np.array(plotdata[x1][srt]), np.array(plotdata['cs2'][srt]), color="red", linestyle = ':', linewidth=2.5, label='cs2')
#plt.plot(np.array(plotdata[x1][srt]), np.array(plotdata['gamma'][srt]), color="red", linestyle = ':', linewidth=2.5, label='gamma')
#plt.plot(np.array(plotdata[x1][srt]), np.array(plotdata['veloc1'][srt]), color="red", linestyle = ':', linewidth=2.5, label='v1', marker='.')

plt.xlabel('$R$')
plt.ylabel('$D$')
plt.xlim(0,1)
#plt.ylim(-1e-2,1e-2)
#plt.legend(loc='best')
plt.grid(True)
#plt.show()
plt.savefig('results.png')

