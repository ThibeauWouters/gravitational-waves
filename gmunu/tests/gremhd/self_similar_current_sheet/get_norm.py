import yt
import numpy as np
import pandas as pd
import glob
import os
from scipy.interpolate import interp1d
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot
from matplotlib import rc_context

highest_res = 8192
ds_f = yt.load('output_res_8192_limiter_woodward_0001.dat')
plotdata_f = ds_f.ortho_ray( 0, (0, 0) )
# Sort the ray values by 'x' so there are no discontinuities in the line plot
srt_f = np.argsort(plotdata_f['x'])
By_ref = interp1d(plotdata_f['x'][srt_f], plotdata_f['Bvec2'][srt_f], kind = 'cubic')

nx = np.array([])
L1_norm = np.array([])
rate = np.array([])

for runname in glob.glob('output*.log'):
    print(os.path.splitext(runname)[0].split("_"))
    nx_p = int(os.path.splitext(runname)[0].split("_")[2])
    if nx_p == highest_res:
        continue
    nx = np.append(nx,[nx_p])
    # Load the dataset
    ds_f = yt.load(os.path.splitext(runname)[0]+'0001.dat')
    T_final = float(ds_f.current_time)
    # cutting the x-axis through the y=0,z=0 
    plotdata_f = ds_f.ortho_ray( 0, (0, 0) )
    # Sort the ray values by 'x' so there are no discontinuities in the line plot
    srt_f = np.argsort(plotdata_f['x'])
    By_f = np.array(plotdata_f['Bvec2'][srt_f])
    #rho_sol = np.array(rho_sol(np.array(plotdata_f['x'][srt_f]),0.0,T_final))
    By_solution = By_ref(np.array(plotdata_f['x'][srt_f]))
    By_diff = np.abs(By_solution - By_f) / np.max(By_solution)
    L1_norm_p = np.sum(By_diff)
    L1_norm_p = L1_norm_p / nx_p
    L1_norm = np.append(L1_norm,[L1_norm_p])
    print(nx_p,L1_norm_p)  

#simname=os.path.splitext()
#print(simname)

# Sort the values by nx
srt = np.argsort(nx)

for i,ind in enumerate(nx[srt]):    
    if (i>0):
        rate_p = np.log( L1_norm[srt[i-1]] / L1_norm[srt[i]] ) / np.log(2.0)
    else:
        rate_p = 0.0
    rate = np.append( rate,[rate_p] )

df = pd.DataFrame({'nx' : nx[srt], 
                   'L1_norm' : L1_norm[srt],
                   'rate' : rate}  )
df['nx'] = df['nx'].astype('int')
df['L1_norm'] = df['L1_norm'].astype('float')
df['rate'] = df['rate'].astype('float')

df = df[['nx','L1_norm','rate']]

df.to_csv("norms_"+ os.path.splitext(runname)[0].split("_")[4] +".csv", index = False, sep = '\t')

