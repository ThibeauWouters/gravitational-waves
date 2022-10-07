import yt
import numpy as np
import pandas as pd
import glob
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot
from matplotlib import rc_context

A = 0.2
v0 = 0.2
theta = 0.0

def rho_sol_pt(x,t):
    return 1.0 + A * np.sin(2.0 * np.pi * ( ( x * np.cos(theta) + 0.0 * np.sin(theta) ) - ( v0 * np.cos(theta) ) * t ) )

rho_sol = np.vectorize(rho_sol_pt)

nx = np.array([])
L1_norm = np.array([])

for runname in glob.glob('output*.log'):
    print(os.path.splitext(runname)[0].split("_"))
    nx_p = int(os.path.splitext(runname)[0].split("_")[2])
    nx = np.append(nx,[nx_p])
    # Load the dataset
    ds_f = yt.load(os.path.splitext(runname)[0]+'0001.dat')
    T_final = float(ds_f.current_time)
    # cutting the x-axis through the y=0,z=0 
    plotdata_f = ds_f.ortho_ray( 0, (0, 0) )
    # Sort the ray values by 'x' so there are no discontinuities in the line plot
    srt_f = np.argsort(plotdata_f['x'])
    rho_f = np.array(plotdata_f['rho'][srt_f])
    #rho_sol = np.array(rho_sol(np.array(plotdata_f['x'][srt_f]),0.0,T_final))
    rho_solution = rho_sol(np.array(plotdata_f['x'][srt_f]),T_final)
    rho_diff = np.abs(rho_solution - rho_f)
    Vol = np.array(plotdata_f['cell_volume'][srt_f])
    L1_norm_p = np.sum(rho_diff*Vol)
    L1_norm_p = L1_norm_p / np.sum(Vol)
    L1_norm = np.append(L1_norm,[L1_norm_p])
    print(nx_p,L1_norm_p)  

#simname=os.path.splitext()
#print(simname)

# Sort the values by nx
srt = np.argsort(nx)

df = pd.DataFrame({'nx' : nx[srt], 
                   'L1_norm' : L1_norm[srt]}  )
df['nx'] = df['nx'].astype('int')

df = df[['nx','L1_norm']]

df.to_csv("norms_"+ os.path.splitext(runname)[0].split("_")[4] +".csv", index = False, sep = '\t')

