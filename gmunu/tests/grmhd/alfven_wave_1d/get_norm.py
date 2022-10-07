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

nx = np.array([])
L1_norm = np.array([])
rate = np.array([])

for runname in glob.glob('output*.log'):
    print(os.path.splitext(runname)[0].split("_"))
    nx_p = int(os.path.splitext(runname)[0].split("_")[2])
    nx = np.append(nx,[nx_p])
    # Load the dataset
    ds_i = yt.load(os.path.splitext(runname)[0]+'0000.dat')
    ds_f = yt.load(os.path.splitext(runname)[0]+'0001.dat')
    # cutting the x-axis through the y=0,z=0 
    plotdata_i = ds_i.ortho_ray( 0, (0, 0) )
    plotdata_f = ds_f.ortho_ray( 0, (0, 0) )
    # Sort the ray values by 'x' so there are no discontinuities in the line plot
    srt_i = np.argsort(plotdata_i['x'])
    srt_f = np.argsort(plotdata_f['x'])
    Bz_i = np.array(plotdata_i['Bvec3'][srt_i])
    Bz_f = np.array(plotdata_f['Bvec3'][srt_f])
    Bz_diff = np.abs(Bz_f - Bz_i)
    Vol = np.array(plotdata_i['cell_volume'][srt_i])
    L1_norm_p = np.sum(Bz_diff*Vol)
    L1_norm_p = L1_norm_p / np.sum(Vol)
    L1_norm = np.append(L1_norm,[L1_norm_p])
    #print(nx_p,L1_norm_p)  

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
    print(i, ind, nx[i],L1_norm[i],rate[i])  
    

df = pd.DataFrame({'nx' : nx[srt], 
                   'L1_norm' : L1_norm[srt],
                   'rate' : rate}  )

df['nx'] = df['nx'].astype('int')
df['L1_norm'] = df['L1_norm'].astype('float')
df['rate'] = df['rate'].astype('float')

df = df[['nx','L1_norm','rate']]

df.to_csv("norms_"+ os.path.splitext(runname)[0].split("_")[4] +".csv", index = False, sep = '\t')
