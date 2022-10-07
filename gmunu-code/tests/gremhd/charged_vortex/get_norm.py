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

for runname in glob.glob('output*woodward*.log'):
    print(os.path.splitext(runname)[0].split("_"))
    nx_p = int(os.path.splitext(runname)[0].split("_")[2])
    nx = np.append(nx,[nx_p])
    T, L1_norm_p = np.loadtxt(runname, unpack=True ,dtype='float',usecols=(1,5),skiprows=1)
    print(nx_p,L1_norm_p[1])  
    L1_norm = np.append(L1_norm,[L1_norm_p[1]])

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
    print(i,rate_p)  

df = pd.DataFrame({'nx' : nx[srt], 
                   'L1_norm' : L1_norm[srt],
                   'rate' : rate}  )
df['nx'] = df['nx'].astype('int')
df['L1_norm'] = df['L1_norm'].astype('float')
df['rate'] = df['rate'].astype('float')

df = df[['nx','L1_norm','rate']]

df.to_csv("norms_"+ os.path.splitext(runname)[0].split("_")[4] +".csv", index = False, sep = '\t')

