import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import scipy.signal as signal
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d

from matplotlib import rc_context
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':15})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



filepath = "../"


# read data
#it global_time dt rho_max alp_min M_rest c1_real c1_img c2_real c2_img c3_real c3_img c4_real c4_img rho22_real rho22_img J_rot T_rot I11_dot I12_dot I13_dot I22_dot I23_dot I33_dot
data_time, data_M_rest, data_c1_real, data_c1_img, data_c2_real, data_c2_img, data_c3_real, data_c3_img, data_c4_real, data_c4_img = np.loadtxt(filepath+"output.log", unpack=True ,dtype='float',usecols=(1,5,6,7,8,9,10,11,12,13),skiprows=1)

df_TD = pd.read_csv(filepath+'output.log', index_col = False, delim_whitespace=True, usecols = ['global_time','M_rest','c1_real','c1_img','c2_real','c2_img','c3_real','c3_img','c4_real','c4_img'])
df_TD = df_TD[['global_time','M_rest','c1_real','c1_img','c2_real','c2_img','c3_real','c3_img','c4_real','c4_img']]
data_time = df_TD['global_time'].to_numpy()
data_M_rest = df_TD['M_rest']
data_c1_real = df_TD['c1_real']
data_c2_real = df_TD['c2_real']
data_c3_real = df_TD['c3_real']
data_c4_real = df_TD['c4_real']
data_c1_img = df_TD['c1_img']
data_c2_img = df_TD['c2_img']
data_c3_img = df_TD['c3_img']
data_c4_img = df_TD['c4_img']

M_rest_interp = interp1d(data_time, data_M_rest, kind = 'cubic')
c1_real_interp = interp1d(data_time, data_c1_real, kind = 'cubic')
c1_img_interp = interp1d(data_time, data_c1_img, kind = 'cubic')

c2_real_interp = interp1d(data_time, data_c2_real, kind = 'cubic')
c2_img_interp = interp1d(data_time, data_c2_img, kind = 'cubic')

c3_real_interp = interp1d(data_time, data_c3_real, kind = 'cubic')
c3_img_interp = interp1d(data_time, data_c3_img, kind = 'cubic')

c4_real_interp = interp1d(data_time, data_c4_real, kind = 'cubic')
c4_img_interp = interp1d(data_time, data_c4_img, kind = 'cubic')

t_min = data_time[1]
t_max = data_time[-1]
t = np.linspace(t_min,t_max,len(data_time))
dt = t[1]-t[0]


M_rest = M_rest_interp(t)
c1_real = c1_real_interp(t)
c1_img = c1_img_interp(t)

c2_real = c2_real_interp(t)
c2_img = c2_img_interp(t)

c3_real = c3_real_interp(t)
c3_img = c3_img_interp(t)

c4_real = c4_real_interp(t)
c4_img = c4_img_interp(t)

c1 = np.sqrt(c1_real**2 + c1_img**2)
c2 = np.sqrt(c2_real**2 + c2_img**2)
c3 = np.sqrt(c3_real**2 + c3_img**2)
c4 = np.sqrt(c4_real**2 + c4_img**2)

c1_norm = c1/M_rest
c2_norm = c2/M_rest
c3_norm = c3/M_rest
c4_norm = c4/M_rest

t = t / 2.03001708E05 # code unit to sec 
dt = dt / 2.03001708E05 # code unit to sec 
# data frame for GW in time domain, in code unit
df_TD = pd.DataFrame( {'t': t,'M_rest': M_rest, 
            'c1_real': c1_real, 'c1_img': c1_img, 'c1': c1, 'c1_norm': c1_norm,      
            'c2_real': c2_real, 'c2_img': c2_img, 'c2': c2, 'c2_norm': c2_norm,
            'c3_real': c3_real, 'c3_img': c3_img, 'c3': c3, 'c3_norm': c3_norm,
            'c4_real': c4_real, 'c4_img': c4_img, 'c4': c4, 'c4_norm': c4_norm                   
             } )
df_TD = df_TD[['t','M_rest','c1_real','c1_img','c1','c1_norm','c2_real','c2_img','c2','c2_norm','c3_real','c3_img','c3','c3_norm','c4_real','c4_img','c4','c4_norm',]]
df_TD.to_csv('cm_TD_xy.csv', index = False, sep = '\t')


# next is to do fft
window = signal.tukey(len(t), alpha = 0.5)
fft_c2 = fft(c2) * 2.0 / t.size
freq = np.linspace(start=0.0, stop=1.0/(2.0 * dt), num=int(t.size/2) )
psd_c2 = np.abs(fft_c2)**2
df_FD = pd.DataFrame( {'f': freq, 'fft_c2': np.abs(fft_c2[:t.size//2]), 'psd_c2': psd_c2[:t.size//2]} )
df_FD = df_FD[['f','fft_c2','psd_c2']]
df_FD.to_csv('c2_FD_xy.csv', index = False, sep = '\t')

# Save the line plot
plt.plot(t,c1_norm,label='$m=1$',color='blue')
plt.plot(t,c2_norm,label='$m=2$',color='orange')
plt.plot(t,c3_norm,label='$m=3$',color='green')
plt.plot(t,c4_norm,label='$m=4$',color='red')
plt.yscale('log')

plt.legend(loc='best')
plt.xlabel('$t$')
plt.ylabel('$|C_m|/C_0$')
plt.xlim(t[0],t[-1])
#plt.ylim(-1e-2,1e-2)
plt.grid(True)
#plt.show()
plt.savefig('cm_t.png')
