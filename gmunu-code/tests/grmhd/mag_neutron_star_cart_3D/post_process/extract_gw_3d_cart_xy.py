import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import scipy.signal as signal
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d


filepath = "../"


# read data
df_TD = pd.read_csv(filepath+'output.log', index_col = False, delim_whitespace=True, usecols = ['global_time','I11_dot','I22_dot','I12_dot'])
df_TD = df_TD[['global_time','I11_dot','I22_dot','I12_dot']]
data_time = df_TD['global_time'].to_numpy()
data_I11_dot = df_TD['I11_dot'].to_numpy()
data_I22_dot = df_TD['I22_dot'].to_numpy()
data_I12_dot = df_TD['I12_dot'].to_numpy()

#data_time, data_I13_dot = np.loadtxt(filepath+"output.log", unpack=True ,dtype='float',usecols=(1,8),skiprows=1)

I11_dot_interp = interp1d(data_time, data_I11_dot, kind = 'cubic')
I12_dot_interp = interp1d(data_time, data_I12_dot, kind = 'cubic')
I22_dot_interp = interp1d(data_time, data_I22_dot, kind = 'cubic')


t_min = data_time[0]
t_max = data_time[-1]
#print(t_max/2.03001708E05)
#t_max = 2.03001708E03 * 2.
t = np.linspace(t_min,t_max,len(data_time))
dt = t[1]-t[0]

I11_dot = I11_dot_interp(t)
I22_dot = I22_dot_interp(t)
I11_ddot = np.gradient(I11_dot, dt)
I22_ddot = np.gradient(I22_dot, dt)

r = (6.7706*10**(-6))*3.08567758*10**26
#r=1
h_plus = (I11_ddot -  I22_ddot)/r # fixme: this is not real hp


I12_dot = I12_dot_interp(t)
I12_ddot = np.gradient(I12_dot, dt)
h_cross = 2*I12_ddot/r # fixme: this is not real hp


t = t / 2.03001708E05 # code unit to sec 
dt = dt / 2.03001708E05 # code unit to sec 
# data frame for GW in time domain, in code unit
df_TD = pd.DataFrame( {'t': t, 
            'I11_dot': I11_dot, 'I11_ddot': I11_ddot,
            'I22_dot': I22_dot, 'I22_ddot': I22_ddot,
            'I12_dot': I12_dot, 'I12_ddot': I12_ddot,
            'hp': h_plus,'hc': h_cross
             } )
df_TD = df_TD[['t','I11_dot','I11_ddot','I22_dot','I22_ddot','I12_dot','I12_ddot','hp','hc']]
df_TD.to_csv('GW_TD_xy.csv', index = False, sep = '\t')

# next is to do fft
window = signal.tukey(len(t), alpha = 0.01)
#window = np.hanning(len(t))
#window = np.blackman(len(t))
#window = np.kaiser(len(t), beta = 14)

# might be useful
#h_plus -= h_plus_max 
#h_plus -= np.mean(h_plus)
#h_plus = h_plus * window

fft_h_plus = fft(h_plus) * 2.0 / t.size
fft_h_cross = fft(h_cross) * 2.0 / t.size

freq = np.linspace(start=0.0, stop=1.0/(2.0 * dt), num=int(t.size/2) )
psd_h_plus = np.abs(fft_h_plus)**2

psd_h_cross = np.abs(fft_h_cross)**2

df_FD = pd.DataFrame( {'f': freq, 'fft_hp': np.abs(fft_h_plus[:t.size//2]), 'psd_hp': psd_h_plus[:t.size//2],
                                  'fft_hc': np.abs(fft_h_cross[:t.size//2]), 'psd_hc': psd_h_cross[:t.size//2]} )
df_FD = df_FD[['f','fft_hp','psd_hp','fft_hc','psd_hc']]
df_FD.to_csv('GW_FD_xy.csv', index = False, sep = '\t')

