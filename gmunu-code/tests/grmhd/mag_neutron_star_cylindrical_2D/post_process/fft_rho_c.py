import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as signal
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d
from matplotlib import rc

filepath = "../"

simname = "output.log"

df_TD = pd.read_csv("point_data.csv", delimiter = '\t', usecols = ['t','rho','psi','alp','vel1','vel2','vel3'])
df_TD = df_TD[['t','rho','psi','alp','vel1','vel2','vel3']]

#y = y / max(y)
df_TD['t'] = df_TD['t'] / 2.03001708E05 # code unit to sec 

t_min = 0.0
#t_min = 1.0E-4
t_max = float(df_TD.nlargest(1,'t')['t'])
#t_max = 1.0E-2
print (t_max)
time_step = 1.0/2.0**12 * (t_max - t_min)
print (time_step)
time = np.arange(0.0, float(t_max), float(time_step))

window = signal.tukey(len(time), alpha = 0.5)
#window = np.hanning(len(time))
#window = np.blackman(len(time))
#window = np.kaiser(len(time), beta = 14)

#print(np.array(df_TD['t']))
#print(df_TD['rho'])

data_t1 = interp1d(df_TD['t'], df_TD['vel1'], kind = 'cubic')
data_t2 = interp1d(df_TD['t'], df_TD['vel2'], kind = 'cubic')
data_temp1 = data_t1(time)# * 1.0e5
data_temp2 = data_t2(time)# * 1.0e5

#data_max1 = np.max(data_temp1)
#data_max2 = np.max(data_temp2)
#data_temp1 -= data_max1
#data_temp2 -= data_max2
data_temp1 -= np.mean(data_temp1)
data_temp2 -= np.mean(data_temp2)

data_temp1 = data_temp1 * window
data_temp2 = data_temp2 * window

data_f1 = fft(data_temp1) * 2.0/time.size
data_f2 = fft(data_temp2) * 2.0/time.size
freq = np.linspace(start=0.0, stop=1.0/(2.0 * time_step), num=int(time.size/2) )
data_psd1 = np.abs(data_f1)**2
data_psd2 = np.abs(data_f2)**2

df_FD = pd.DataFrame( {'f': freq, 'fft1': np.abs(data_f1[:time.size//2]), 'psd1': data_psd1[:time.size//2], 'fft2': np.abs(data_f2[:time.size//2]), 'psd2': data_psd2[:time.size//2]} )
df_FD = df_FD[['f','fft1','psd1','fft2','psd2']]
df_FD.to_csv('W_vel_FD.csv', index = False, sep = '\t')
