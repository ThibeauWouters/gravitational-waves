import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as signal
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d
from matplotlib import rc

filepath = "../"

simname = "output.log"

x, y = np.loadtxt(open(filepath+simname,'rt').readlines()[:-1], unpack=True ,dtype='float',usecols=(1,3),skiprows=1)
y = y / max(y)
x = x / 2.03001708E05 # code unit to sec 
df_TD = pd.DataFrame({'t' : x, 'rho' : y})
df_TD = df_TD[['t','rho']]
df_TD.to_csv("rho_c_TD.csv", index = False, sep = '\t')

t_min = 0.0
t_max = float(df_TD.nlargest(1,'t')['t'])
#t_max = 1.0E-2
print (t_max)
time_step = 1.0/2.0**12 * t_max
print (time_step)
time = np.arange(0.0, float(t_max), float(time_step))

window = signal.tukey(len(time), alpha = 0.5)
#window = np.hanning(len(time))
#window = np.blackman(len(time))
#window = np.kaiser(len(time), beta = 14)

print(np.array(df_TD['t']))
print(df_TD['rho'])
data_t = interp1d(df_TD['t'], df_TD['rho'], kind = 'cubic')
#data_t = interp1d(np.array(df_TD['t']), np.array(df_TD['rho']), kind = 'cubic')
data_temp = data_t(time)# * 1.0e5
data_max = df_TD['rho'][1]
print (data_max)
data_temp -= data_max 
data_temp -= np.mean(data_temp)
data_temp = data_temp * window

data_f = fft(data_temp) * 2.0/time.size
freq = np.linspace(start=0.0, stop=1.0/(2.0 * time_step), num=int(time.size/2) )
data_psd = np.abs(data_f)**2

df_FD = pd.DataFrame( {'f': freq, 'fft': np.abs(data_f[:time.size//2]), 'psd': data_psd[:time.size//2]} )
df_FD = df_FD[['f','fft','psd']]
df_FD.to_csv('rho_c_FD.csv', index = False, sep = '\t')
