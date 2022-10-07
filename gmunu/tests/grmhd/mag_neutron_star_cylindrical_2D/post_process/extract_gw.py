import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import scipy.signal as signal
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':10})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)


filepath = "../"


# read data
df_TD = pd.read_csv(filepath+'output.log', index_col = False, delim_whitespace=True, usecols = ['global_time','I11_dot','I22_dot'])
df_TD = df_TD[['global_time','I11_dot','I22_dot']]
data_time = df_TD['global_time'].to_numpy()
data_I11_dot = df_TD['I11_dot'].to_numpy()
data_I22_dot = df_TD['I22_dot'].to_numpy()

I11_dot_interp = interp1d(data_time, data_I11_dot, kind = 'cubic')
I22_dot_interp = interp1d(data_time, data_I22_dot, kind = 'cubic')

t_min = data_time[0]
t_max = data_time[-1]
t = np.linspace(t_min,t_max,len(data_time))
dt = t[1]-t[0]

I11_dot = I11_dot_interp(t)
I22_dot = I22_dot_interp(t)
I11_ddot = np.gradient(I11_dot, dt)
I22_ddot = np.gradient(I22_dot, dt)
h_plus = I22_ddot - I11_ddot # fixme: this is not real hp

t = t / 2.03001708E05 # code unit to sec 
dt = dt / 2.03001708E05 # code unit to sec 
# data frame for GW in time domain, in code unit
df_TD = pd.DataFrame( {'t': t, 
            'I11_dot': I11_dot, 'I11_ddot': I11_ddot,
            'I22_dot': I22_dot, 'I22_ddot': I22_ddot,
            'hp': h_plus
             } )
df_TD = df_TD[['t','I11_dot','I11_ddot','I22_dot','I22_ddot','hp']]
df_TD.to_csv('GW_TD.csv', index = False, sep = '\t')

# next is to do fft
#window = signal.tukey(len(t), alpha = 0.5)
window = np.hanning(len(t))
#window = np.blackman(len(t))
#window = np.kaiser(len(t), beta = 14)

# might be useful
#h_plus -= h_plus_max 
h_plus -= np.mean(h_plus)
h_plus = h_plus * window

fft_h_plus = fft(h_plus) * 2.0 / t.size
freq = np.linspace(start=0.0, stop=1.0/(2.0 * dt), num=int(t.size/2) )
psd_h_plus = np.abs(fft_h_plus)**2

df_FD = pd.DataFrame( {'f': freq, 
               'fft': np.abs(fft_h_plus[:t.size//2]), 'psd': psd_h_plus[:t.size//2]} )
df_FD = df_FD[['f','fft','psd']]
df_FD.to_csv('GW_FD.csv', index = False, sep = '\t')

#plt.plot(t,I11_dot)
#plt.plot(t,I11_ddot-I22_ddot)
#plt.plot(freq, df_FD['fft'])

# Save the line plot
#plt.legend(loc='best')
#plt.xlabel('$R$')
#plt.ylabel('$D$')
#plt.xlim(0,8000)
#plt.ylim(-1e-2,1e-2)
#plt.grid(True)
#plt.show()
#plt.savefig('test.png')
