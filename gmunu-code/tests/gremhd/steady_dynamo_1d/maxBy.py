import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab # to show the plot

from matplotlib import rc_context

import yt
import numpy as np

eta = 1.0E-2

x, y = np.loadtxt("./output.log", unpack=True ,dtype='float',usecols=(1,4),skiprows=1)

plt.plot(x, y)
plt.yscale('log')
plt.ylabel('$\\max(B^y)$')
plt.xlabel('$t$')
plt.xlim(0,100)
plt.ylim(1E-2,1E5)



#plt.legend(loc='best')
with rc_context({'mathtext.fontset': 'stix'}):
  plt.savefig('Bymax.png')
