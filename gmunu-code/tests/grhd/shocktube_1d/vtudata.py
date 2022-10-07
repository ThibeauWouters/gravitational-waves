import sys
sys.path.append('$HOME/gitlab/gmunu/tools/python')

from amrvac_pytools.vtkfiles import read, amrplot
import matplotlib.pyplot as plt
import pylab

d1=read.load_vtkfile(0,file='shocktube',type='vtu')

plt.plot(d1.getCenterPoints(),d1.rho/10, label='rho0/10', marker='.')

plt.legend(loc='best')
plt.grid(True)
plt.savefig('hydro.png')

