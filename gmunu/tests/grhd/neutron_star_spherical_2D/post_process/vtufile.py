import pylab
from gmunu_pytools.vtkfiles import read, amrplot

ds = read.load_vtkfile(100,file='output',type='vtu')

p1 = amrplot.polyplot(ds.rho, ds)
#p1.show(ds.rho, ds)
