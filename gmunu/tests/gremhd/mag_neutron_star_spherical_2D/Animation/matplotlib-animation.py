import yt
#import matplotlib
#matplotlib.use('Agg')
from matplotlib.animation import FuncAnimation
from matplotlib import rc_context
#import pylab # to show the plot

ts = yt.load('../output*.dat', geometry_override='spherical', unit_system='code')

plot = yt.SlicePlot(ts[0], 2, 'rho')
plot.set_cmap(field="rho", cmap='hot')
plot.set_zlim('rho', 1.4E-3, 1e-10)

Radius = ts[0].domain_width[0]
plot.set_center((Radius/2.0,Radius/2.0))
plot.set_width(Radius,Radius)
plot.annotate_timestamp(corner='upper_right',draw_inset_box=True)

fig = plot.plots['rho'].figure


# animate must accept an integer frame number. We use the frame number
# to identify which dataset in the time series we want to load
def animate(i):
    ds = ts[i]
    plot._switch_ds(ds)
    #plot.save('hydro'+str(i)+'.png')

animation = FuncAnimation(fig, animate, frames=len(ts))

# Override matplotlib's defaults to get a nicer looking font
with rc_context({'mathtext.fontset': 'stix'}):
    animation.save('animation.mp4',fps=20)
