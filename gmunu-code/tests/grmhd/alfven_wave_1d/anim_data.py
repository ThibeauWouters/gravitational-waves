import yt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(-1.2, 1.2))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    ds = yt.load('output'+str(i).zfill(4)+'.dat')
    # cutting the x-axis through the y=0,z=0 
    plotdata = ds.ortho_ray( 0, (0, 0) )
    srt = np.argsort(plotdata['x'])
    x = np.array(plotdata['x'][srt])
    y = np.array(plotdata['Bvec3'][srt])
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=21, interval=1, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('basic_animation.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

plt.show()
