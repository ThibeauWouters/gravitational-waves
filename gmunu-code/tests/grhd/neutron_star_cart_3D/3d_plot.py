import matplotlib
matplotlib.use('Agg')
import pylab # to show the plot

import yt
import numpy as np

savepath = '/users/chi-kit.cheong/public_html/'

units_override = dict( 
length_unit=(1.0, 'cm'),
time_unit=(1.0, 's'),
mass_unit=(1.0, 'g'),
)

# Load the dataset
ds = yt.load('output0028.dat', unit_system='code')

sc = yt.create_scene(ds, lens_type = 'perspective')
# Get a reference to the VolumeSource associated with this scene
# It is the first source associated with the scene, so we can refer to it
# using index 0.
source = sc[0]

# Set the bounds of the transfer function
#source.tfh.set_bounds((3e-31, 5e-27))

# set that the transfer function should be evaluated in log space
source.tfh.set_log(True)

# Make underdense regions appear opaque
source.tfh.grey_opacity = False

# Plot the transfer function, along with the CDF of the density field to
# see how the transfer function corresponds to structure in the CDF
source.tfh.plot(savepath + 'transfer_function.png', profile_field='rho')

# save the image, flooring especially bright pixels for better contrast
#sc.save(savepath + 'rendering.png')
sc.save(savepath + 'rendering.png', sigma_clip=6.0)
