import gmunu_pytools as apt
import matplotlib.pyplot as plt

ds = apt.load_datfile('output0100.dat')
ds.get_info()
data = ds.load_all_data()

ds.get_coordinate_arrays()

x = ds.get_coordinate_arrays()
rho_data = data['rho']
#r_max, theta_max = ds.domain_width


#p = ds.amrplot('rho', polar=True)

p = ds.amrplot('rho', draw_mesh=True, mesh_linewidth=0.5, mesh_color='white',
#               mesh_linestyle='solid', mesh_opacity=0.8)
#p.ax.set_title("Density plot")
#p.fig.tight_layout()

#ds.show()

ax = plt.subplot(111, polar=True)
ctf = ax.contourf(x[1],x[0],p)
plt.colorbar(ctf)
plt.show()
