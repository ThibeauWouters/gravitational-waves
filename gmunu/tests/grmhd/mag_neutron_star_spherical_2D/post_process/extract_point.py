import yt
import numpy as np
import pandas as pd
import glob

# extraction location
r_cut = 4.0
theta_cut = np.pi/4.0
phi_cut = 0.0

data_sets = glob.glob('../output*.dat')

time = np.array([])
rho_p = np.array([])
psi_p = np.array([])
alp_p = np.array([])
W_vel1_p = np.array([])
W_vel2_p = np.array([])
W_vel3_p = np.array([])
vel1_p = np.array([])
vel2_p = np.array([])
vel3_p = np.array([])
v1_p = np.array([])
v2_p = np.array([])
v3_p = np.array([])


for i, data_set in enumerate(data_sets):
  # Load the dataset
  ds = yt.load(data_set, geometry_override='spherical', unit_system='code')
  # get time
  # get 1D data set
  slice_data = ds.ortho_ray( 'r', (theta_cut, phi_cut) )
  # locate where is the data
  data_loc = min(range(len(slice_data['r'])), key=lambda i: abs(float(slice_data['r'][i])-r_cut))
  time = np.append(time,[ds.current_time])
  rho_p = np.append(rho_p, [slice_data['rho'][data_loc]])
  psi_p = np.append(psi_p, [slice_data['psi'][data_loc]])
  alp_p = np.append(alp_p, [slice_data['alp'][data_loc]])
  W_vel1_p = np.append(W_vel1_p, [slice_data['W_vel1'][data_loc]])
  W_vel2_p = np.append(W_vel2_p, [slice_data['W_vel2'][data_loc]])
  W_vel3_p = np.append(W_vel3_p, [slice_data['W_vel3'][data_loc]])
  # now workout the veloc1 to 3
  wv1 = float(slice_data['W_vel1'][data_loc])
  wv2 = float(slice_data['W_vel2'][data_loc])
  wv3 = float(slice_data['W_vel3'][data_loc])
  psi4 = float(slice_data['psi'][data_loc]**4)
  w = wv1**2 + r_cut**2 * wv2**2 + r_cut**2 * np.sin(theta_cut)**2 * wv3**2
  w = psi4 * w
  w = np.sqrt( 1.0 + w )
  vel1_p = np.append(vel1_p, [slice_data['W_vel1'][data_loc]/w])
  vel2_p = np.append(vel2_p, [slice_data['W_vel2'][data_loc]/w])
  vel3_p = np.append(vel3_p, [slice_data['W_vel3'][data_loc]/w])
  psi2 = float(slice_data['psi'][data_loc]**2)
  v1_p = np.append(v1_p, [ psi2 * wv1 / w ])
  v2_p = np.append(v2_p, [ psi2 * r_cut * wv2 / w ])
  v3_p = np.append(v3_p, [ psi2 * r_cut * np.sin(theta_cut) * wv3 / w ])


# Sort the values by 't'
srt = np.argsort(time)

df_TD = pd.DataFrame({'t' : time[srt], 
                      'rho' : rho_p[srt],
                      'psi' : psi_p[srt],
                      'alp' : alp_p[srt],
                      'W_vel1' : W_vel1_p[srt],
                      'W_vel2' : W_vel2_p[srt],
                      'W_vel3' : W_vel3_p[srt],
                      'vel1' : vel1_p[srt],
                      'vel2' : vel2_p[srt],
                      'vel3' : vel3_p[srt],
                      'v1' : v1_p[srt],
                      'v2' : v2_p[srt],
                      'v3' : v3_p[srt]   }  )

df_TD = df_TD[['t','rho','psi','alp','W_vel1','W_vel2','W_vel3','vel1','vel2','vel3','v1','v2','v3']]
df_TD.to_csv("point_data.csv", index = False, sep = '\t')

