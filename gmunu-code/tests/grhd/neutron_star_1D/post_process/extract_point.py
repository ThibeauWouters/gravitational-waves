import yt
import numpy as np
import pandas as pd
import glob

# extraction location
r_cut = 5.0
theta_cut = 0.0
phi_cut = 0.0

data_sets = glob.glob('../output*.dat')

time = np.array([])
rho_p = np.array([])
psi_p = np.array([])
alp_p = np.array([])
W_vel1_p = np.array([])
vel1_p = np.array([])


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
  # now workout the veloc1 to 3
  wv1 = float(slice_data['W_vel1'][data_loc])
  psi4 = float(slice_data['psi'][data_loc]**4)
  w = wv1**2
  w = psi4 * w
  w = np.sqrt( 1.0 + w )
  vel1_p = np.append(vel1_p, [slice_data['W_vel1'][data_loc]/w])


# Sort the values by 't'
srt = np.argsort(time)

df_TD = pd.DataFrame({'t' : time[srt], 
                      'rho' : rho_p[srt],
                      'psi' : psi_p[srt],
                      'alp' : alp_p[srt],
                      'W_vel1' : W_vel1_p[srt],
                      'vel1' : vel1_p[srt]  }  )

df_TD = df_TD[['t','rho','psi','alp','W_vel1','vel1']]
df_TD.to_csv("point_data.csv", index = False, sep = '\t')

