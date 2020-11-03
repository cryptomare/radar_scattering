## Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
## These MATLAB-based computer codes are made available to the remote
## sensing community with no restrictions Users may download them and
## use them as they see fit The codes are intended as educational tools
## with limited ranges of applicability, so no guarantees are attached to
## any of the codes.

import numpy as np
from S2RTR_contributions import S2RTR_Contributions
import matplotlib.pyplot as plt
import pandas as pd

f = 5.5 #frequency (GHz)
s1 = 0.001 # rms surface height of top interface (m)
s3 = 0.01 # rms surface height of bottom interface (m)

eps2 = 1.5 # effective dielectric of layer
eps3 = complex(6.28,-1.53) # corresponding to mv = 0.16

a = 0.15 #albedo must be < 0.2

kappa_e = 0.15  # extinction coefficient in Np/m

d = 1 #layer thickness in meters


theta = np.arange(5,60,5) #incidence angle (deg)
nt = len(theta)

#--calculate sigma_0 for the Rayleigh layer
output_list = []

for i in range(1,nt+1):
    output = S2RTR_Contributions(eps2,eps3,f,s1,s3,a,kappa_e,d, theta[i-1], ss_mode = 'iiem')
    output_list.append(output)

col_names = [
    "ice_surf",    "vol_scat",    "refl_bis",    "refl_bak_refl",
    "snow_surf",    "cross_ice_surf",    "cross_snow_surf",

    ]
df = pd.DataFrame(output_list,columns=col_names)

df.index = theta

print(df.head())

for col_name in col_names:
    if 'cross' in col_name:
        ls = '--'
    else:
        ls = '-'
    plt.plot(theta,df[col_name],label=col_name,ls=ls)
plt.legend()
plt.show()

#
# figure(1)
#
# plot(theta, sig_0_vv, theta, sig_0_hh, theta, sig_0_hv)
# xlabel('\theta (deg)')
# ylabel('\sigma^0')
# legend('vv', 'hh', 'hv')
# grid
