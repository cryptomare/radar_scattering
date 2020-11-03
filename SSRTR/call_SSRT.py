## Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
## These MATLAB-based computer codes are made available to the remote
## sensing community with no restrictions Users may download them and
## use them as they see fit The codes are intended as educational tools
## with limited ranges of applicability, so no guarantees are attached to
## any of the codes.

import numpy as np
from S2RTR_DistinctUB import S2RTR_DistinctUB
import matplotlib.pyplot as plt

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
vv, hh, hv = [], [], []

for i in range(1,nt+1):
    sig_0_vv, sig_0_hh, sig_0_hv = S2RTR_DistinctUB(eps2,eps3,f,s1,s3,a,kappa_e,d, theta[i-1])
    vv.append(sig_0_vv)
    hh.append(sig_0_hh)
    # hv.append(sig_0_hv)

plt.plot(theta, vv)
plt.plot(theta, hh)
plt.plot(theta, hv)
# plt.ylim(-26,-6)
plt.show()

#
# figure(1)
#
# plot(theta, sig_0_vv, theta, sig_0_hh, theta, sig_0_hv)
# xlabel('\theta (deg)')
# ylabel('\sigma^0')
# legend('vv', 'hh', 'hv')
# grid
