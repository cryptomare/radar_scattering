# Code 11.2: S2RT / R Backscattering from Rayleigh Layer with Distinct Upper Boundary

# Description: Code computes sigma_0_vv and sigma_0_hh for a weakly
# scattering of a Rayleigh layer with albedo a < 0.2.The layer is above a ground
# surface, and its upper and lower interfaces are characterized by the PRISM
# model(code 10.5).

# Input Variables:
# eps2: complex dielectric constant of middle layer
# eps3: complex dielectric constant of ground surface
# f: frequency(GHz)
# s1: rms height of top surface(m)
# s3: rms height of ground surface(m)
# a: Single - scattering albedo of Rayleigh layer(unitless)
# kappa_e: extinction coefficient of the Rayleigh layer(Np / m)
# d: thickness of the Rayleigh layer(m)
# theta: Incidence angle(degrees)

# Output Variables:
# sigma_0_vv(dB)
# sigma_0_hh(dB)

# Book reference: Section 11 - 2.2 and eq 11.79

from PRISM1_ForwardModel import PRISM1_ForwardModel
import math
import cmath

import sys
sys.path.append('../IIEM')
from IIEM.I2EM import backscatter

def dB(x):
    return 10 * math.log10(x)

def S2RTR_Contributions(eps2, eps3, f, s1, s3, a, kappa_e, d, theta, ss_mode,L=0.05):

    theta_r = theta * math.pi / 180 # transform to radian

    thetapr = math.asin(math.sqrt(1 / eps2.real) * math.sin(theta_r)) # angle inside the layer  # 2
    thetapr = thetapr * 180 / math.pi # transform to degrees

    costhetapr = math.sqrt(1 - 1 / eps2.real * math.sin(theta_r) ** 2)

    kappa_s = a * kappa_e # scattering coefficient

    # ---- call the PRISM - 1 surface scattering model
    # --  scattering due to top interface


    if 'prism' in ss_mode.lower():
        sig_s_vv, sig_s_hh, sig_s_hv = PRISM1_ForwardModel(eps2, theta, s1, f)
    elif 'iiem' in ss_mode.lower():

        block_crosspol = False if 'cross' in ss_mode.lower() else True

        sig_s_vv, sig_s_hh, sig_s_hv = backscatter(f,
                                                   s1,
                                                   L,
                                                   theta,
                                                   eps2,
                                                   'math.exponential',
                                                   3.5,
                                                   block_crosspol=block_crosspol)

    sig_12_vv = 10 ** (sig_s_vv / 10) # transform to linear
    sig_12_hh = 10 ** (sig_s_hh / 10)
    sig_12_hv = 10 ** (sig_s_hv / 10)

    # --  scattering due to bottom interface
    # note here we have assumed that PRISM1 model can be still used by
    # considering the ratio between the dielectric constants across the surface
    # boundary as the new dielectric constant input(i.e.used the dielectric
    # contrast as input).

    if 'prism' in ss_mode.lower():
        sig_s_vv, sig_s_hh, sig_s_hv = PRISM1_ForwardModel(eps3 / eps2, thetapr, s3, f)
    elif 'iiem' in ss_mode.lower():
        sig_s_vv, sig_s_hh, sig_s_hv = backscatter(f,
                                                   s3,
                                                   L,
                                                   thetapr,
                                                   eps3/eps2,
                                                   'math.exponential',
                                                   3.5,
                                                   block_crosspol=block_crosspol)

    sig_23_vv = 10 ** (sig_s_vv / 10)
    sig_23_hh = 10 ** (sig_s_hh / 10)
    sig_23_hv = 10 ** (sig_s_hv / 10)

    # -- calculate transmissivity in layer
    tau = kappa_e * d / costhetapr
    T =  math.exp(-tau)

    # -- calculate reflectivity of upper interface
    t1, t2, t3, t4, t5, t6, Th_12, Tv_12 = ReflTransm_PlanarBoundary(1, eps2, theta)

    # -- calculate reflectivity of lower interface
    # t1 & t2 redefined here but neither of them get used

    t1, t2, gammah_23, gammav_23, t3, t4, t5, t6, = ReflTransm_PlanarBoundary(eps2, eps3, thetapr)

    # -- calculate the total backscattering coefficient according to eq 11.79


    # VV components

    vv_ice_surf = (Tv_12 ** 2) * (T ** 2) * sig_23_vv
    vv_vol_scat =  (Tv_12 ** 2) * 0.75 * a * costhetapr * (1 - T ** 2)
    vv_refl_bis =  (Tv_12 ** 2) * 0.75 * a * costhetapr * (1 - T ** 2) * (gammav_23 ** 2 * T ** 2)
    vv_refl_bak_refl = (Tv_12 ** 2) * 6 * kappa_s * d * gammav_23 * T ** 2
    vv_snow_surf =  sig_12_vv
    vv_total = vv_ice_surf + vv_vol_scat + vv_refl_bis + vv_refl_bak_refl + vv_snow_surf

    vv_db_dict = {'ice_surf': dB(vv_ice_surf),
                  'vol_scat': dB(vv_vol_scat),
                  'refl_bis': dB(vv_refl_bis),
                  'refl_bak_refl':dB(vv_refl_bak_refl),
                  'snow_surf':dB(vv_snow_surf),
                  'total':dB(vv_total),
                  }

    # HV compnents

    hv_cross_ice_surf = Tv_12 * Th_12 * (T ** 2 * sig_23_hv)
    hv_cross_snow_surf =  sig_12_hv
    hv_total = hv_cross_snow_surf + hv_cross_ice_surf

    hv_db_dict = {'ice_surf': dB(hv_cross_ice_surf),
                  'snow_surf':dB(hv_cross_snow_surf),
                  'total':dB(hv_total),
                  }

    # VV components

    hh_ice_surf = (Th_12 ** 2) * (T ** 2) * sig_23_hh
    hh_vol_scat =  (Th_12 ** 2) * 0.75 * a * costhetapr * (1 - T ** 2)
    hh_refl_bis =  (Th_12 ** 2) * 0.75 * a * costhetapr * (1 - T ** 2) * (gammav_23 ** 2 * T ** 2)
    hh_refl_bak_refl = (Th_12 ** 2) * 6 * kappa_s * d * gammav_23 * T ** 2
    hh_snow_surf =  sig_12_hh
    hh_total = hh_ice_surf + hh_vol_scat + hh_refl_bis + hh_refl_bak_refl + hh_snow_surf

    hh_db_dict = {'ice_surf': dB(hh_ice_surf),
                  'vol_scat': dB(hh_vol_scat),
                  'refl_bis': dB(hh_refl_bis),
                  'refl_bak_refl':dB(hh_refl_bak_refl),
                  'snow_surf':dB(hh_snow_surf),
                  'total':dB(hh_total),
                  }

    return({'hh':hh_db_dict,
            'hv':hv_db_dict,
            'vv':vv_db_dict})


def ReflTransm_PlanarBoundary(eps1, eps2, theta1d):

    theta1 = math.radians(theta1d)

    sin_theta2 = cmath.sqrt(eps1) / cmath.sqrt(eps2) * math.sin(theta1)
    cos_theta2 = cmath.sqrt(1 - sin_theta2 ** 2)

    rhoh = (cmath.sqrt(eps1) * math.cos(theta1) - cmath.sqrt(eps2) * cos_theta2) / (cmath.sqrt(eps1) * math.cos(theta1)
                                                                                  + cmath.sqrt(eps2) * cos_theta2)

    rhov = (cmath.sqrt(eps1) * cos_theta2 - cmath.sqrt(eps2) * math.cos(theta1)) / (cmath.sqrt(eps1) * cos_theta2
                                                                                  + cmath.sqrt(eps2) * math.cos(theta1))

    tauh = 1 + rhoh
    tauv = (1 + rhov) * (math.cos(theta1) / cos_theta2)

    gammah = abs(rhoh) ** 2
    gammav = abs(rhov) ** 2

    Th = 1 - gammah
    Tv = 1 - gammav

    return(rhoh, rhov, gammah, gammav, tauh, tauv, Th, Tv)

