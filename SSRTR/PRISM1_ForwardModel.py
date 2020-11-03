# Code 10.5: PRISM - 1 - Forward Model

# Description: Code computes sigma_0 for all three polarization
# combinations, given the surface parameters.

# Input Variables:
# eps = eps' - j eps'': Complex dielectric constant of the scattering
# medium
# theta: Incidence angle(degrees)
# s: rms height(m)
# f: Frequency(GHz)

# Output Products:
# sigma_0_vv(dB)
# sigma_0_hh(dB)
# sigma_0_hv(dB)

# Book Reference: Section 10 - 5

import math, cmath

def PRISM1_ForwardModel(eps, theta, s, f):

    ks = s * (2 * math.pi * f / 0.3) # calculate roughness parameter

    theta = theta * (math.pi / 180.0) # incidence angle in radian

    gamma0 = Fresn_Refl0(eps) # reflectivity(normal incidence)

    [gammav, gammah] = Fresn_Refl(eps, theta)

    p = (1 - (2 * theta / math.pi) ** (1 / (3 * gamma0)) *  math.exp(-ks)) ** 2
    q = 0.23 * math.sqrt(gamma0) * (1 -  math.exp(-ks))

    g = 0.70 * (1 -  math.exp(-0.65 * ks ** 1.8))

    sigvv = g * (math.cos(theta)) ** 3 / math.sqrt(p) * (gammav + gammah)

    sig_0_vv = 10 * math.log10(sigvv)
    sig_0_hh = 10 * math.log10(sigvv * p)
    sig_0_hv = 10 * math.log10(sigvv * q)

    return (sig_0_vv, sig_0_hh, sig_0_hv)


def Fresn_Refl0(eps):

    # calculates Fresnel reflectivity at normal incidence.
    gamma0 = (abs((1 - cmath.sqrt(eps)) / (1 + cmath.sqrt(eps)))) ** 2
    return(gamma0)

# -----------------------------------------------------------------
def Fresn_Refl(eps, theta):

    # calculates Fresnel reflectivities of v and h - polarizations at given set
    # of incidence angles.

    [rho_v, rho_h] = refl_coef(theta, 1, eps)
    gammav = (abs(rho_v)) ** 2
    gammah = (abs(rho_h)) ** 2

    return(gammav, gammah)

# ----------------------------------------------------------------
def refl_coef(the1, eps1, eps2):

    # calculates the v and h - polarized reflection coefficients of a plane
    # dielectric surface

    n1 = cmath.sqrt(eps1)
    n2 = cmath.sqrt(eps2)
    costh2 = cmath.sqrt(1 - (n1 * math.sin(the1) / n2) ** 2)

    rho_v = -(n2 * math.cos(the1) - n1 * costh2) / (n2 * math.cos(the1) + n1 * costh2)
    rho_h = (n1 * math.cos(the1) - n2 * costh2) / (n1 * math.cos(the1) + n2 * costh2)
    return(rho_v, rho_h)