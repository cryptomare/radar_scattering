from utils import multi_roughness_spectrum
import math
import mpmath
from scipy.integrate import dblquad
import numpy as np
from math import sqrt

def multiscale_X(fr, sig1, L1, sig2, L2, sig3, L3, thi, er):

    er = er.conjugate()

    sig1 = sig1 * 100  # change from m to cm scale
    L1 = L1 * 100

    sig2 = sig2 * 100  # change from m to cm scale
    L2 = L2 * 100

    sig3 = sig3 * 100  # change from m to cm scale
    L3 = L3 * 100

    sig12 = sig1 ** 2
    sig22 = sig2 ** 2
    sig32 = sig3 ** 2

    sigs2 = sig12 + sig22 + sig32

    # - fr: frequency in GHz
    # - sig: rms height of surface in cm
    # - L: correlation length of surface in cm
    # - theta_d: incidence angle in degrees
    # - er: relative permittivity
    # - sp: type of surface correlation function

    error = 1.0e8

    k = 2 * math.pi * fr / 30 # wavenumber in free space.Speed of light is in cm / sec
    theta = thi * math.pi / 180 # transform to radian


    ks = k * math.sqrt(sigs2) # roughness parameter
    kl = k * L

    ks2 = ks * ks
    kl2 = kl ** 2

    cs = math.cos(theta)
    s = math.sin(theta + 0.001)

    s2 = s ** 2

    # -- calculation of reflection coefficints
    rt = sqrt(er - s2)


    rv = (er * cs - rt) / (er * cs + rt)


    rh = (cs - rt) / (cs + rt)

    rvh = (rv - rh) / 2

    # print(rt, rv, rh, rvh)
    # exit()


    # -- rms slope values
    sig_l = sig / L

    wn = multi_roughness_spectrum(sig12, sig22, sig32, sigs2, L1, L2, L3, k, s, Ts)

    rss = sqrt(2) / sigs2 * (sig12 * sig1 / L1 + sig22 * sig2 / L2 + sig32 * sig3 / L3)

    # --- Selecting number of spectral components of the surface roughness
    # if auto == 0:

    n_spec = 15 # numberofterms to include in the surface roughness spectra

    #
    #
    # if auto == 1:
    #     n_spec = 1
    #     while error > 1.0e-8:
    #         n_spec = n_spec + 1
    #         error = (ks2 * (2 * cs) ** 2) ** n_spec / math.factorial(n_spec)


    # -- calculating shadow consideration in single scat(Smith, 1967)

    ct = mpmath.cot(theta + 0.001)
    farg = (ct / sqrt(2) / rss).real


    gamma = 0.5 * (math.exp(-float(farg.real) ** 2) / 1.772 / float(farg.real) - erfc(float(farg.real)))

    Shdw = 1 / (1 + gamma)

    # -- calculating multiple scattering contribution
    # ------ a double integration function

    math.factorials = {}
    for number in np.arange(1,n_spec+1):
        math.factorials[number] = math.factorial(number)

    svh = dblquad(lambda phi, r : xpol_integralfunc_vec(r, phi, sp, xx, ks2, cs, s,
                                                            kl2, L, er, rss, rvh, n_spec, math.factorials),
                          0.1, 1, lambda x : 0, lambda x : math.pi)[0]




    svh = svh * 1.e-5 # # un - scale after rescalingin the integrand function.
    sigvh = 10 * math.log10(svh * Shdw)

    return(sigvh)