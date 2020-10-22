import math
import cmath
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import ipywidgets as wg


Np_to_dB = 8.685889638

def calculate_e_eff(e1,e2,nu):
    
    first_bit = 2*e1 - e2 + 3*nu*(e2-e1)
    
    second_bit = (2*e1 - e2 + 3*nu*(e2-e1))**2 + (8*e1*e2)
    
    e_eff = (first_bit + cmath.sqrt(second_bit))/4
    
    return e_eff

def k_IBA(nu, corr_len, wavenumber, e1, e2):
    
    """
    
    nu: volume fraction of ice
    corr_len: Correlation Length
    wavenumber: 
    
    """
    
    e_eff = calculate_e_eff(e1,e2,nu)
    
    k_a = nu * wavenumber * e2.imag * abs((2*e_eff + e1)/(2*e_eff+e2))**2
    
    
    pre_mod = (2* corr_len**3 * wavenumber**4 * nu * (1-nu))/32
    
    mod_term = abs(((e2-e1)*(2*e_eff + e1))/(2*e_eff + e2))**2
    
    k_s = pre_mod * mod_term
    
    k_e = k_a + k_s
    

    return_dict = {'k_s' : k_s,
                   'k_a' : k_a,
                   'k_e' : k_e,
                   'CL':corr_len,
                   'wavenumber':wavenumber}
    
    return return_dict

def wn_from_wl(wavelength):

    return (2*math.pi/wavelength)

def calc_e_ice(T_C, f):
    
    """Takes temperature in C and returns real part of ice permittivity"""
    
    # Calculate e_prime
    
    e_prime = 3.1884 + (9.81e-3 * T_C)
    
    # Calculate e_double_prime
    
    T_K = T_C + 273.15
    
    theta = (300/T_K) - 1
    
    alpha = (0.00504+0.0062*theta) * math.exp(-22.1*theta)
    
    B1 = 0.0207
    B2 = 1.16e-11
    b = 335
    
    beta = (B1/T_K) * (math.exp(b/T_K) / ((math.exp(b/T_K) - 1)**2)) + B2*(f**2) + \
           math.exp(-9.963 + 0.0372*(T_K-273.16))
    
    
    e_double_prime = (alpha/f) + beta*f
    
    e = complex(e_prime, e_double_prime)
    
    return e

def make_plot(T, corr_len, nu, ylim, yscale):
    e1 = complex(1,0) # Air
    c = 3e8
    
    corr_len = corr_len*1e-3
    nu = nu/1023

    df_list = []

    for frequency in np.linspace(1,40,50):

        freq_r = frequency * 1e9
        wavelength = c/freq_r

        output_dict = k_IBA(nu,
                            corr_len = corr_len,
                            wavenumber= wn_from_wl(wavelength),
                            e1=e1,
                            e2=calc_e_ice(T, frequency))

        output_dict['frequency'] = freq_r
        output_dict['wavelength'] = wavelength

        df_list.append(output_dict)

    df = pd.DataFrame(df_list)
    
    plt.figure(figsize=(8,5))

    plt.plot(df['frequency']/1e9, df['k_a'] * Np_to_dB, label='Absorption')
    plt.plot(df['frequency']/1e9, df['k_s'] * Np_to_dB, label='Scattering')
    plt.plot(df['frequency']/1e9, df['k_e'] * Np_to_dB, label='Total Extinction', ls='--')

    if yscale.lower() == 'log':
        plt.yscale('log')
        plt.legend(loc='lower right', fontsize='x-large')
    elif yscale.lower() == 'linear':
        plt.ylim(0,ylim)
        plt.legend(loc='upper left', fontsize='x-large')
        
        
    plt.ylabel('Power loss (dB/m)',fontsize='x-large')
    plt.xlabel('Frequency (GHz)',fontsize='x-large')


    
    plt.show()

def make_interactive_plot():

    temp_slide = wg.FloatSlider(value=-4,
                                   min=-40,
                                   max=0,
                                   step=2,
                                   description=r'T ($^{\circ}$C)')

    CL_slide = wg.FloatSlider(value=0.1,
                               min=0.1,
                               max=1,
                               step=0.05,
                               description=r'CL (mm)')

    dens_slide = wg.FloatSlider(value=450,
                                min=150,max=600,
                                step=25,
                                description=r'$\rho$, (kgm$^{-3}$)')


    ylim_slide = wg.FloatSlider(value=10,
                               min=1,
                               max=30,
                               step=3,
                               description=r'y axis lim')

    yscale_button = wg.RadioButtons(options=['Linear', 'Log'],
                                   description='Y Scale:')

    wg.interact(make_plot,
                T = temp_slide,
                corr_len = CL_slide,
                nu = dens_slide,
                ylim = ylim_slide,
                yscale = yscale_button)
