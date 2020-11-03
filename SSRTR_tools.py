import sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import ipywidgets as wg
import matplotlib.gridspec as gridspec
import math
from IBA_tools import k_IBA, calc_e_ice, wn_from_wl

sys.path.append('SSRTR')
sys.path.append('IIEM')
from SS_tools import dB_to_lin
from SSRTR.S2RTR_contributions import S2RTR_Contributions

def S2RTR_to_lin(r_dict):
    lin_results = defaultdict(lambda: defaultdict())
    
    for pol_key in r_dict.keys():
        contributions = r_dict[pol_key]
        for con_key in contributions.keys():
            lin_results[pol_key][con_key] = dB_to_lin(contributions[con_key])
    
    return(lin_results)

def make_interactive_SSRTR_plot():
    
    corr_len_slide = wg.FloatSlider(value=0.2,
                               min=0.1,
                               max=1.5,
                               step=0.1,
                               description=r'CL (mm)')
    
    snow_roughness_slide = wg.FloatSlider(value=0.5,
                           min=0.1,
                           max=3,
                           step=0.2,
                           description=r'$\sigma_{snow}$ (mm)')

    ice_roughness_slide = wg.FloatSlider(value=0.5,
                           min=0.1,
                           max=3,
                           step=0.2,
                           description=r'$\sigma_{ice}$ (mm)')

    
    depth_slide = wg.FloatSlider(value=0.1,
                               min=0.05,
                               max=1,
                               step=0.05,
                               description=r'Dep$_{snow}$ (m)')
    
    density_slide = wg.FloatSlider(value=300,
                               min=200,
                               max=750,
                               step=50,
                               description=r'$\rho_{s}$ (kgm$^{-3}$)')
    
    temp_slide = wg.FloatSlider(value=-1,
                               min=-20,
                               max=0,
                               step=1,
                               description=r'Temp ($^{\circ}$C)')
    
    incidence_angle_slide = wg.FloatSlider(value=30,
                               min=5,
                               max=80,
                               step=5,
                               description=r'Inc Angle ($^{\circ}$)')

    wg.interact(make_SSRTR_bar_plot,
                corr_len = corr_len_slide,
                snow_depth = depth_slide,
                snow_density = density_slide,
                snow_roughness = snow_roughness_slide,
                ice_roughness = ice_roughness_slide,
                temperature_deg = temp_slide,
                inc_angle_deg = incidence_angle_slide)
    
def e_dry_snow(e_ice, nu):
    
    if nu < 0.45:
        e_dash = 1 + 0.4667*nu + 1.435*(nu**3)
    else:
        e_dash = (1+0.4795*nu)**3
    
    e_double_dash_numerator = 9 * nu * e_ice.imag
    
    e_double_dash_denom_in_brackets = (2 + nu) + (e_ice.real*(1-nu))
    
    e_double_dash = e_double_dash_numerator/(e_double_dash_denom_in_brackets**2)
    
    return complex(e_dash, e_double_dash)

def make_SSRTR_bar_plot(corr_len,
                        snow_depth,
                        snow_density,
                        snow_roughness,
                        ice_roughness,
                        temperature_deg,
                        inc_angle_deg):
    
    e_air = 1
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    constituent_names = ['Ice Surf', 'Snow Vol BS', 'Snow Bi 1R', 'Snow BS 2R', 'Snow Surf']
    
    bar_centers = [1,1.4]
    level_centers = [1,1.2]
    bw = 0.1
    
    bar_xticks, lev_xticks = [], []
    for bar_center, lev_center in zip(bar_centers, level_centers):
        for offset in [-1,0,1]:
            bar_xticks.append(bar_center+(offset*bw/2))
            lev_xticks.append(lev_center+(offset*bw/2))
            
    xticklabels = ['HH','\n\nKu', 'VV', 'HH', '\n\nKa', 'VV']
    
    fig = plt.figure(figsize=(12,5))
    widths = [1,0.7]
    spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=widths)
    ax = fig.add_subplot(spec[0])
    ax2 = fig.add_subplot(spec[1])
    
    for i in [0,2,3,5]:
        ax2.axvline(lev_xticks[i], ls='--', c='k', alpha=0.5)
    ax2.set_ylim(-40,0)
    
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")

    ax2.set_xticks(lev_xticks)
    ax2.set_xticklabels(xticklabels, fontsize='x-large')    
    ax2.set_ylabel(r'Total Backscatter, $\sigma_0$ (dB)', fontsize = 'x-large')
    ax2.set_xlim(level_centers[0]-bw, level_centers[1]+bw)
    
    ax.set_xticks(bar_xticks)
    ax.set_xticklabels(xticklabels, fontsize='x-large')
    ax.set_ylim(0,1)
    ax.set_yticks(np.linspace(0,1,11))
    ax.set_yticklabels(np.arange(0,101,10))
    ax.set_ylabel('Backscatter Contribution (%)', fontsize = 'x-large')
    
    constituents = ['ice_surf', 'vol_scat', 'refl_bis', 'refl_bak_refl', 'snow_surf']
    
    for frequency, bar_position, level_position in zip([13.5,30],
                                                       [bar_centers[0],bar_centers[1]], 
                                                       [level_centers[0], level_centers[1]]):
        
        wavenumber = wn_from_wl(3e8/(frequency*1e9))
        
        e_ice = calc_e_ice(temperature_deg,
                                frequency)
        
        e_snow = e_dry_snow(e_ice, snow_density/917)

        IBA_results = k_IBA(nu=snow_density/917,
                  corr_len=corr_len/1000,
                  wavenumber = wavenumber,
                  e1=e_air,
                  e2=e_ice)

        albedo = IBA_results['k_s']/IBA_results['k_e']

        r_dict = S2RTR_Contributions(e_snow,
                                    e_ice,
                                    f=frequency,
                                    s1=snow_roughness/1000,
                                    s3=ice_roughness/1000,
                                    a=albedo,
                                    kappa_e=IBA_results['k_e'],
                                    d=snow_depth,
                                    theta=inc_angle_deg,
                                    ss_mode='IIEM',
                                    L=0.05)
        
        lin_r_dict = S2RTR_to_lin(r_dict)
        
        db_totals = {'vv':r_dict['vv']['total'],'hh':r_dict['hh']['total']}
        
        totals = {'vv':lin_r_dict['vv']['total'],'hh':lin_r_dict['hh']['total']}
        
        
        for pol, bar_offset in zip(['vv', 'hh'],[-bw/2,bw/2]):
            
            ax2.scatter(level_position+bar_offset, db_totals[pol], marker = '_', s=3800, c='k')
            
            cons = [0]
            for con in constituents:
                cons.append(lin_r_dict[pol][con]/totals[pol])

            cumsums = np.cumsum(cons)

            for i, name, color in zip(range(len(cons)-1), constituent_names, colors):
                bottom = cumsums[i]
                height = cumsums[i+1]-bottom

                if (bar_position > 1) or (pol == 'hh'): name = None
                bar = ax.bar(bar_position+bar_offset,height,bottom=bottom,
                        label=name, width = 0.1, color=color, edgecolor='black')

    ax.legend(loc='center')
    plt.show()
