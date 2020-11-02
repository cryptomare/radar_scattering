import math
from math import sin, cos, radians, asin
import numpy as np
import cmath
from matplotlib.widgets import Slider, Button, RadioButtons
import ipywidgets as wg
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
sys.path.append('IIEM')
from IIEM.I2EM import I2EM_Bistat_model, backscatter



def fresnel_R(angle, n1, n2, pol):
    
    # S polarization is H pol
    # p polarization is V pol
    
    angle_i = math.radians(angle)
    
    angle_t = np.arcsin((n1/n2)*math.sin(angle_i))
    
    if pol.lower() == 'h':
        numerator = (n1 * math.cos(angle_i)) - (n2 * math.cos(angle_t))
    elif pol.lower() == 'v':
        numerator = (n1 * math.cos(angle_t)) - (n2 * math.cos(angle_i))

    
    denominator = (n1 * math.cos(angle_i)) + (n2 * math.cos(angle_t))
    
    frac = numerator/denominator
    mod_frac = abs(frac)
    
    R = mod_frac**2
    
    return(R)

def fresnel_T(angle, n1, n2, pol):
    
    T = 1 - fresnel_R(angle, n1, n2, pol)
    
    return T

def make_fresnel_plot(n2):
    f_R_V = []
    f_R_H = []

    angles = range(0,91)

    for angle in angles:

        f_R_V.append(fresnel_R(angle,1,n2,'v'))
        f_R_H.append(fresnel_R(angle,1,n2,'h'))
        
    brewsters_angle = angles[np.argmin(np.array(f_R_V))]
    
    bas = f'{brewsters_angle}'

    plt.plot(angles,np.array(f_R_V)*100,label='V',c='crimson')
    plt.plot(angles,np.array(f_R_H)*100,label='H',c='darkblue')
    plt.scatter(brewsters_angle, np.min(f_R_V),marker='|',s=200, c='g', label= f"Brewster's Angle ({brewsters_angle}"+r"$^{\circ}$)")
    plt.ylabel('Reflected Power (%)', fontsize='x-large')
    plt.xlabel('Angle ($^{\circ}$)', fontsize='x-large')
    plt.ylim(0,100)
    plt.xlim(0,90)
    plt.legend()
    plt.show()
    
    
def make_interactive_fresnel_plot():
    
    
    n2_slider = wg.FloatSlider(value=2,
                               min=1.2,
                               max=4,
                               step=0.2,
                               description=r'Refr Index, n')
        
    wg.interact(make_fresnel_plot,
                    n2 = n2_slider)

def map_n_to_alpha(n):
    gradient = 1/6
    intercept = -1/6
    alpha = (gradient * n) + intercept
    return(alpha)

def make_snell_plot(inc_angle, n2, pol):

    inc_line_magnitude = 5
    n1 = 1
    intersection_point = [0,0]

    inc_line_start = [-inc_line_magnitude*sin(radians(inc_angle)),inc_line_magnitude*cos(radians(inc_angle))]
    inc_x_coords = [intersection_point[0], inc_line_start[0]]
    inc_y_coords = [intersection_point[1], inc_line_start[1]]

    ref_line_magnitude = 5*fresnel_R(inc_angle,n1,n2,pol)
    ref_line_start = [ref_line_magnitude*sin(radians(inc_angle)),ref_line_magnitude*cos(radians(inc_angle))]
    ref_x_coords = [intersection_point[0], ref_line_start[0]]
    ref_y_coords = [intersection_point[1], ref_line_start[1]]

    tra_line_magnitude = 5*fresnel_T(inc_angle,n1,n2,pol)
    tra_angle = asin((n1/n2)*sin(radians(inc_angle)))

    tra_line_start = [tra_line_magnitude*sin(tra_angle),-tra_line_magnitude*cos(tra_angle)]
    tra_x_coords = [intersection_point[0], tra_line_start[0]]
    tra_y_coords = [intersection_point[1], tra_line_start[1]]

    ref_prop = dict(arrowstyle=f"-|>,head_width={ref_line_magnitude/3},head_length={ref_line_magnitude/1.5}",
                shrinkA=0,shrinkB=0, color='green')

    tra_prop = dict(arrowstyle=f"-|>,head_width={tra_line_magnitude/3},head_length={tra_line_magnitude/1.5}",
                shrinkA=0,shrinkB=0, color='crimson')

    inc_prop = dict(arrowstyle=f"-|>,head_width={inc_line_magnitude/3},head_length={inc_line_magnitude/1.5}",
                shrinkA=0,shrinkB=0, color='darkblue')
    
    fig = plt.figure(figsize=(8,8))  # create a figure object
    ax = fig.add_subplot(1, 1, 1)

    ax.annotate("", xy=(tra_line_start[0], tra_line_start[1]),
                 xytext=(0,0), arrowprops=tra_prop, c = 'crimson')

    ax.annotate("", xy=(ref_line_start[0], ref_line_start[1]),
                 xytext=(0,0), arrowprops=ref_prop, c = 'green')

    ax.annotate("", xy=(inc_line_start[0]/2, inc_line_start[1]/2),
                 xytext=(inc_line_start[0], inc_line_start[1]), arrowprops=inc_prop, c = 'darkblue')
    
    ax.annotate(f"Reflected Power: {int(100*ref_line_magnitude/5)}%",
                 xy=(0.02,0.1), xycoords= 'axes fraction',fontsize='xx-large')
    ax.annotate(f"Transmit Power: {100-int(100*ref_line_magnitude/5)}%",
                 xy=(0.02,0.03), xycoords= 'axes fraction',fontsize='xx-large')
    
    ax.annotate(f"Air",
                 xy=(0.02,0.52), xycoords= 'axes fraction',fontsize='25')
    ax.annotate(f"Dielectric",
                 xy=(0.02,0.45), xycoords= 'axes fraction',fontsize='25')
    
    ax.set_aspect('equal')
    
    custom_lines = [Line2D([0], [0], color='darkblue', lw=2),
                    Line2D([0], [0], color='crimson', lw=2),
                    Line2D([0], [0], color='green', lw=2)]

    ax.set_ylim(-5,5)
    ax.set_xlim(-5,5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.plot(inc_x_coords, inc_y_coords,c='darkblue')
    ax.axhline(0,c='k')
    ax.axvline(0,c='k',ls='--')
    ax.legend(custom_lines, ['Incident', 'Transmitted', 'Reflected'], fontsize='xx-large', loc='upper right')
    
    alpha = map_n_to_alpha(n2)
    
    ax.axhspan(0,-5,alpha=alpha,color='gray')
    plt.show()

def make_interactive_snell_plot():

    inc_angle_slider = wg.FloatSlider(value=40,
                                   min=0,
                                   max=90,
                                   step=5,
                                   description=r'Inc. Angle ($^{\circ}$)')

    n2_slider = wg.FloatSlider(value=2,
                               min=1,
                               max=4,
                               step=0.2,
                               description=r'Refr Index, n')

    pol_button = wg.RadioButtons(options=['H', 'V'],
                                   description='Polarization')

    wg.interact(make_snell_plot,
                inc_angle = inc_angle_slider,
                n2 = n2_slider,
                pol = pol_button)


def dB_to_lin(decibels):
    lin = np.power(10,decibels/10)
    return(lin)

def make_IIEM_plot(inc_angle,
                   sigma_h,
                   frequency):
    
    angles = np.arange(-81, 81, 2)
    vv_list, hh_list = [], []

    r_min = -40
    r_max = 20

    for angle in angles:

        (vv, hh) = I2EM_Bistat_model(frequency=frequency,
                                      sigma_h=sigma_h/1000,
                                      CL=0.1,
                                      theta_i=inc_angle,
                                      theta_scat=angle,
                                      phi_scat=180,
                                      er=complex(3,0), # Dielectric constant
                                      sp='math.exponential',
                                      xx=1)
        
        vv_list.append(vv)
        hh_list.append(hh)

    # vv_list = dB_to_lin(np.array(vv_list))

    angles = np.radians(angles)



    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111, polar=True)

    ax.plot(angles, vv_list, color = 'darkblue', label = 'VV')
    ax.plot(angles, hh_list, color = 'crimson', label = 'HH')
    ax.plot([0,np.radians(inc_angle)], # theta values
            [r_min,r_max], color = 'green') # R values
    
    ref_prop = dict(arrowstyle=f"-|>,head_width=0.5,head_length=1",
                shrinkA=0,shrinkB=0, color='green')
    
    ax.annotate('',
            xy=(np.radians(inc_angle), (r_min+r_max)/2), #theta, radius
            xytext=(np.radians(inc_angle),20),  # theta, radius
            textcoords='data',
                        arrowprops=ref_prop,
               )
    

    top_tick = math.ceil(np.max(vv_list))
    bottom_tick = math.floor(np.min(vv_list))
    ax.set_rmax(r_max)
    ax.set_rmin(r_min)
    ax.set_theta_zero_location('N')
    ax.set_thetamin(90)
    ax.set_thetamax(-90)
    ax.set_rticks(np.arange(r_min,r_max+1,10))  # Less radial ticks
    ax.set_xlabel('Normalised Bistatic Radar Cross Section (dB)',labelpad=-40, fontsize='x-large')  # Move radial labels away from plotted line

    ax.legend(fontsize='x-large',ncol=2)

    plt.show()
    
def make_interactive_IIEM_plot1():

    inc_angle_slider = wg.FloatSlider(value=40,
                                   min=0,
                                   max=90,
                                   step=5,
                                   description=r'Inc. Angle ($^{\circ}$)')
    
    freq_slider = wg.FloatSlider(value=13,
                                   min=2,
                                   max=30,
                                   step=2,
                                   description=r'Freq (GHz)')
    
    sigma_slider = wg.FloatSlider(value=1,
                                   min=1,
                                   max=10,
                                   step=1,
                                   description=r'$\sigma_h$ (mm)')


    wg.interact(make_IIEM_plot,
                inc_angle = inc_angle_slider,
                sigma_h= sigma_slider,
                frequency= freq_slider,
#                 n2 = n2_slider,
               )
    

def make_IIEM_ka_ku_plot(inc_angle,
                   sigma_h,
                   n2,
                   CL,
                   pol):
    
    angles_deg = np.arange(-80, 81, 5)
    ku_hh_list, ka_hh_list, ku_vv_list, ka_vv_list = [], [], [], []

    r_min = -40
    r_max = 20

    for angle in angles_deg:
        
        if angle == 0:
            angle = 0.1
        
        for frequency in [13.5, 30]:

            (vv, hh) = I2EM_Bistat_model(frequency=frequency,
                                          sigma_h=sigma_h/1000,
                                          CL=CL/1000,
                                          theta_i=inc_angle,
                                          theta_scat=angle,
                                          phi_scat=180,
                                          er=complex(math.sqrt(n2),0), # Dielectric constant
                                          sp='math.exponential',
                                          xx=1)

            if frequency == 13.5:
                ku_hh_list.append(hh)
                ku_vv_list.append(vv)
            elif frequency == 30:
                ka_hh_list.append(hh)
                ka_vv_list.append(vv)

    angles = np.radians(angles_deg)
        
    index = np.where(angles_deg==inc_angle)[0][0]

    results = []
    for bistatic_list in [ku_hh_list, ku_vv_list, ka_hh_list, ka_vv_list]:
        results.append(bistatic_list[index])
        
        
    ###########################################################################################
    
    
    fig = plt.figure(figsize=(15,5))
    
    widths = [1.8,0.2]
    spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=widths)
    
    ##############################################################
    
    ax = fig.add_subplot(spec[0], polar=True)
    
    if pol == 'HH':
        ax.plot(angles, ku_hh_list, color = 'darkblue', label = f'Ku {pol}')
        ax.plot(angles, ka_hh_list, color = 'crimson', label = f'Ka {pol}')
    if pol == 'VV':
        ax.plot(angles, ku_vv_list, color = 'darkblue', label = f'Ku {pol}')
        ax.plot(angles, ka_vv_list, color = 'crimson', label = f'Ka {pol}')

    ax.plot([0,np.radians(inc_angle)], # theta values
            [r_min,r_max], color = 'green') # R values
    
    ref_prop = dict(arrowstyle=f"-|>,head_width=0.5,head_length=1",
                shrinkA=0,shrinkB=0, color='green')
    
    ax.annotate('',
            xy=(np.radians(inc_angle), (r_min+r_max)/2), #theta, radius
            xytext=(np.radians(inc_angle),20),  # theta, radius
            textcoords='data',
                        arrowprops=ref_prop,
               )
    

    ax.set_rmax(r_max)
    ax.set_rmin(r_min)
    ax.set_theta_zero_location('N')
    ax.set_thetamin(90)
    ax.set_thetamax(-90)
    ax.set_rticks(np.arange(r_min,r_max+1,10))  # Less radial ticks
    ax.set_xlabel('Normalised Bistatic Radar Cross Section (dB)',labelpad=-40, fontsize='x-large')  # Move radial labels away from plotted line

    ax.legend(fontsize='x-large',ncol=2)

    ##############################################################################

    # Get backscatter_values for incidence angle
   
    ax1 = fig.add_subplot(spec[1])
    
    bw = 0.1
    
    tick_positions = [1,1.35]
    x = np.array(tick_positions)
    
    ax1.bar(x-bw/2, [results[0],results[2]], width=bw, label = 'HH')
    ax1.bar(x+bw/2, [results[1],results[3]], width=bw, label = 'VV')
    
    ax1.set_ylabel('Backscatter (dB)', fontsize='xx-large')
    ax1.set_xlabel('Frequency', fontsize='xx-large')
    
    ax1.set_ylim(-40,0)
    
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(['Ku','Ka'],fontsize=20)
    
    ax1.legend(loc='lower right', fontsize='large', bbox_to_anchor=(-0.2, 0.85))
    
    plt.subplots_adjust(wspace=-0.4)
    plt.show()
    
def make_interactive_IIEM_plot2():

    inc_angle_slider = wg.FloatSlider(value=40,
                                   min=0,
                                   max=90,
                                   step=5,
                                   description=r'Inc. Angle ($^{\circ}$)')
    
    n2_slider = wg.FloatSlider(value=1.6,
                               min=1,
                               max=4,
                               step=0.2,
                               description=r'Refr Index, n')

    
    sigma_slider = wg.FloatSlider(value=1,
                                   min=1,
                                   max=10,
                                   step=1,
                                   description=r'$\sigma_h$ (mm)')
    
    CL_slider = wg.FloatSlider(value=10,
                                   min=3,
                                   max=20,
                                   step=1,
                                   description='CL (mm)')
    
    pol_button = wg.RadioButtons(options=['HH', 'VV'],
                                   description='Polarization')


    wg.interact(make_IIEM_ka_ku_plot,
                inc_angle = inc_angle_slider,
                sigma_h= sigma_slider,
                n2= n2_slider,
                CL=CL_slider,
                pol = pol_button,
               )
