# Import modules

import pygyre as pg

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


import numpy as np
import pandas as pd
import fortranformat as ff



import yaml


import os
from glob import glob

from scipy.optimize import curve_fit, minimize
from scipy.signal import savgol_filter, argrelextrema
from scipy import interpolate

from bokeh.palettes import Category20c


import json
import re
import vaex


cgrav = 6.67430e-8
msun = 1.3271244e26 / 6.67430e-8
rsun = 6.957e10

figdir = '/Users/eoin/Documents/Snapshot_Seismic/figs'
data_dir = '/Users/eoin/Documents/Snapshot_Seismic/data'







def plot_modified_N2(pf, N2_new, plot_opt=2, N2_old=None, xlims=None, ylims=None):

    # br_n = np.insert(N2_new[:-1], 0, 0)
    # r_n = np.insert(pf.radius_cm[::-1][:-1], 0, 0)
        
    buoy_r_new = get_buoy_radius(np.insert(N2_new[:-1], 0, 0), np.insert(pf.radius_cm[::-1][:-1], 0, 0))
    buoy_r = get_buoy_radius(np.insert(pf.brunt_N2[::-1][:-1], 0, 0), np.insert(pf.radius_cm[::-1][:-1], 0, 0))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 3))

    ax1.plot(buoy_r, 1000*np.sqrt(np.abs(pf.brunt_N2[::-1])), color='orange', linestyle='--', alpha=0.3)
    # ax1.plot(buoy_r, 1000*np.sqrt(np.abs(pf.brunt_N2_composition_term[::-1])), color='grey', alpha=0.4, label='comp')
    ax1.plot(buoy_r_new, 1000*np.sqrt(np.abs(N2_new)), color='green', zorder=9, label='Best profile')
    ax1.scatter(buoy_r_new, 1000*np.sqrt(np.abs(N2_new)), color='black', zorder=10, s=1)
    ax1.plot(buoy_r_new, 1000*np.sqrt(np.abs(N2_new - pf.brunt_N2_structure_term[::-1])),
             zorder=-10, color='blue', alpha=1, label='comp. term')
    # ax1.plot(buoy_r, 1000*np.sqrt(np.abs(N2_new)), color='green', label='new')
    # ax1.scatter(buoy_r, 1000*np.sqrt(np.abs(N2_new)), color='green', label='new')
       
    if N2_old is not None:
        ax1.plot(get_buoy_radius(N2_old, pf.radius_cm[::-1]), 1000*np.sqrt(np.abs(N2_old)),
                 color='grey', alpha=0.4, label='other', zorder=-2)
        # ax1.scatter(get_buoy_radius(N2_old, pf.radius_cm[::-1]), 1000*np.sqrt(np.abs(N2_old)),
        #          color='green', alpha=0.4, zorder=-2)        

    xlim_test = 1000*np.sqrt(np.abs(pf.brunt_N2_composition_term[::-1]))
    xlim1 = next((x[0] for x in enumerate(xlim_test) if x[1] > 0.05), 0)
    xlim2 = len(pf.dm) - next((x[0] for x in enumerate(xlim_test[::-1]) if x[1] > 0.05), 0) - 1

    if plot_opt == 1:
        ax2.plot(pf.radius_cm[::-1], 1000*np.sqrt(np.abs(pf.brunt_N2[::-1])), color='red', label='old')
        ax2.plot(pf.radius_cm[::-1], 1000*np.sqrt(np.abs(pf.brunt_N2_composition_term[::-1])), color='grey', alpha=0.4, label='old')
        ax2.plot(pf.radius_cm[::-1], 1000*np.sqrt(np.abs(N2_new)), color='green', label='new')
        
        ax2.set_xlim(pf.radius_cm[::-1][xlim1 - 20], pf.radius_cm[::-1][xlim2 + 20])

    elif plot_opt == 2:
        ax2.plot(1000*np.sqrt(np.abs(pf.brunt_N2[::-1])), color='red', label='old')
        ax2.plot(1000*np.sqrt(np.abs(pf.brunt_N2_composition_term[::-1])), color='grey', alpha=0.4, label='old')
        ax2.plot(1000*np.sqrt(np.abs(N2_new)), color='green', label='new')
        
        ax2.set_xlim(xlim1 - 20, xlim2 + 20)

    ax1.legend()
    if xlims is not None:
        ax1.set_xlim(xlims)
    else:
        ax1.set_xlim(0, 1)
        
    if ylims is not None:
        ax1.set_ylim(ylims)
    else:
        ax1.set_ylim(0, 5)
    
    return fig, ax1, ax2


def general_nprof(pf, br_vals, n_vals, rcurve_vals=None, scurve_vals=None):
    # makes general shape of brunt-vaisalla profile
    
    ncomp = pf.brunt_N2_composition_term[::-1]
    nstruct = pf.brunt_N2_structure_term[::-1]
    
    if rcurve_vals is None:
        rcurve_vals = [0] * len(br_vals)
    
    if scurve_vals is None:
        scurve_vals = [0.1] * len(br_vals)
    
    br = get_buoy_radius(pf.brunt_N2[::-1], pf.radius_cm[::-1])
    
    n_vals = [(0.001 * x) ** 2 for x in n_vals]
    
    for k in range(25):
    
        i_vals = [np.argmin(np.abs(br - x)) for x in br_vals]

        for j, x in enumerate(i_vals[1:], 1):
            if x - i_vals[j - 1] == 0:
                i_vals[j] = x + 1

        # Start of composition profile
        ncomp_start = np.full(i_vals[0] + 1, n_vals[0])
        ncomp_new = ncomp_start

        for ival_a, ival_b, nval_a, nval_b in zip(i_vals[:-1], i_vals[1:], n_vals[:-1], n_vals[1:]):
            # if ival_b > ival_a + 1:
            # rad_int = pf.radius_cm[::-1][ival_a:ival_b+1]
            rad_int = 10**pf.logP[::-1][ival_a:ival_b+1]
            # rad_int = 10**pf.logRho[::-1][ival_a:ival_b+1]
            rad_int = (rad_int - rad_int[0])/(rad_int[-1] - rad_int[0])
            n_int = nval_a - (nval_a - nval_b) * rad_int
            ncomp_new = np.append(ncomp_new, n_int[1:])


        # End of composition profile
        ncomp_end = np.full(len(pf.dm) - i_vals[-1] - 1, n_vals[-1])
        ncomp_new = np.append(ncomp_new, ncomp_end)
         

        N2_new = ncomp_new + nstruct
        # N2_new = ncomp + nstruct

    
        # add curves
        for rcurve, scurve, ival, nval in zip(rcurve_vals, scurve_vals, i_vals, n_vals):
            if rcurve > 0:
                N_new = 1000*np.sqrt(np.abs(N2_new))

                cum_dist_before = np.cumsum(np.sqrt(np.diff(br[:ival+1][::-1])**2 + np.diff(N_new[:ival+1][::-1])**2))
                cum_dist_after = np.cumsum(np.sqrt(np.diff(br[ival:])**2 + np.diff(N_new[ival:])**2))
                ia = ival - 1 - next((x[0] for x in enumerate(cum_dist_before) if x[1] > scurve), ival)
                ib = ival + 1 + next((x[0] for x in enumerate(cum_dist_after) if x[1] > scurve), len(pf.dm) - ival - 1)
                # ia = np.maximum(0, ival - 10)
                # ib = np.maximum(0, ival + 10)

                brsect = br[ia:ib+1]
                nsect = N2_new[ia:ib+1]

                brsect_norm = (brsect - brsect[0])/(brsect[-1] - brsect[0])

                x = np.linspace(0, 1, 10000)
                alpha = rcurve
                y = alpha/(x + alpha) - (alpha/(1 + alpha))*x


                theta1 = np.arctan((N_new[ib] - N_new[ival])/(br[ib] - br[ival]))
                theta2 = np.pi + np.arctan((N_new[ival] - N_new[ia])/(br[ival] - br[ia]))

                l2 = np.sqrt((N_new[ia] - N_new[ival])**2 + (br[ia] - br[ival])**2)
                l1 = np.sqrt((N_new[ib] - N_new[ival])**2 + (br[ib] - br[ival])**2)

                rot = np.array([[l1*np.cos(theta1), l2*np.cos(theta2)], [l1*np.sin(theta1), l2*np.sin(theta2)]])

                xnew, ynew = np.dot(rot, np.column_stack((x, y)).T)
                xnew = xnew + br[ival]
                ynew = ynew + N_new[ival]

                f_br_ynew = interpolate.interp1d(xnew, ynew, fill_value='extrapolate')

                n_new = f_br_ynew(brsect)

                n_new = (0.001 * n_new) ** 2

                N2_new[ia:ib+1] = n_new
            
        br_n = np.insert(N2_new[:-1], 0, 0)
        r_n = np.insert(pf.radius_cm[::-1][:-1], 0, 0)
        br = get_buoy_radius(br_n, r_n)
           
    return N2_new




def read_gyre_model(moddir):
    """
    Reads GYRE model and parses data into DataFrame
    
    Input:
        moddir (str) : directory of model
    
    Returns:
        df (pd.DataFrame) : Contains model data
    """
    
    columns = ['k', 'radius', 'mass', 'luminosity', 'P', 'T', 'rho', 'gradT', 'N2', 'gamma1', 'grad_ad',
              'v_T', 'opacity', 'dopacity_dt', 'dopacity_drho', 'eps_nuc', 'depsnuc_dt', 'depsnuc_drho',
              'omega']
    
    df = pd.read_csv(moddir, skiprows=1, header=None, delim_whitespace=True)
    
    df.columns = columns
    
    header = pd.read_csv(moddir, nrows=1, header=None, delim_whitespace=True)
    
    header_cols = ['Ngrid', 'Mass', 'R_photosphere', 'L_photosphere', 'version']
    
    header.columns = header_cols
    
    return df, header


def write_gyre_model(df, header, moddir_out):

    # Now save in MESA .data format
    ff.config.RECORD_SEPARATOR = ''
    iformat = ff.FortranRecordWriter('1I6')
    lineformat = ff.FortranRecordWriter('(1E27.16)')
    
    inputFilenew = open(moddir_out, "w")
    
    header_vals = header.iloc[0].values

    h1 = iformat.write([int(header_vals[0])])
    h2 = lineformat.write(header_vals[1:-1])
    h3 = ff.FortranRecordWriter('1I7').write([int(header_vals[4])])
    
    inputFilenew.write(h1 + h2 + h3 + '\n')
    
    # actually write each row to file 
    # this bit takes a few seconds if the file is long
    for index, row in df.iterrows():
        i = iformat.write([row[0]])
        n = lineformat.write(row[1:])
        inputFilenew.write(i + n + '\n')

    inputFilenew.close()


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]




def compute_buoy_radius(N2_vals, r_vals):
    # Compute buoyancy radius
    # Returns normalised buoyancy radius profile and absolute value of buoyancy radius
    n2 = N2_vals[1:]
    rad = r_vals

    istart = next((x[0] for x in enumerate(n2) if x[1] > 1e-10))
    n_div_r = np.sqrt(np.abs(n2[istart:]))/rad[1:][istart:]

    dr = np.diff(rad)[istart:]

    buoy_r = np.cumsum(n_div_r * dr)

    buoy_r = np.concatenate((np.full(istart + 1, 0), buoy_r))

    return buoy_r/buoy_r[-1], buoy_r[-1]





def plot_2d_arrows(figname=None):
    
    modnames = []
    p1_list = [0.365, 0.370, 0.375, 0.380, 0.385, 0.390, 0.395]
    p2_list = [0.395, 0.400, 0.405, 0.410, 0.415, 0.420, 0.425]
    
    p1_list = [0.395, 0.400, 0.405, 0.410, 0.415, 0.420, 0.425]
    p2_list = [0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045]
    
    
    for p1 in p1_list:
        p1_names = []
        for p2 in p2_list:
            modname = 'inv_s2_' + str(int(1000*p1)).zfill(3) + '_' + str(int(1000*p2)).zfill(3)
            p1_names.append(modname)
        modnames.append(p1_names)
        
    
    nmodes_desired = np.arange(2, 50)
    
    nmin = None
    nmax = None
    
    target_label = 'inv_s2_410_030'
    # target_summary = sig_process_period_spacing(target_label, 0, nmin=nmin, nmax=nmax)
    nmode_target, dp_target = get_period_spacing_model(target_label)

    n_min = 20
    n_gap = 18
    n_upto = n_min + n_gap + 1
    f = f2

    target_br_width = 1/get_variation_in_fit_params(nmode_target, dp_target, n_gap=n_gap, n_upto=n_upto, n_min=n_min, f=f)['P'][0]
    target_amp = get_variation_in_fit_params(nmode_target, dp_target, n_gap=n_gap, n_upto=n_upto, n_min=n_min, f=f)['A'][0]
    print(target_br_width, target_amp)

    diffsx_list = []
    diffsy_list = []
    
    for p1_names in modnames:
        br_widths = []
        amps = []
        for test_label in p1_names:
            nmode, dp = get_period_spacing_model(test_label)
            br_width = 1/(get_variation_in_fit_params(nmode, dp, n_gap=n_gap, n_upto=n_upto, n_min=n_min, f=f)['P'][0])
            br_widths.append(br_width)

            amp = get_variation_in_fit_params(nmode, dp, n_gap=n_gap, n_upto=n_upto, n_min=n_min, f=f)['A'][0]
            amps.append(amp)

        print(p1_names)
        print(np.array(br_widths))
        diffsx = [target_br_width - x for x in br_widths]

        diffsy = [x - target_amp for x in amps]
        
        
        print(np.round(diffsx, 3))

        diffsx_list.append(np.array(diffsx))
        diffsy_list.append(np.array(diffsy))
    
    fig, ax = plt.subplots(1, figsize=(2, 2))
    
    diffsx_list = np.array(diffsx_list).T
    diffsy_list = np.array(diffsy_list).T
    
    calx = np.abs(np.mean((diffsx_list[:, 0] - diffsx_list[:, -1])/(0.425 - 0.395)))
    caly = np.mean((diffsy_list[0] - diffsy_list[-1])/0.030)
    
    # calx = np.abs(np.mean((diffsx_list[0] - diffsx_list[-1])))
    # caly = np.mean((diffsy_list[:, 0] - diffsy_list[:, -1]))
    
    ax.quiver(p1_list, p2_list, 1.*diffsx_list/calx, 1*diffsy_list/caly, color='black')
    
    # ax.scatter(0.410, 0.030, marker='x', color='red')
    
    ax.set_xlabel(r'$\Pi_{\mu}$ - Position of Boundary')
    ax.set_ylabel('Width of Boundary in Norm. BR')
    
    fig.subplots_adjust(right=0.99, top=0.99, left=0.15, bottom=0.15)
    
    if figname is None:
        figname = 'arrows_'
    filename = figdir + '/' + figname
    if os.path.exists(filename + '.png'):
        j = 1
        while os.path.exists(filename + str(j) + '.png'):
            j = j + 1
        filename = filename + str(j)
    
    fig.savefig(filename)
    
    plt.close(fig)

    plt.show()
    



# def compute_grid():
#     # Template function to compute grid
#     moddir = '/Users/eoin/Documents/Snapshot_Seismic/models/evolutionary_test_grids'

#     subdir = 'diffov_diff_highres2'
#     gyre_file = moddir + '/' + subdir + '/LOGS_ms/profile14.data.GYRE'
#     profile_file = moddir + '/' + subdir + '/LOGS_ms/profile14.data'
    
#     df, header = read_gyre_model(gyre_file)
#     pf = MR(profile_file)


#     p1_list = [0.395, 0.400, 0.405, 0.410, 0.415, 0.420, 0.425]
#     p2_list = [0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045]

#     for p1, p2 in itertools.product(p1_list, p2_list):
#         modname = 'inv_s2_' + str(int(1000*p1)).zfill(3) + '_' + str(int(1000*p2)).zfill(3)
#         print("Doing", modname)
#         pf = MR(profile_file)
#         br_vals = [0.02, 0.052, p1 - p2, p1, 0.8]
#         n_vals = [0, 3.395, 2.3, 0.1866, 0]
#         rcurve_vals = [0.0, 0.0, 0.02, 0.002, 0]
#         scurve_vals = [0.1, 0.1, 0.1, 0.1, 0.1]
#         print(modname)
#         N2_new = general_nprof(pf, br_vals, n_vals, rcurve_vals, scurve_vals)

#         struct2 = {'N2': np.insert(N2_new[:-1], 0, 0)[::-1]}
#         compute_freq_given_struct(struct2, gyre_file, modname=modname, overwrite=1, template='template_294', log_flag=False)





def assign_params_to_mods(mods, **kwargs):
    for param in kwargs:
        for mod, val in zip(mods, kwargs[param]):
            setattr(mod, param, val)

    return mods



plt.style.use('paper_thin')


def plot_static_mods_overview(static_mods, figname='temp_static_mods', plot_labels=['A', 'B', 'C', 'D', 'E', 'F'],
                              xlims=None, ylims=None, obs=None, show_flag=False):
    '''
    plots N2 profile and period spacing patterns
    as a function of n_g and period
    '''

    fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(2, 3, figsize=(3.5433 * 2, 3.5433 * 4/3))
    
    colors = Category20c[20][::4]

    static_mods = assign_params_to_mods(static_mods, color=colors, label=plot_labels)

    if obs is not None:
        obs.set_nmode(static_mods[0].nmode, static_mods[0].p)
        static_mods.append(obs)

    for mod in static_mods:
        ax2.plot(mod.p, mod.dp/1e3, color=mod.color, label=mod.label, linewidth=1)
        ax2.scatter(mod.p, mod.dp/1e3, color=mod.color, marker='o', s=4, zorder=2)

        p_div_mean_dp = mod.p / np.mean(mod.dp_rot)
        dp_new = np.diff(p_div_mean_dp) * 86400
        ax5.plot(1e3*p_div_mean_dp[:-1], 1e3*dp_new/1e3, color=mod.color, label=mod.label, linewidth=1)
        ax5.scatter(1e3*p_div_mean_dp[:-1], 1e3*dp_new/1e3, color=mod.color, marker='o', s=4, zorder=2)
       
      
        ax3.plot(mod.nmode, mod.dp/1e3, color=mod.color, label=mod.label, linewidth=1)
        ax3.scatter(mod.nmode, mod.dp/1e3, color=mod.color, marker='o', s=4, zorder=2)
    
        ax4.plot(mod.nmode, mod.dp_rot/1e3, color=mod.color, label=mod.label, linewidth=1)
        ax4.scatter(mod.nmode, mod.dp_rot/1e3, color=mod.color, marker='o', s=4, zorder=2) 


        ax1.plot(mod.buoy_r, mod.N_plot, color=mod.color, label=mod.label, linewidth=1)
        ax0.plot(mod.mass, mod.h1, color=mod.color, label=mod.label, linewidth=1)

    xticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
    ax1.set_xticks(xticks) 
    ax1.set_xticklabels(xticks)

    xtick_minor = np.arange(0, 1, 0.05)
    ax1.set_xticks(xtick_minor, minor=True)

    ymax = np.minimum(ax1.get_ylim()[1], 6)

    yticks = np.arange(0, 5)
    ax1.set_yticks(yticks)
    ytick_minor = np.arange(0, ymax, 0.25)
    ax1.set_yticks(ytick_minor, minor=True)
    ax1.set_ylim(0, np.minimum(ax1.get_ylim()[1], 4))

    ax1.set_ylabel(r'N 10$^3$s')
    ax1.set_ylim(0, np.minimum(ax1.get_ylim()[1], 4))
    ax1.set_ylim(0, ymax)
 
    ax1.set_xlabel(r'$\mathcal{\Pi}_{0}/\mathcal{\Pi}_{\mathrm{r}}$')
    ax1.set_xlim(-0.0, 1.01)

    if xlims is not None:
        ax1.set_xlim(xlims)
    if ylims is not None:
        ax1.set_ylim(ylims)

    ax1.legend()


    ax2.set_xlabel(r'P (days)')
    ax2.set_ylabel(r'$\Delta$P / 1000s')

    ax2.yaxis.set_major_locator(MaxNLocator(4))
    ax2.yaxis.set_minor_locator(MaxNLocator(20))

    ax0.set_xlabel(r'Mass (M$_{\mathrm{\odot}}$)')
    ax0.set_ylabel(r'$^1$H')
    ax0.xaxis.set_major_locator(MaxNLocator(4))

    fig.subplots_adjust(right=0.99, top=0.99, left=0.1, bottom=0.12, wspace=0.35)

    custom_savefig(fig, figname)

    if show_flag:
        plt.show()
    else:
        plt.close(fig)

def custom_savefig(fig, figname):
    filename = figdir + '/' + figname
    if os.path.exists(filename + '.png'):
        j = 1
        while os.path.exists(filename + str(j) + '.png'):
            j = j + 1
        filename = filename + str(j)

    fig.savefig(filename + '.png', dpi=600)






def save_parameters_to_yaml(parameters, file_path):
    with open(file_path, 'w') as file:
        yaml.dump(parameters, file)

def load_parameters_from_yaml(file_path):
    with open(file_path, 'r') as file:
        parameters = yaml.safe_load(file)
    return parameters



