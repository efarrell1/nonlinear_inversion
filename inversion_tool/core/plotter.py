# Import modules

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bokeh.palettes import Category20c, Greens, Blues

from .stellar_models import GyreMod, get_obs, Observation, PeriodSpacing
from .mesa_reader import MesaData as MR
from .run_build_star import find_nearest_model

import os



def plot_each_timestep(inverter):
    # log(inverter.log_file, "Saving figure")
    fig, (ax) = plt.subplots(1)

    mod = GyreMod(inverter.gyre_subdir, m=inverter.m)

    ps = mod.set_dp_for_nmin_nmax(inverter.inputs['nmode_fit'], inverter.inputs['nmode_fit'] + len(inverter.obs.dp) - 1)

    ax.plot(inverter.obs.p, inverter.obs.dp, marker='o', color='grey', linewidth=2)
    ax.plot(inverter.obs.p[inverter.nmin:inverter.nmax], inverter.obs.dp[inverter.nmin:inverter.nmax], color='k', linewidth=1.2)

    ax.plot(ps.p, ps.dp, color=Greens[5][0], linewidth=1.2)
    ax.scatter(ps.p, ps.dp, color=Greens[5][0], marker='o', s=5, linewidth=1.2)

    if inverter.n_iter > 1:
        mod = GyreMod(inverter.previous_iterations_gyre_subdir[-1], m=inverter.m)
        ps = mod.set_dp_for_nmin_nmax(inverter.inputs['nmode_fit'], inverter.inputs['nmode_fit'] + len(inverter.obs.dp) - 1)
        ax.plot(ps.p, ps.dp, color=Greens[5][2], marker='o', linewidth=1., zorder=-1)


    xtext = 0.98
    ytext = 0.98

    for key, value in inverter.profile_parameters.items():
        ax.text(xtext, ytext, key + ' ' + "{:.5f}".format(value), transform=ax.transAxes, ha='right', va='top', fontsize=5, color='k')
        ytext = ytext - 0.04

    ax.set_xlabel('P (d)')
    ax.set_ylabel(r'$\Delta$ P')

    fig.subplots_adjust(right=0.99, top=0.99, left=0.15, bottom=0.15)

    custom_savefig(fig, inverter.fig_dir, "period_spacing_" + str(inverter.n_iter).zfill(5))

    plt.close(fig)


def plot_profiles(inverter):
    # log(inverter.log_file, "Saving figure")
    fig, (ax) = plt.subplots(1)

    pf = MR(inverter.static_subdir + '/static.data')

    ax.plot(pf.mass/pf.mass[0], pf.h1, color='green', linewidth=1.5)

    database_dir =  "/Users/eoin/Documents/Snapshot_Seismic/models/snapshot_database"
    lookup_table_dir = database_dir + '/lookup_table'

    nearest_model = find_nearest_model(pf.mass[0],
    							   	   inverter.profile_parameters['zval_lookup'],
    							   	   inverter.profile_parameters['xcore'],
    							   	   lookup_table_dir)
    mod = MR(nearest_model, file_type='model')


    starmass = getattr(mod, 'M/Msun')
    # mass = starmass * (1 - np.cumsum(oldmodel.dq))
    mass = starmass * (np.cumsum(mod.dq[::-1]))[::-1]


    ax.plot(mass/mass[0], mod.h1, color='grey', linewidth=1.5)

    ax.set_xlabel(r'M$_{\mathrm{r}}$')
    ax.set_ylabel(r'H1')

    ax.set_xlim(0, 1.01)
    ax.set_ylim(0, 0.75)

    fig.subplots_adjust(right=0.99, top=0.99, left=0.15, bottom=0.15)

    custom_savefig(fig, inverter.fig_dir, "profiles_" + str(inverter.n_iter).zfill(5))

    plt.close(fig)


def plot_main_iteration_ps(inverter):

    # log(inverter.log_file, "Saving figure")
    fig, (ax) = plt.subplots(1)

    ax.plot(inverter.obs.p, inverter.obs.dp, color='grey', linewidth=2, zorder=-2)
    ax.scatter(inverter.obs.p, inverter.obs.dp, marker='o', s=6, color='grey', linewidth=2, zorder=-2)
    ax.plot(inverter.obs.p[inverter.nmin:inverter.nmax], inverter.obs.dp[inverter.nmin:inverter.nmax], color='k', linewidth=0.4, zorder=3)

    mods = [GyreMod(gyre_subdir, m=inverter.m) for gyre_subdir in inverter.previous_iterations_gyre_subdir]

    colors = Greens[9][:len(mods)+1][1:][::-1]
    colors = Blues[9][:len(mods)+1][::-1]
    labels = np.arange(1, len(mods) + 1)

    for mod, color, label in zip(mods, colors, labels):
        if label == len(mods):
            alpha = 1
            linewidth = 1.2
            linestyle = '-'
            zorder = 10
            s = 5
        else:
            alpha = 1
            linewidth = 0.6
            linestyle = ':'
            zorder = -1
            s = 2
        if label == 1:
            plot_label = 'starting guess'
        else:
            plot_label = "Iter. " + str(label - 1)
        ps = mod.set_dp_for_nmin_nmax(inverter.inputs['nmode_fit'], inverter.inputs['nmode_fit'] + len(inverter.obs.dp) - 1)
        ax.plot(ps.p, ps.dp, color=color, linewidth=linewidth, label=plot_label, zorder=zorder, linestyle=linestyle)
        # if label == len(mods) or len(mods) == 1:
        ax.scatter(ps.p, ps.dp, color=color, marker='o', alpha=0.8, s=s, zorder=zorder)

    xtext = 0.03
    ytext = 0.02

    allowed_vals = inverter.inputs['jacobian_profile_params']

    for key, value in inverter.profile_parameters.items():
        if key == 'mass':
            key = 'mass'
        if key in allowed_vals:
            ax.text(xtext, ytext, key + ': ' + "{:.5f}".format(value), transform=ax.transAxes, ha='left', va='bottom', fontsize=5, color='k')
            ytext = ytext + 0.04

    ax.set_xlabel('P (d)')
    ax.set_ylabel(r'$\Delta$ P')

    legend = ax.legend(loc='upper right', fontsize=6, facecolor='white', frameon=True, edgecolor='white')
    legend.set_zorder(100)


    fig.subplots_adjust(right=0.99, top=0.99, left=0.15, bottom=0.15)

    custom_savefig(fig, inverter.fig_dir, "jacobian_iteration_ps_" + str(inverter.n_jacob_iter).zfill(5))

    plt.close(fig)





def plot_fit_params(ps, fig_dir, n_iter, nmin, nmax):
    df = ps.get_fit_params(nmin, nmax)

    # fit_params = (df['P'].values[0], df['A'].values[0], df['Phi'].values[0], df['c'].values[0])
    fit_params = (df['P'].values[0], df['A'].values[0], df['Phi'].values[0], df['c'].values[0])

    fit_vals = ps.fit_fun(ps.nmode, *fit_params)
    nfit_continuous = np.linspace(ps.nmode[0], ps.nmode[-1], 1000)

    fit_vals_cont = ps.fit_fun(nfit_continuous, *fit_params)

    fig, (ax) = plt.subplots(1)

    ax.plot(ps.nmode, ps.dp_rot, marker='o', color='black', linewidth=1.2)
    ax.plot(ps.nmode, fit_vals, color='green')
    ax.plot(nfit_continuous, fit_vals_cont, color='blue')

    ax.set_xlabel('P (d)')
    ax.set_ylabel(r'$\Delta$ P')

    fig.subplots_adjust(right=0.99, top=0.99, left=0.15, bottom=0.15)

    custom_savefig(fig, fig_dir, "period_spacing_fit_" + str(n_iter).zfill(5))

    plt.close(fig)



def plot_slope(ps, fig_dir, n_iter, nmin, nmax):
    
    popt, pcov = curve_fit(flinear, ps.p, np.array(ps.dp))

    fig, (ax) = plt.subplots(1)

    xplot = np.linspace(ps.p[0], ps.p[-1], 1000)
    yplot = popt[1] + xplot * popt[0]

    ax.plot(ps.p, ps.dp, marker='o', color='black', linewidth=1.2)
    ax.plot(xplot, yplot, color='green')

    ax.set_xlabel('P (d)')
    ax.set_ylabel(r'$\Delta$ P')

    fig.subplots_adjust(right=0.99, top=0.99, left=0.15, bottom=0.15)

    custom_savefig(fig, fig_dir, "slope_fit_" + str(n_iter).zfill(5))

    plt.close(fig)





def custom_savefig(fig, figdir, figname):
    filename = figdir + '/' + figname
    if os.path.exists(filename + '.png'):
        j = 1
        while os.path.exists(filename + str(j) + '.png'):
            j = j + 1
        filename = filename + str(j)

    fig.savefig(filename + '.png', dpi=600)


