# Import modules

import matplotlib.pyplot as plt
from bokeh.palettes import Category20c, Greens, Blues


import numpy as np
import pandas as pd

import fortranformat as ff
import pygyre as pg

import os
import subprocess
import json
import re
import vaex
from glob import glob
import itertools

from scipy.optimize import curve_fit, minimize
from scipy.signal import savgol_filter, argrelextrema
from scipy import interpolate
from natsort import natsorted

from .mesa_reader import MesaData as MR
from .run_build_star import generate_snapshot_model, find_nearest_model, generate_snapshot_name
from .jacobians import jacob_param_options
from .stellar_models import GyreMod, get_obs, Observation, PeriodSpacing
from .plotter import plot_each_timestep, plot_main_iteration_ps, plot_fit_params, plot_slope, plot_profiles


from ..utils.tools import load_parameters_from_yaml, save_parameters_to_yaml, read_gyre_model, compute_buoy_radius
from ..utils.rename_inlist import edit_inlist
from ..utils.logger import log_parameters, debug, log, log_separator, log_initial_status

import shutil



cgrav = 6.67430e-8
msun = 1.3271244e26 / 6.67430e-8
rsun = 6.957e10


inv_tool_dir = '/Users/eoin/Documents/Snapshot_Seismic/nonlinear_inversion'




class Inverter():
    def __init__(self, run_dir):
        self.run_dir = run_dir

        # Load inputs
        self.gyre_dir = self.run_dir + '/gyre_mods'
        self.static_dir = self.run_dir + '/static_mods'
        
        self.input_file = self.run_dir + '/inputs.yml'

        try:
            self.inputs = load_parameters_from_yaml(self.input_file)
        except OSError:
            raise OSError("No input file in the run directory", run_dir)

        default_inputs = load_parameters_from_yaml(inv_tool_dir + '/inversion_tool/data/inputs/defaults/inputs_default.yml')

        self.inputs.update(default_inputs)

        # Set output variable names
        self.model_name = self.inputs['model_name']
        self.output_dir = self.run_dir + '/outputs/' + self.model_name
        
        if os.path.exists(self.output_dir):
            j = 2
            while os.path.exists(self.output_dir + '_v' + str(j)):
                j = j + 1
        
            self.output_dir = self.output_dir + '_v' + str(j)

        # Make output directories if they don't already exist
        os.makedirs(self.gyre_dir, exist_ok=True)
        os.makedirs(self.static_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)

        self.fig_dir = self.output_dir
        self.log_file = self.output_dir + '/log.txt'
        self.output_file = self.output_dir + '/outputs.csv'
        self.archive_dir = self.run_dir + '/archive'

        self.build_fun = self.inputs['build_fun']
        self.m = self.inputs['m']
        self.nmode_fit = self.inputs['nmode_fit']
        self.n_iter = self.inputs['n_iteration']
        self.nmin = self.inputs['nmin']
        self.nmax = self.inputs['nmax']
        self.profile_parameters = self.inputs['profile_parameters']

        self.observation = self.inputs['observation']
        self.obs = get_obs(self.observation)
        if self.nmax == -1:
            self.nmax = len(self.obs.p)

        self.obs_vals = get_params_given_period_spacing(self.obs, self.nmin, self.nmax)
        self.obs_slope = self.obs_vals['slope']
        self.obs_intercept = self.obs_vals['intercept']
        self.obs_t = self.obs_vals['T']
        self.obs_a = self.obs_vals['A']

        self.jacobian = None
        self.b = None

        try:
            self.df_out = pd.read_csv(self.output_file, index_col='modname')
        except OSError:
            self.create_output_df()
        
        self.delta_xvals = {}

    
    def run_grid(self):
        # Runs a different mode to the inversion where we just try different parameters
        vary_params = self.inputs['vary_profile_parameters']

        ref_profile_parameters = self.profile_parameters
        
        keys = vary_params.keys()
        values = vary_params.values()


        combinations = [dict(zip(vary_params.keys(), combination)) for combination
                    in itertools.product(*vary_params.values())]
    
        model_names = []

        for combo in combinations:
            temp_p_params = ref_profile_parameters
            temp_p_params.update(combo)
            self.profile_parameters = temp_p_params

            self.do_inversion_step()
            self.profile_parameters = ref_profile_parameters
            self.previous_iterations_gyre_subdir = [self.gyre_subdir]

     
    def run(self):
        # Main mode to run the inversion
        self.do_before_inversion()
        self.do_inversion()


    def do_before_inversion(self):
        # Things to do before an inversion

        # Copy input file to output directory for future reference       
        shutil.copy(self.input_file, self.output_dir + '/inputs.yml')

        # Create output file or read from csv if already exists
        try:
            self.df_out = pd.read_csv(self.output_file, index_col='modname')
        except OSError:
            self.create_output_df()

        # Log initial parameters and status
        log_initial_status(self.log_file, self.nmode_fit, self.m)

    
    def do_inversion(self):
        # Iteratively tries to converge to best model

        # Compute model with initial parameters
        log(self.log_file, 'Computing starting model')
        self.do_inversion_step()

        self.previous_iterations_gyre_subdir = [self.gyre_subdir]

        self.n_jacob_iter = 0

        if self.inputs['plot_outputs']:
            plot_main_iteration_ps(self)
        
        log_separator(self.log_file)

        for j in self.inputs['jacobian_iterations']:

            self.jacob_option = j
            self.do_jacobian_step()
            self.previous_iterations_gyre_subdir.append(self.gyre_subdir)
            self.n_jacob_iter = self.n_jacob_iter + 1

            if self.inputs['plot_outputs']:
                plot_main_iteration_ps(self)

        log(self.log_file, "Final Parameters:", gap=False)
        log_parameters(self.profile_parameters)


    def do_inversion_step(self):
        self.create_static_dir(False)

        self.run_static_model()
        self.create_gyre_dir()
        self.run_gyre_model()
        self.append_output_to_df()

        if self.inputs['plot_outputs']:
            plot_each_timestep(self)
            plot_profiles(self)

        self.n_iter = self.n_iter + 1


    def do_after_step(self):
        pass

    def do_jacobian_step(self):

        self.compute_ref_params()
        self.compute_jacobian()
        self.compute_b_for_jacobian()
        self.solve_jacob_update_params()
        self.do_inversion_step()
        log(self.log_file, '')
        log_separator(self.log_file)


        # self.compute_ref_params()
        # self.compute_jacobian()
       
        # k_iter_max = self.inputs['k_iter_max']
        # k_iter = 0

        # large_delta_xval = True

        # while large_delta_xval and k_iter < k_iter_max:
        #     # for j in range(k_iter):
        #     # self.compute_ref_params()
        #     self.compute_b_for_jacobian()
        #     self.solve_jacob_update_params()
        #     self.do_inversion_step()
        #     log(self.log_file, '')
        #     log_separator(self.log_file)

        #     k_iter = k_iter + 1

        #     large_delta_xval = np.abs(self.delta_xvals['xcore']) > 0.001




    def get_jacobian_column(self, profile_param, ref_params, non_zero_params):
        # Computes column for Jacobian for a given profile parameter

        self.profile_param = profile_param

        # For the initial guess
        if self.jacobian is None:
            new_mass = np.maximum(0.1, ((self.obs_vals['P0'] - ref_params['P0'])/0.05) * 0.1)

            initial_delta_vals = {'rot':    0.02,
                                  'xcore':  + 0.02,
                                  'alpha':  ((self.obs_a - ref_params['A'])/1000) * 0.08,
                                  'mass':   0.02,
                                  'mabund': 0.002,
                                  'slope_ratio': 0.1}  

            delta_val = initial_delta_vals[profile_param]
        
        # In all other cases, the Jacobian from the previous timestep is used to calculate the delta value
        else:
            # Solve Jacobian with new b and get value
            self.compute_b_for_jacobian()
            x, residuals, rank, singular_values = np.linalg.lstsq(self.jacobian, self.b, rcond=None)
            param_index = self.inputs['jacobian_profile_params'].index(profile_param)
            delta_val = x[param_index]
        
        # Make sure the change in value is greater than minimum and less than some maximum
        # so we don't get stupid derivatives
        min_delta_val = self.inputs['profile_param_min_del_val'][profile_param]
        max_delta_val = self.inputs['profile_param_max_del_val'][profile_param]

        if delta_val == 0:
            delta_val = min_delta_val
        elif np.abs(delta_val) < min_delta_val:
            delta_val = min_delta_val * delta_val/np.abs(delta_val)
        elif np.abs(delta_val) > max_delta_val:
            delta_val = max_delta_val * delta_val/np.abs(delta_val)

        old_val = self.profile_parameters[profile_param] 
        new_val = old_val + delta_val
        
        # Make sure the value is >= min and <= max allowed values so we don't go to stupid regions of parameter space
        min_allowed_val = self.inputs['profile_param_min_vals'][profile_param]
        max_allowed_val = self.inputs['profile_param_max_vals'][profile_param]
        
        new_val = np.maximum(min_allowed_val, new_val)
        new_val = np.minimum(max_allowed_val, new_val)

        # Check to make sure we have not stayed on the limit from the previous timestep
        if new_val == min_allowed_val and new_val == old_val:
            new_val = min_allowed_val + min_delta_val
        elif new_val == max_allowed_val and new_val == old_val:
            new_val = max_allowed_val - min_delta_val

        # Set to floats
        x1 = round(float(old_val), 5)
        x2 = round(float(new_val), 5) 
        dx = x2 - x1

        # Log output
        log(self.log_file, "Computing " + profile_param + " component of Jacobian".ljust(43), "{:.5f}".format(x1) + ' -> ' + "{:.5f}".format(x2), gap=False)
        debug("Old " + profile_param.ljust(15), "New " + profile_param.ljust(15), "delta val".ljust(15), gap=False)
        debug(str(x1).ljust(15), str(x2).ljust(15), str(dx).ljust(15))

        # Compute derivative
        self.profile_parameters[profile_param] = new_val
        self.do_inversion_step()
        self.profile_parameters[profile_param] = ref_params[profile_param]
        fit_params = dict(self.df_out.iloc[-1])

        # Compute Jacobian column based on outputs from model
        delta_val = fit_params[profile_param] - ref_params[profile_param]

        jacob_vals = self.get_jacob_vals(fit_params, ref_params, delta_val)

        # Return dictionary with entries for Jacobian
        return_vals = {key:value if key in non_zero_params else 0 for key, value in jacob_vals.items()}

        return return_vals


    def get_jacob_vals(self, fit_params, ref_params, dx):
        jacob_vals = {}


        for param in fit_params:
            if param in ref_params:
                jacob_vals[param] = (fit_params[param] - ref_params[param])/dx

        jacob_vals['PlastP0'] = (jacob_vals['Plast'] - jacob_vals['P0'])
        jacob_vals['dn0'] = (jacob_vals['nmin1'] - jacob_vals['nmin0'])
        try:
            jacob_vals['dn1'] = (jacob_vals['nmin2'] - jacob_vals['nmin1'])
        except KeyError:
            jacob_vals['dn1'] = 0
        try:
            jacob_vals['dn2'] = (jacob_vals['nmin3'] - jacob_vals['nmin2'])
        except KeyError:
            jacob_vals['dn2'] = 0

        return jacob_vals


    def compute_ref_params(self):
        self.ref_params = dict(self.df_out.iloc[-1])


    def setup_jacobian(self):

        jacobian_profile_params = self.inputs['jacobian_profile_params']

        self.nonzero_jacobian_entries = jacob_param_options[self.jacob_option]

        for param in jacobian_profile_params:
            if param not in self.nonzero_jacobian_entries:
                self.nonzero_jacobian_entries[param] = []

        self.jacobian_profile_params = [param for param in jacobian_profile_params if len(self.nonzero_jacobian_entries[param]) > 0]

        self.jacobian_period_spacing_params = natsorted(list(set(value for sublist in self.nonzero_jacobian_entries.values() for value in sublist)))


    def compute_jacobian(self, columns=['rot', 'alpha', 'mass']):
        log(self.log_file, "Computing Jacobian Iteration", self.n_jacob_iter)
        log_parameters(self.profile_parameters)

        self.setup_jacobian()

        jacobian_columns = []

        for param in self.jacobian_profile_params:
            column = self.get_jacobian_column(param, self.ref_params, self.nonzero_jacobian_entries[param])
            jacobian_columns.append(column)

        # Convert columns to DataFrame
        df = pd.DataFrame(jacobian_columns)

        # Reorder the rows to alphabetical order for ease of viewing
        df = df[self.jacobian_period_spacing_params]

        # For logging purposes
        self.jacob_log = df.T
        self.jacob_log.columns = self.jacobian_profile_params

        jacob = df.values.T
        self.jacobian = jacob
        self.jacobian_period_spacing_params = self.jacobian_period_spacing_params


    def compute_b_for_jacobian(self):

        jacob_vals = self.get_jacob_vals(self.obs_vals, self.ref_params, 1)

        b = [jacob_vals[key] for key in self.jacobian_period_spacing_params]


        self.b = b

    def log_jacobian(self):
        log(self.log_file, "")
        log(self.log_file, "Jacobian Matrix:", gap=False)
        log(self.log_file, self.jacob_log)

        vals_mod = [self.obs_vals[key] for key in self.jacobian_period_spacing_params]
        vals_ref = [self.ref_params[key] for key in self.jacobian_period_spacing_params]

        bvector_print = np.column_stack((vals_ref, vals_mod, self.b))
        dfb = pd.DataFrame(bvector_print, index=self.jacobian_period_spacing_params, columns=['Observation', 'Model', 'Difference'])
        log(self.log_file, "Difference in period parameters:", gap=False)
        log(self.log_file, dfb)

        log(self.log_file, self.df_diff_profile_params)
        log(self.log_file, "")

        # debug("Residuals", residuals, gap=False)
        # debug("Rank", rank)




    def solve_jacob_update_params(self):
        # Use least squares to find the solution

        self.x, residuals, rank, singular_values = np.linalg.lstsq(self.jacobian, self.b, rcond=None)

        table_row = []

        for delta_val, param in zip(self.x, self.inputs['jacobian_profile_params']):
            self.delta_xvals[param] = delta_val
            old_val = self.profile_parameters[param]
            new_val = old_val + delta_val
            
            min_allowed_val = self.inputs['profile_param_min_vals'][param]
            max_allowed_val = self.inputs['profile_param_max_vals'][param]
            
            new_val = np.maximum(min_allowed_val, new_val)
            new_val = np.minimum(max_allowed_val, new_val)
            
            self.profile_parameters[param] = new_val

            table_row.append([new_val, old_val, delta_val])

        self.df_diff_profile_params = pd.DataFrame(table_row, columns=[ 'New', 'Old', 'Change'], index=self.jacobian_profile_params)

        self.log_jacobian()


        # log(self.log_file, "new " + param + ":".ljust(15), "{:.5f}".format(new_val), gap=False)


    def create_static_dir(self, overwrite=True):
        starting_model = self.get_starting_model()
        self.set_modname()

        if not os.path.isfile(self.mod_path) or overwrite:
            # log(self.log_file, "Creating MESA subdirectory for static model", self.modname)

            modname = generate_snapshot_model(starting_model, self.static_dir, self.build_fun, 
                                              change_mass_flag=True, **self.snapshot_kwargs)

            self.create_inlists()
        else:
            # log(self.log_file, "MESA subdirectory already exists", self.modname)
            pass



    def set_modname(self):

        self.snapshot_kwargs = {key: value for key, value in self.profile_parameters.items() if key != "rot"}

        mod_path = generate_snapshot_name(self.static_dir, self.snapshot_kwargs)

        self.mod_path = mod_path
        
        self.modname = self.mod_path.split('/')[-1][:-4]

        self.static_subdir = self.static_dir + '/' + self.modname



    def get_starting_model(self):
        if True or len(self.df_out) == 0:
            database_dir =  "/Users/eoin/Documents/Snapshot_Seismic/models/snapshot_database"
            lookup_table_dir = database_dir + '/lookup_table'

            params = self.profile_parameters
            starting_model = find_nearest_model(params['mass'],
                                        params['zval_lookup'],
                                        params['xcore'],
                                        lookup_table_dir)
            return starting_model

    def get_best_model(self):
        df = self.df_out
        return df.loc[df['chi_sq'].idxmin()]


    def create_inlists(self):
        debug("Editing inlists", self.modname)

        # Copy template directory to subdirectory for static model
        template = inv_tool_dir + '/inversion_tool/data/inputs/template_buildastar/static_template'
        shutil.copytree(template, self.static_subdir, dirs_exist_ok=True)

        # Copy model into new subdirectory for static model
        new_modpath = self.static_subdir + '/' + self.modname + '.mod'

        if not os.path.exists(new_modpath):
            shutil.copy(self.mod_path, new_modpath)

        dotmod_name = "'" + self.modname + ".mod'"
        list_dir = inv_tool_dir + '/inversion_tool/data/inputs/template_buildastar'
        hist_file = list_dir + '/history_columns.list'
        prof_file = list_dir + '/profile_columns.list'

        template_inlist = template + '/inlist_static'
        modified_inlist = self.static_subdir + '/' + 'inlist_static'
        dict_changes = {
            'load_model_filename': dotmod_name,
            'history_columns_file': hist_file,
            'profile_columns_file': prof_file
        }

        quotes = ['history_columns_file', 'profile_columns_file']

        edit_inlist(template_inlist, dict_changes, modified_inlist, quotes=quotes)


    def run_static_model(self):
        # use subprocess.run() to execute the script and wait for it to complete
        run_script = self.static_subdir + '/rew'
        cwd = self.static_subdir

        if not os.path.exists(self.static_subdir + '/static.data.GYRE'):
            debug("Running MESA static model ...", gap=False)
            process = subprocess.run(['bash', run_script], cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if process.returncode == 0:
                debug("Successfully computed MESA static model")
        else:
            debug("MESA static model already computed")


    def set_gyre_subdir(self):
        rot_val = self.profile_parameters['rot']
        self.gyre_subdir = self.gyre_dir + '/' + self.modname + '_r' + str(int(1000*rot_val)).zfill(4)

    def get_gyre_subdir(self, modname, rot_val):
        return self.gyre_dir + '/' + modname + '_r' + str(int(1000*rot_val)).zfill(4)


    def create_gyre_dir(self, overwrite=True):
        self.set_gyre_subdir()

        if not os.path.exists(self.gyre_subdir + '/evol.mesa') or overwrite:

            template = inv_tool_dir + '/inversion_tool/data/inputs/template_buildastar/gyre_template'

            shutil.copytree(template, self.gyre_subdir, dirs_exist_ok=True)

            dict_changes = {
                'coriolis_method': 'TAR',
                'Omega_rot': self.profile_parameters['rot'],
                'm': self.m
            }

            quotes = ['coriolis_method']

            inlist_name = self.gyre_subdir + '/gyre.in'

            debug('Creating GYRE Inlists for ' + self.modname, 'rot', dict_changes['Omega_rot'])

            edit_inlist(inlist_name, dict_changes, inlist_name, quotes=quotes)


            static_mod_path = self.static_subdir + '/static.data.GYRE'
            evol_mesa_path = self.gyre_subdir + '/evol.mesa'

            if not os.path.exists(evol_mesa_path):
                shutil.copy(static_mod_path, evol_mesa_path)


    def run_gyre_model(self):

        self.set_gyre_subdir()

        run_script = self.gyre_subdir + '/run_gyre'
        cwd = self.gyre_subdir

        if not os.path.exists(self.gyre_subdir + '/summary.h5'):

            debug("Running GYRE Model ...", gap=False)
            process = subprocess.run(['bash', run_script], cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if process.returncode == 0:
                debug("Successfully completed GYRE model")



    def append_output_to_df(self):

        rot_val = self.profile_parameters['rot']
        mod_saved = self.modname + '_r' + str(int(1000*rot_val)).zfill(4)
        if True or mod_saved not in list(df['modname'].values):

            debug("Saving to dataframe")

            mod = GyreMod(self.gyre_subdir, m=self.m)

            ps = mod.set_dp_for_nmin_nmax(self.inputs['nmode_fit'], self.inputs['nmode_fit'] + len(self.obs.dp) - 1)

            nmin = self.nmin
            nmax = self.nmax

            fit_vals = get_params_given_period_spacing(ps, nmin, nmax)

            obs = self.obs
            obs_vals = self.obs_vals

            diff_vals = {}
            for fit_name in fit_vals:
                try:
                    diff_vals['diff_' + fit_name] = obs_vals[fit_name] - fit_vals[fit_name]
                except:
                    diff_vals['diff_' + fit_name] = np.NaN

            # chi_sq = np.sum((ps.p[nmin:nmax] - obs.p[nmin:nmax])**2 + (ps.dp[nmin:nmax] - obs.dp[nmin:nmax])**2)

            row = {}
            row['modname'] = self.modname + '_r' + str(int(1000*rot_val)).zfill(4)

            row['n_iteration'] = self.n_iter

            row.update(self.profile_parameters)
            row.update(fit_vals)
            row.update(diff_vals)
            row['nrandom'] = 0

            try:
                df = pd.read_csv(self.output_file)

                debug("Saving to output csv")
                df = pd.concat([df, pd.DataFrame(row, index=[0])], ignore_index=True)
            except OSError:
                df = pd.DataFrame(row)

            df.to_csv(self.output_file, index=False)

            self.df_out = pd.read_csv(self.output_file, index_col='modname')


        
    def create_output_df(self):
        if not os.path.isfile(self.output_file):
            struct_params = list(self.inputs['profile_parameters'].keys())
            other_cols = ['chi_sq', 'nrandom']
            name_cols = ['modname', 'n_iteration']

            df_cols = name_cols + struct_params + other_cols


            df = pd.DataFrame([], columns=df_cols)
            df = df.set_index('modname')

            df.to_csv(self.output_file)

            self.df_out = df





def get_params_given_period_spacing(ps, nmin, nmax):

    slope, intercept = ps.fit_slope()

    df = ps.get_fit_params(nmin, nmax)

    params = {'slope': round(float(slope), 0),
              'intercept': round(float(intercept), 0),
              'T': round(float(df['P'].values[0]), 4),
              'A': round(float(df['A'].values[0]), 0),
              'Plast': round(float(ps.p[-1]), 5)}

    for i in range(len(ps.p)):
        params['P' + str(i)] = round(float(ps.p[i]), 5)
        params['dP' + str(i)] = round(float(ps.dp[i]), 5)
        params['dP' + str(i) + '_rot'] = round(float(ps.dp_rot[i]), 5)


    dp = ps.dp_rot

    jmins = argrelextrema(dp, np.less)[0]

    jfits = []


    nmin = 0

    for jmin in jmins:
        
        r1 = (dp[jmin + 1] - dp[jmin])/dp[jmin]
        r2 = (dp[jmin - 1] - dp[jmin])/dp[jmin]
        ratio = r1/r2
        if ratio < 1:
            fratio = np.cos(r1/r2 * (np.pi/2))
        elif ratio >= 1:
            fratio = -np.cos(r2/r1 * (np.pi/2))

        jfit = jmin + 0.5 * fratio
        params['nmin' + str(nmin)] = jfit
        nmin = nmin + 1

    
    for j in np.arange(1, len(dp) - 2):
        
        r1 = (dp[j + 1] - dp[j])/dp[j]
        r2 = (dp[j - 1] - dp[j])/dp[j]
        ratio = r1/r2
        if ratio < 1:
            fratio = np.cos(r1/r2 * (np.pi/2))
        elif ratio >= 1:
            fratio = -np.cos(r2/r1 * (np.pi/2))

        jfit = j + 0.5 * fratio
        params['jshape' + str(j)] = jfit

 
    return params



def create_evol_model(mass, z, min_D_mix, xc_min, moddir, modname):
    # Creates evolutionary model given M, Z, min_D_mix

    # Copy template directory to subdirectory for static model
    template = inv_tool_dir + '/inversion_tool/data/inputs/template_buildastar/static_template'
    new_moddir = moddir + '/' + modname
    shutil.copytree(template, new_moddir, dirs_exist_ok=True)

    list_dir = inv_tool_dir + '/inversion_tool/data/inputs/template_buildastar'
    hist_file = list_dir + '/history_columns.list'
    prof_file = list_dir + '/profile_columns.list'

    prems_inlist = template + '/inlist_prems'
    ms_inlist = template + '/inlist_ms'

    dict_changes = {
        'initial_mass': mass,
        'initial_z': z,
        'min_D_mix': min_D_mix,
        'history_columns_file': hist_file,
        'profile_columns_file': prof_file,
        'xa_central_lower_limit(1)': xc_min
    }

    quotes = ['history_columns_file', 'profile_columns_file']

    edit_inlist(prems_inlist, dict_changes, prems_inlist, quotes=quotes)
    edit_inlist(ms_inlist, dict_changes, ms_inlist, quotes=quotes)

    shutil.copy(new_moddir + '/rn_evol', new_moddir + '/rn')



def create_static_model(dotmod_file_path, new_moddir):
    # Creates evolutionary model given M, Z, min_D_mix

    # Copy template directory to subdirectory for static model
    template = inv_tool_dir + '/inversion_tool/data/inputs/template_buildastar/static_template_adiabatic_pz'
    # new_moddir = moddir + '/' + subdir
    shutil.copytree(template, new_moddir, dirs_exist_ok=True)

    list_dir = inv_tool_dir + '/inversion_tool/data/inputs/template_buildastar'
    hist_file = list_dir + '/history_columns.list'
    prof_file = list_dir + '/profile_columns.list'

    new_static_inlist = new_moddir + '/inlist_static'

    modname = dotmod_file_path.split('/')[-1]

    dict_changes_static = {
        'load_model_filename': modname,
        'history_columns_file': hist_file,
        'profile_columns_file': prof_file}

    quotes = ['load_model_filename', 'history_columns_file', 'profile_columns_file']

    edit_inlist(new_static_inlist, dict_changes_static, new_static_inlist, quotes=quotes)

    shutil.copy(dotmod_file_path, new_moddir + '/' + modname) 



def create_gyre_tests(evol_mesa_file, moddir, subdir, omega=0.3, m=1):

    template = inv_tool_dir + '/inversion_tool/data/inputs/template_buildastar/gyre_template'

    gyre_dir = moddir + '/' + subdir

    shutil.copytree(template, gyre_dir, dirs_exist_ok=True)

    dict_changes = {
        'coriolis_method': 'TAR',
        'Omega_rot': omega,
        'm': m
    }

    quotes = ['coriolis_method']

    inlist_name = gyre_dir + '/gyre.in'

    edit_inlist(inlist_name, dict_changes, inlist_name, quotes=quotes)

    evol_mesa_path = gyre_dir + '/evol.mesa'

    if not os.path.exists(evol_mesa_path):
        shutil.copy(evol_mesa_file, evol_mesa_path)




if __name__ == '__main__':
    pass
 