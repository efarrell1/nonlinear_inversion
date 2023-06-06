import pygyre as pg
from .mesa_reader import MesaData as MR

import numpy as np
import pandas as pd
import fortranformat as ff

import os
from glob import glob

from scipy.optimize import curve_fit, minimize
from scipy.signal import savgol_filter, argrelextrema
from scipy import interpolate

import subprocess

import json
import re
import vaex

from ..utils.tools import load_parameters_from_yaml, save_parameters_to_yaml, read_gyre_model, compute_buoy_radius


def convert_pedersen_data_to_csv():
    df_list = []

    for data_file in data_files:
        if 'spacings' in data_file:
            df = pd.read_csv(data_file, delimiter='\t')
        else:
            df = pd.read_csv(data_file, delimiter='\t')
            df['per'] = 1/df['freq']

        name = data_file.split('KIC')[-1].split('_')[0]
        df['KIC'] = name
        df_list.append(df)


    df_all = pd.concat(df_list, ignore_index=True)
    df_all.to_csv('/Users/eoin/Documents/Snapshot_Seismic/data/pedersen_data.csv')


def get_obs_df_kic(kic_id):
    csv_filename = '/Users/eoin/Documents/Snapshot_Seismic/data/pedersen_data.csv'
    df = pd.read_csv(csv_filename)
    dfnew = df[df['KIC'] == kic_id].copy()
    return dfnew




def get_obs(kic='3459297'):
    if type(kic) == int:
        return get_obs2(kic)
    else:
        inv_tool_dir = '/Users/eoin/Documents/Snapshot_Seismic/inversion_tool'

        moddir = inv_tool_dir + '/inversion_tool/data/inputs/observations/kepler_' + str(kic) + '.csv'

        obs = Observation(moddir)
        obs.color='black'
        obs.label = kic
        obs.buoy_r = np.array([])
        obs.N_plot = np.array([])

        obs.mass = np.array([])
        obs.h1 = np.array([])

        return obs


def get_kic_name(kic_id):
    input_file_dir = '/Users/eoin/Documents/Snapshot_Seismic/inversion_tool/inversion_tool/data/inputs/observations'

    file_path = input_file_dir + '/' + str(int(kic_id)) + '.json'

    return file_path



def get_obs2(kic_id):
    file_path = get_kic_name(kic_id)

    with open(file_path, "r") as file:
        data = json.load(file)   

    p = data['per']
    dp = np.diff(p)*86400
    obs = Observation2(p[:-1], dp)

    obs.color='black'
    obs.label = kic_id
    obs.buoy_r = np.array([])
    obs.N_plot = np.array([])

    obs.mass = np.array([])
    obs.h1 = np.array([])

    return obs





class Observation():
    def __init__(self, moddir):
        self.moddir = moddir
        df = pd.read_csv(self.moddir, header=None)
        
        try:
            df.columns = ['p', 'dp']
        except:
            df.columns = ['n', 'p', 'f', 'a', 'theta', 'sn']
            df['dp'] = np.diff(df['p'].values)*86400
        
        self.p = df['p'].values
        self.dp = df['dp'].values
        self.periods = self.p

        self.nmode = np.arange(len(self.dp))

        self.set_period_spacing()

    def set_period_spacing(self):
        self.period_spacing = PeriodSpacing(self.dp, self.periods, nmode=self.nmode)

    def set_nmode(self, n_match, p_match, i_match=0):
        istart = np.argmin(np.abs(p_match - self.p[i_match]))
        nstart = n_match[istart]
        self.nmode = np.arange(nstart, len(self.dp) + nstart, 1)
        self.set_period_spacing()


    def __getattr__(self, method_name):
        try:
            return getattr(self.period_spacing, method_name)
        except AttributeError:
            raise AttributeError("could not find parameter")    





class Observation2():
    def __init__(self, p, dp):
        # self.moddir = moddir
        # df = pd.read_csv(self.moddir, header=None)
        
        # try:
        #     df.columns = ['p', 'dp']
        # except:
        #     df.columns = ['n', 'p', 'f', 'a', 'theta', 'sn']
        #     df['dp'] = np.diff(df['p'].values)*86400

        self.p = p
        self.dp = dp
        self.periods = self.p

        self.nmode = np.arange(len(self.dp))

        self.set_period_spacing()

    def set_period_spacing(self):
        self.period_spacing = PeriodSpacing(self.dp, self.periods, nmode=self.nmode)

    def set_nmode(self, n_match, p_match, i_match=0):
        istart = np.argmin(np.abs(p_match - self.p[i_match]))
        nstart = n_match[istart]
        self.nmode = np.arange(nstart, len(self.dp) + nstart, 1)
        self.set_period_spacing()


    def __getattr__(self, method_name):
        try:
            return getattr(self.period_spacing, method_name)
        except AttributeError:
            raise AttributeError("could not find parameter")    



class StellarStructure():
    def __init__(self, pfdir, rot=None):
        self.pfdir = pfdir
        self.pf = MR(self.pfdir)

        self.rot = rot

        gyre_dirs = glob(self.pfdir[:-5] + '_r*')

        self.gyre_dirs = sorted([x for x in gyre_dirs if os.path.isdir(x)])

        try:
            self.rot_vals = [int(x[-4:])/1e3 for x in gyre_dirs]
        except:
            self.rot_vals = None

        if rot is not None:
            self.set_gyre_mod(rot=rot)

    def set_gyre_mod(self, rot):
        moddir = self.pfdir[:-5] + '_r' + str(int(1000 * rot)).zfill(4)
        self.gyre_mod = GyreMod(moddir)

    def get_all_gyre_mods(self):
        moddirs = [self.pfdir[:-5] + '_r' + str(int(1000 * rot)).zfill(4)
                   for rot in self.rot_vals]
        gyre_mods = [GyreMod(moddir) for moddir in moddirs]

        return gyre_mods


    def __getattr__(self, method_name):
        try:
            return getattr(self.pf, method_name)
        except AttributeError:
            pass

        try:
            return getattr(self.gyre_mod, method_name)
        except AttributeError:
            raise AttributeError("could not find parameter")    

def fquartic(x, a, b, c, d, e):
    return e * x ** 4 + d * x ** 3 + a * x ** 2 + b * x + c

def flinear(x, b, c):
    return  b * x + c

def fquadratic(x, a, b, c):
    return  a*x**2 + b * x + c

def fcubic(x, a, b, c, d):
    return  d * x ** 3 + a*x**2 + b * x + c


def fcubic(x, a, b, c, d):
    return  d * x ** 3 + a*x**2 + b * x + c

def fit_abs_sine(x, period, amplitude, phi, yoffset):
    return amplitude * np.abs(np.sin((np.pi * (x - x[0]) / period) + phi)) + yoffset

def fit_sine(x, period, amplitude, phi, yoffset):
    return 0.5*amplitude * np.sin((np.pi * (x) * 2 / period) + phi) + yoffset

def fit_abs_sine_slope(x, period, amplitude, phi, yoffset, slope):
    return amplitude * np.abs(np.sin((np.pi * (x - x[0]) / period) + phi)) + yoffset + slope * x

def objective_abs_sine(params, x, y):
    A, B, C, D = params
    y_fit = fit_abs_sine(x, A, B, C, D)
    residuals = y_fit - y
    return np.sum(residuals**2)


def objective_abs_sine_slope(params, x, y):
    A, B, C, D, E = params
    y_fit = fit_abs_sine_slope(x, A, B, C, D, E)
    residuals = y_fit - y
    return np.sum(residuals**2)



class PeriodSpacing():
    def __init__(self, dp, p, nmode=None):
        self.dp = dp
        self.p = p

        if nmode is None:
            self.nmode = np.arange(len(dp))
        else:
            self.nmode = nmode

        self.frot = fcubic
        self.fit_fun = fit_abs_sine
        self.objective = objective_abs_sine


        self.dp_fit = self.dp
        self.nmode_fit = self.nmode

        self.correct_for_rotation_slope()

        self.set_pguess()


    def correct_for_rotation_slope(self, xdata='nmode'):
        # takes period spacing pattern of rotating model
        # and corrects for the slope

        if xdata == 'periods':
            xvals = self.p
        elif xdata == 'nmode':
            xvals = self.nmode
           
        self.popt_rot, self.pcov_rot = curve_fit(self.frot, np.array(xvals), 1/np.array(self.dp))

        self.dp_rot = (self.dp * self.frot(xvals, *self.popt_rot))*3e3


    def fit_slope(self):

        # print(self.p, self.dp)

        self.popt_slope, self.pcov_slope = curve_fit(flinear, self.p, np.array(self.dp))

        return [-self.popt_slope[0], self.popt_slope[1]]


    def set_pguess(self):
        self.maxima_indices = argrelextrema(self.dp_rot, np.greater)[0]
        if len(self.maxima_indices) != 0:
            self.p_guess = len(self.dp_rot)/len(self.maxima_indices)
        else:
            self.p_guess = 5


    def quick_fit(self):
        # Quickly fit dp vs. nmode
        self.set_guess_params()
        self.compute_fit()


    def set_main_period(self):
        self.quick_fit()
        self.pmain = self.result.x[0]


    def set_guess_params(self, phi_guess=3.0, p_guess=None):

        # self.set_pguess(self.dp_fit)
       
        # Amplitude:
        # a_guess = np.abs(1.3 * (np.max(self.dp_fit) - np.min(self.dp_fit)))
        # amin = a_guess * 0.5
        # amax = a_guess * 1.3

        a_guess = np.abs(1.3 * (np.max(self.dp_fit) - np.min(self.dp_fit)))
        amin = a_guess * .9
        amax = a_guess * 1.5
        
        # Offset:
        if self.fit_fun == fit_abs_sine:
            c_guess = np.abs(np.min(self.dp_fit)*0.9)
            cmin = c_guess * 0.8
            cmin = c_guess * 1
            cmax = np.abs(np.min(self.dp_fit))*.95
        elif self.fit_fun == fit_sine: 
            c_guess = np.mean(self.dp_fit)
            cmin = c_guess * 0.9
            cmax = c_guess * 1.8
        else:
            c_guess = np.abs(np.min(self.dp_fit)*0.9)
            cmin = c_guess * 0.8
            cmax = np.min(self.dp_fit)*.95

        # Period:
        if p_guess is None:
            p_guess = self.p_guess

        pmin = p_guess * 0.85
        pmax = p_guess * 1.4

        if self.objective == objective_abs_sine:
            self.params0 = [self.p_guess, a_guess, phi_guess, c_guess]
            self.bounds = [(pmin, pmax), (amin, amax), (None, None), (cmin, cmax)]
        elif self.objective == objective_abs_sine_slope:
            self.params0 = [p_guess, a_guess, phi_guess, c_guess, -100]
            self.bounds = [(pmin, pmax), (amin, amax), (None, None), (cmin, cmax), (-1000, 100)]


                    
    def compute_fit(self):
        try:
            self.result = minimize(self.objective, self.params0, args=(self.nmode_fit, self.dp_fit), bounds=self.bounds)
        except ValueError:
            print(self.bounds)
            return 
        
        self.resid = np.sum((self.fit_fun(self.nmode_fit, *self.result.x) - self.dp_fit)**2)

        self.params_fit = self.result.x


    def phase_cycle(self):
        
        phi_guesses = np.arange(1, 5)
        
        resid_list = []
        params_fit_list = []
        
        for phi_guess in phi_guesses:
            self.set_guess_params(phi_guess=phi_guess)
            self.compute_fit()

            resid_list.append(self.resid)
            params_fit_list.append(self.params_fit)
        
        ibest = np.argmin(resid_list)
        self.params_best = params_fit_list[ibest]


    def phase_period_cycle(self):
        
        phi_guesses = np.arange(1, 5)
        p_guesses = np.arange(2.1, 2.2, 0.01)
        
        resid_list = []
        params_fit_list = []
        
        for phi_guess in phi_guesses:
            for p_guess in p_guesses:
                self.set_guess_params(phi_guess=phi_guess, p_guess=p_guess)
                self.compute_fit()

                resid_list.append(self.resid)
                params_fit_list.append(self.params_fit)
            
        ibest = np.argmin(resid_list)
        self.params_best = params_fit_list[ibest]

        a_guess = np.abs(1.3 * (np.max(self.dp_fit) - np.min(self.dp_fit)))
        amin = a_guess * .9
        amax = a_guess * 1.5

        c_guess = np.abs(np.min(self.dp_fit)*0.9)
        cmin = c_guess * 1
        cmax = np.abs(np.min(self.dp_fit))*.95

        # print("a_guess", a_guess, amin, amax)
        # print("c_guess", c_guess, cmin, cmax)




    def set_covariance(self):
        self.set_guess_params(phi_guess=self.params_best[2])
        params0 = self.params_best
        bounds = ([x[0] for x in self.bounds], [x[1] for x in self.bounds])
        bounds[0][2] = -100
        bounds[1][2] = 100
        
        popt, pcov = curve_fit(self.fit_fun, self.nmode_fit, self.dp_fit, bounds=bounds, p0=params0)
        
        self.popt_cfit = popt
        self.pcov = np.sqrt(np.abs(np.diag(pcov)))
     
        return self.popt_cfit, self.pcov

    def limit_period_spacing(self, nmin=None, nmax=None, rot=True):
        # Returns nmode, dp for nmode >= nmin and <= nmax
        if rot:
            dp = self.dp_rot
        else:
            dp = self.dp

        if nmin is not None:
            dp = np.array([d for d, n in zip(dp, self.nmode) if n >= nmin])
            nmode = np.array([n for n in self.nmode if n >= nmin])
        
        if nmax is not None:
            dp = np.array([d for d, n in zip(dp, nmode) if n <= nmax])
            nmode = np.array([n for n in nmode if n <= nmax])
        
        self.nmode_fit = nmode
        self.dp_fit = dp

        return self.nmode_fit, self.dp_fit


    def get_fit_params(self, nmin, nmax, cov_flag=False, rot=True):
        self.limit_period_spacing(nmin=nmin, nmax=nmax, rot=rot)

        # if len(self.nmode_fit) != nmax - nmin + 1:
        #     print('error - likely missing modes somewhere', 'nmin', nmin, 'nmax', nmax)
        #     return None

        # Sets self.params_best
        self.phase_cycle()

        df_columns = ['n', 'P', 'A', 'Phi', 'c']
        if self.objective == objective_abs_sine_slope:
            df_columns.append('m')

        if not cov_flag:
            df_vals = np.concatenate(([nmin], self.params_best))
        elif cov_flag:
            if self.objective == objective_abs_sine:
                df_columns = df_columns + ['P2', 'A2', 'Phi2', 'c2', 'P_err', 'A_err', 'Phi_err', 'c_err']
            elif self.objective == objective_abs_sine_slope:
                df_columns = df_columns + ['P2', 'A2', 'Phi2', 'c2', 'm2', 'P_err', 'A_err', 'Phi_err', 'c_err', 'm_err']
            self.set_covariance()
            df_vals = np.concatenate(([nmin], self.params_best, self.popt_cfit, self.pcov))

        data = {column: val for column, val in zip(df_columns, df_vals)}
        df = pd.DataFrame(data, index=[0])

        return df

    def get_fit_params2(self, nmin, nmax, cov_flag=False, rot=True):
        self.limit_period_spacing(nmin=nmin, nmax=nmax, rot=rot)

        # if len(self.nmode_fit) != nmax - nmin + 1:
        #     print('error - likely missing modes somewhere', 'nmin', nmin, 'nmax', nmax)
        #     return None

        # Sets self.params_best
        self.phase_period_cycle()

        df_columns = ['n', 'P', 'A', 'Phi', 'c']
        if self.objective == objective_abs_sine_slope:
            df_columns.append('m')

        if not cov_flag:
            df_vals = np.concatenate(([nmin], self.params_best))
        elif cov_flag:
            if self.objective == objective_abs_sine:
                df_columns = df_columns + ['P2', 'A2', 'Phi2', 'c2', 'P_err', 'A_err', 'Phi_err', 'c_err']
            elif self.objective == objective_abs_sine_slope:
                df_columns = df_columns + ['P2', 'A2', 'Phi2', 'c2', 'm2', 'P_err', 'A_err', 'Phi_err', 'c_err', 'm_err']
            self.set_covariance()
            df_vals = np.concatenate(([nmin], self.params_best, self.popt_cfit, self.pcov))

        data = {column: val for column, val in zip(df_columns, df_vals)}
        df = pd.DataFrame(data, index=[0])

        return df

    def get_fit_params_slope(self, nmin, nmax, cov_flag=False, rot=False):
        self.limit_period_spacing(nmin=nmin, nmax=nmax, rot=rot)

        # if len(self.nmode_fit) != nmax - nmin + 1:
        #     print('error - likely missing modes somewhere', 'nmin', nmin, 'nmax', nmax)
        #     return None

        # Sets self.params_best
        self.phase_cycle()

        df_columns = ['n', 'P', 'A', 'Phi', 'c']
        if self.objective == objective_abs_sine_slope:
            df_columns.append('m')

        if not cov_flag:
            df_vals = np.concatenate(([nmin], self.params_best))
        elif cov_flag:
            if self.objective == objective_abs_sine:
                df_columns = df_columns + ['P2', 'A2', 'Phi2', 'c2', 'P_err', 'A_err', 'Phi_err', 'c_err']
            elif self.objective == objective_abs_sine_slope:
                df_columns = df_columns + ['P2', 'A2', 'Phi2', 'c2', 'm2', 'P_err', 'A_err', 'Phi_err', 'c_err', 'm_err']
            self.set_covariance()
            df_vals = np.concatenate(([nmin], self.params_best, self.popt_cfit, self.pcov))

        data = {column: val for column, val in zip(df_columns, df_vals)}
        df = pd.DataFrame(data, index=[0])

        # print(self.dp)

        return df


    def get_variation_in_fit_params(self, nmin=None, nmax=None, ngap=10, cov_flag=False, rot=True):
        
        self.set_main_period()

        if nmin is None:
            nmin = np.min(self.nmode)
        if nmax is None:
            nmax = np.max(self.nmode)
        
        nmin_list = np.arange(nmin, nmax - ngap, 1)
        nmax_list = np.array([x + ngap - 1 for x in nmin_list])
        
        df_list = [self.get_fit_params(nmin, nmax, cov_flag=cov_flag, rot=rot) for nmin, nmax in zip(nmin_list, nmax_list)]
    
        df = pd.concat(df_list)

        self.df_fit_params = df

        return df



class GyreMod():
    def __init__(self, moddir, summary_name='summary.h5', detail_name='detail', read_mod=True, m=1):

        self.moddir = moddir
        self.table_summary = pg.read_output(self.moddir + '/' + summary_name)
        self.df = self.table_summary.to_pandas()

        self.detail_name = detail_name

        self.m = m

        # Loads hydrostatic structure model
        if read_mod:
            mod, mod_header = read_gyre_model(self.moddir + '/evol.mesa')

            self.mod = mod
            self.N2 = self.mod['N2']
            self.N_plot = np.sqrt(np.abs(self.N2)) * 1e3
            self.mod_header = mod_header

            self.compute_buoyancy_radius()
        else:
            self.mod = self.df
        
        self.set_dp()

        self.period_spacing()

        self.gyre_file = self.moddir + '/gyre.in'

        inputFile = open(self.gyre_file, "r+")
        inFile = inputFile.readlines()
        inputFile.close()

        infile_split = [line.split(' = ') for line in inFile]
        mesa_params = [l[0] for l in infile_split]
        mesa_params = [x.strip() for x in mesa_params]
        mesa_vals = [l[-1] for l in infile_split]
        mesa_vals = [x.strip() for x in mesa_vals]

        dinputs = dict(zip(mesa_params, mesa_vals))

        self.rot = float(dinputs['Omega_rot'])


    def compute_buoyancy_radius(self):
        N2_vals = self.mod['N2']
        r_vals = self.mod['radius']
        br_profile, br_mag = compute_buoy_radius(N2_vals, r_vals)
        self.buoy_r = br_profile
        self.bouy_r_mag = br_mag

        return self.buoy_r
    

    def pf(self, n, l=1, m=1):
        if m >= 0:
            m = '+' + str(m)
        else:
            m = '-' + str(m)

        if n > 0:
            n = '+' + str(n)
        else:
            n = str(n)

        try:
            return pg.read_output(self.moddir + '/' + self.detail_name + '.l' + str(l) + '.n' + n + '.h5')
        except:
            return pg.read_output(self.moddir + '/' + self.detail_name + '.l' + str(l) + '.n' + n + '.m' + str(m) + '.h5')

    def set_dp(self):

        df =  self.df.query('l == 1 & n_p == 0')
        df = df[df['m'] == self.m]
        periods = 1/np.real(df['freq'])
        dp = -1 * np.diff(periods) * 86400

        self.dp = dp[::-1]
        self.nmode = df['n_g'].values[1:][::-1]

        self.periods = periods[1:][::-1]

    def set_dp_for_obs(self, obs, i_match=0):
        # Sets Dp, nmode, periods given observation class
        
        istart = np.argmin(np.abs(obs.p[i_match] - self.p))
        imax = np.minimum(istart + len(obs.p), len(self.p))
        if imax - istart > 7:
            dp_new = self.dp[istart:imax]
            periods_new = self.p[istart:imax]
            self.ps_obs = PeriodSpacing(dp_new, periods_new)

            return self.ps_obs

        else:
            self.ps_obs = self.period_spacing
            return self.ps_obs


    def set_dp_for_nmin_nmax(self, nmin, nmax):

        istart = np.where(self.nmode == nmin)[0][0]
        iend = np.where(self.nmode == nmax)[0][0]

        p = self.p[istart:iend+1]
        dp = self.dp[istart:iend+1]

        return PeriodSpacing(dp, p)


    def period_spacing(self):
        self.set_dp()
        self.period_spacing = PeriodSpacing(self.dp, self.periods, nmode=self.nmode)

        return self.period_spacing



    def __getattr__(self, method_name):
        try:
            return getattr(self.mod, method_name)
        except AttributeError:
            pass

        try:
            return getattr(self.period_spacing, method_name)
        except AttributeError:
            pass
    
        try:
            return getattr(self.df, method_name)
        except AttributeError:
            raise AttributeError("could not find parameter")


