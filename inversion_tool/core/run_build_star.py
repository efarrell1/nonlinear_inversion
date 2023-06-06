import numpy as np

import itertools
import subprocess
import json

from .buildastar import build_new_h_gradient, change_mass, profileplot

from .mesa_reader import MesaData as MR

import glob
import subprocess

import pandas as pd
import os


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def generate_snapshot_name(run_dir, kwargs):

    # Generate filename based on inputs values
    names = []

    if 'mass' in kwargs:
        mass_name = 'm' + str(int(1000 * kwargs['mass'])).zfill(4)
        names.append(mass_name)
    
    if 'dx_core' in kwargs:
        xcore_name = 'dx' + str(int(1000 * kwargs['dx_core'])).zfill(4)
        names.append(xcore_name)
    
    if 'xcore' in kwargs:
        xcore_name = 'xc' + str(int(1000 * kwargs['xcore'])).zfill(4)
        names.append(xcore_name)

    if 'mcore' in kwargs:
        core_name = 'c' + str(int(1000 * kwargs['mcore'])).zfill(4)
        names.append(core_name)

    if 'mabund' in kwargs:
        core_name = 'ma' + str(int(1000 * kwargs['mabund'])).zfill(4)
        names.append(core_name)

    if 'mslope' in kwargs:
        slope_name = 's' + str(int(1000 * kwargs['mslope'])).zfill(4)
        names.append(slope_name)
    
    if 'slope_ratio' in kwargs:
        core_name = 's' + str(int(1000 * kwargs['slope_ratio'])).zfill(4)
        names.append(core_name)

    if 'alpha' in kwargs:
        xsurface_name = 'a' + str(int(10000 * kwargs['alpha'])).zfill(4)
        names.append(xsurface_name)

    if 'fa' in kwargs:
        xsurface_name = 'fa' + str(int(1000 * kwargs['fa'])).zfill(4)
        names.append(xsurface_name)

    if 'fb' in kwargs:
        xsurface_name = 'fb' + str(int(1000 * kwargs['fb'])).zfill(4)
        names.append(xsurface_name)

    if 'xsurface' in kwargs:
        xsurface_name = 'xs' + str(int(1000 * kwargs['xsurface'])).zfill(4)
        names.append(xsurface_name)

    if 'zval_lookup' in kwargs:
        xsurface_name = 'z' + str(int(10000 * kwargs['zval_lookup'])).zfill(4)
        names.append(xsurface_name)



    special_kwargs = ['mass', 'mcore', 'mslope', 'xcore', 'xsurface',
                      'alpha', 'zval_lookup', 'fa', 'fb', 'mabund',
                      'slope_ratio', 'dx_core', 'change_mass']

    other_kwargs = [key for key in kwargs if key not in special_kwargs]

    other_names = ["{:04.0f}".format(1000*kwargs[key]) for key in sorted(other_kwargs)]

    names = names + other_names
    new_name = '_'.join(names) + '.mod'
    new_name = run_dir + '/' + new_name
    
    return new_name



def generate_snapshot_model(starting_model, run_dir, build_fun, yscale='linear',
                            keep_original=False, i_original=0, change_mass_flag=False,
                            **kwargs):

    # print("Creating", kwargs)
    
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)

    new_name = generate_snapshot_name(run_dir, kwargs)

    # Create new model file and save to new_name
    build_new_h_gradient(starting_model, new_name, build_fun, **kwargs)

    # Modify the mass in the new model file if necessary
    if change_mass_flag:
        # print(kwargs['mass'])
        change_mass(new_name, new_name, kwargs['mass'])

    # Make a plot of the changes
    profileplot(starting_model, run_dir, new_name, yscale=yscale, **kwargs)

    # Copy the original model if desired
    if keep_original:
        subprocess.call(["cp", starting_model, run_dir + '/original.mod'])
        subprocess.call(["touch", run_dir + '/.original'])
        with open(run_dir + '/.original', 'w') as f:
            f.write(str(i_original))

    return new_name
    

def generate_snapshot_model_given_kwargs(lookup_table_dir, run_dir, build_fun, yscale='linear',
                                         keep_original=False, i_original=0, change_mass=False,
                                         **kwargs):
    """
    Uses generate_snapshot_model to generate single snapshot models based on set of kwargs
    """

    starting_model = find_nearest_model(kwargs['mass'],
                                kwargs['zval_lookup'],
                                kwargs['xcore'],
                                lookup_table_dir)

    model_name = generate_snapshot_model(starting_model, run_dir, build_fun, yscale=yscale,
                            keep_original=keep_original, i_original=i_original,
                            change_mass=change_mass, **kwargs)

    print("creating", model_name)

    return model_name



def generate_multiple_snapshot_models(lookup_table_dir, run_dir, build_fun, yscale='linear',
                                      keep_original=False, i_original=0, change_mass=False,
                                      **kwargs):

    """
    Uses generate_snapshot_model to generate many snapshot models based on list of kwargs
    """

    # Combinations of modifications to make
    keys = kwargs.keys()
    values = kwargs.values()
    combinations = [dict(zip(kwargs.keys(), combination)) for combination
                    in itertools.product(*kwargs.values())]
    
    model_names = []

    for combo in combinations:

        starting_model = find_nearest_model(combo['mass'],
                                    combo['zval_lookup'],
                                    combo['xcore'],
                                    lookup_table_dir)

        model_name = generate_snapshot_model(starting_model, run_dir, build_fun, yscale=yscale,
                                keep_original=keep_original, i_original=i_original,
                                change_mass=change_mass, **combo)
        model_names.append(model_name)
        print("creating", model_name)

    return model_names




def get_all_moddirs(starting_grid):
    """
    Needs to take starting models and label by M, Z and Xcore and 
    put into a single folder
    Also need to generate a file for quick searching of nearest model
    """

    
    masses = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    zvals = [0.003, 0.01, 0.015, 0.020]

    masses = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5]
    zvals = [0.003, 0.009, 0.014, 0.020]

    def get_name(mass, z):
        name = 'M' + str(int(10*mass)).zfill(3) + '_Z' + str(int(1000*z)).zfill(3)
        return name

    main_moddir = '/Users/eoin/Documents/Snapshot_Seismic/models'
    moddirs = [main_moddir + '/' + starting_grid + '/' + get_name(m, z)
               for m, z in itertools.product(masses, zvals)]

    return moddirs


    

def process_single_model(moddir, database_dir):
    """
    Post-process evolutionary model by copying .mod files
    to 
    """

    models = sorted(glob.glob(moddir + '/model*.mod'))

    for model in models:
        mod = MR(model, file_type='model')
        xcore = mod.h1[-1]

        xcore = str(int(10000*xcore)).zfill(4)
        
        mod_prefix = model.split('/')[-2]

        mod_name = mod_prefix + '_X' + xcore + '.mod'

        new_modname = database_dir + '/' + mod_name

        subprocess.call(["cp", model, new_modname])
        print("Copied", new_modname)




def make_lookup_table(database_dir, output_table_name='lookup_table',
                      overwrite=False):
    """
    Make lookup table to find nearest starting model
    """
    output_dir = database_dir + '/' + output_table_name

    if os.path.isfile(output_dir) and not overwrite:
        print("Make Lookup Table exited without writing a new file.\nFile exists or overwrite is set to False")
        return None


    models = sorted(glob.glob(database_dir + '/*.mod'))

    table = {}

    for model in models:

        table[model] = {}

        mod = MR(model, file_type='model')

        xcore = mod.h1[-1]
        mass = getattr(mod, 'M/Msun')
        zval = getattr(mod, 'initial_z')

        table[model]['xcore'] = xcore
        table[model]['mass'] = mass
        table[model]['zval'] = zval

    df = pd.DataFrame(table).T

    df['moddir'] = df.index.values

    print(df)

    if not os.path.isfile(output_dir) or overwrite:
        df.to_csv(output_dir)



def processs_model_grid(starting_grid):
    """
    process grid after it's computed
    """
    database_dir =  "/Users/eoin/Documents/Snapshot_Seismic/models/snapshot_database"

    moddirs = get_all_moddirs(starting_grid)
    for moddir in moddirs:
        process_single_model(moddir, database_dir)

    make_lookup_table(database_dir, output_table_name='lookup_table')

    print('---------- DONE ----------')



def find_nearest_model(mass, zval, xcore, lookup_table_dir):
    """
    Nearest model given by lookup table
    """
    df = pd.read_csv(lookup_table_dir)

    # Go in order of mass, then Z, then xcore

    mval_nearest = find_nearest(df['mass'], mass)
    df = df[df['mass'] == mval_nearest]

    zval_nearest = find_nearest(df['zval'], zval)
    df = df[df['zval'] == zval_nearest]

    xcore_nearest = find_nearest(df['xcore'], xcore)
    df = df[df['xcore'] == xcore_nearest]

    nearest_moddir = df['moddir'].values[0]

    return nearest_moddir



def compute_parameters_const_hemass(mcore, mslope, mass, xcore, xsurface, mslope_factors=[0.8, 1.2]):
    """
    Given a certain set of parameters describing a MS star, provide other parameters
    that describe another star with the same total H/He mass
    """

    hcore_mass = xcore * mcore * mass
    
    hsurface_mass = xsurface * (1 - mcore - mslope) * mass

    hslope_mass = xcore * mslope * mass + 0.5 * (xsurface - xcore) * (mslope * mass)

    h_mass = hcore_mass + hsurface_mass + hslope_mass

    # mslope_factors = [0.8, 1.2]

    mcores = []
    mslopes = []
    xcores = []
    xsurfaces = []

    for mslope_factor in mslope_factors:

        mslope_new = mslope * mslope_factor

        mcore_new = mcore + 0.5 * (mslope - mslope_new)

        xsurface_new = xsurface

        # print(mslope_new, mcore_new, xsurface_new)

        xcore_new = xsurface_new + (h_mass - xsurface_new*mass)/((mass) * (mcore_new + 0.5*mslope_new))

        mcores.append(mcore_new)
        mslopes.append(mslope_new)
        xcores.append(xcore_new)
        xsurfaces.append(xsurface_new)


    return mcores, mslopes, xcores, xsurfaces


# def fit_profile_with_params(modname):
    
#     from scipy.optimize import minimize
#     from mesa_reader import MesaData as MR

#     mod = MR(modname, file_type='model')
#     h1_vals = {'h1': mod.h1}
#     starmass = getattr(mod, 'M/Msun')
#     dq = mod.dq
#     mass = starmass * (np.cumsum(dq[::-1]))[::-1]
#     dm = starmass * dq
#     f_type = 'mod_ms_7param'

#     def objective(params, x, y):
#         # for scipy
#         param_names = ['mcore', 'mslope', 'xcore', 'xsurface', 'fa', 'fb', 'alpha']
#         kwargs = dict(zip(param_names, params))
#         new_h1 = functions[f_type](h1_vals, ['h1'], mass, dm, starmass, **kwargs)

#         return np.sum((new_h1 - h1_vals['h1']) ** 2)
#         # return new_h1

#     xcore_guess = mod.h1[-1]
#     xsurface_guess = mod.h1[0]

#     mcore_guess = next((m for m, h1 in zip(mass, mod.h1) if h1 < xcore_guess + 0.01))/starmass
#     mslope_guess = next((m for m, h1 in zip(mass, mod.h1) if h1 < xsurface_guess - 0.01)) - mcore_guess
#     mslope_guess = 0.2*mslope_guess/starmass

#     params0 = [mcore_guess, mslope_guess, xcore_guess, xsurface_guess, 0.5, 0.5, 0.00]
#     bounds = [(mcore_guess, mcore_guess),
#               (mslope_guess * 0.5, mslope_guess * 1.5),
#               (xcore_guess, xcore_guess),
#               (xsurface_guess, xsurface_guess),
#               (0, 1), (0, 1), (0, 0.01)]


#     # print(params0)
#     # print(bounds)
#     # print(mass)
#     # print(mod.h1)
#     # print()
#     result = minimize(objective, params0, args=(mass, mod.h1), bounds=bounds)
#     # result = minimize(get_profile, params0, args=(mass, mod.h1))

#     return params0, result.x




    return kwargs





if __name__ == '__main__':
    database_dir =  "/Users/eoin/Documents/Snapshot_Seismic/models/snapshot_database"
    lookup_table_dir = database_dir + '/lookup_table'
    work_dir = "/Users/eoin/Documents/Snapshot_Seismic/models/inversion_space"

    output_dir = 'test_1'
    run_dir =  work_dir + '/' + output_dir

    build_fun = 'mod_ms_7param'

    kwargs = {
    'mcore' : 0.15,
    'mslope' : 0.10,
    'mass' : 4.0,
    'xcore' : 0.20,
    'xsurface' : 0.70,
    'zval_lookup' : 0.020,
    'fa': 0.8,
    'fb': 0.2,
    'alpha': 0.002}

    starting_model = run_dir + '/M040_Z020_X2060_m0400_c0014_s0100_xc0200_xs0700_0002_0800_0200_0020.mod'
    
    generate_snapshot_model(starting_model, run_dir, build_fun, yscale='linear',
                                         keep_original=False, i_original=0, change_mass=True,
                                         **kwargs)

    # generate_snapshot_model_given_kwargs(lookup_table_dir, run_dir, build_fun, yscale='linear',
    #                                      keep_original=False, i_original=0, change_mass=True,
    #                                      **kwargs)
