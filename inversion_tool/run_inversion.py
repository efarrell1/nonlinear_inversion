#!/Users/eoin/opt/anaconda3/bin/python

import numpy as np

from run_build_star import generate_snapshot_model, find_nearest_model

import sys

import yaml



def save_parameters_to_yaml(parameters, file_path):
    with open(file_path, 'w') as file:
        yaml.dump(parameters, file)

def load_parameters_from_yaml(file_path):
    with open(file_path, 'r') as file:
        parameters = yaml.safe_load(file)
    return parameters


output_dir = sys.argv[1]

database_dir =  "/Users/eoin/Documents/Snapshot_Seismic/models/snapshot_database"
lookup_table_dir = database_dir + '/lookup_table'
work_dir = "/Users/eoin/Documents/Snapshot_Seismic/models/inversion_space"


run_dir =  work_dir + '/' + output_dir



# Load Log
try:
    log = load_parameters_from_yaml(run_dir + '/log.yaml')
except OSError:
    raise OSError("Could not find yaml file. Tried", run_dir + '/log.yaml')


build_fun = 'mod_ms_7param'
model_kwarg_names = ['mcore', 'mslope', 'new_mass', 'xcore', 'xsurface', 'zval_lookup', 'fa', 'fb', 'alpha']

model_kwargs = {x: log[x] for x in model_kwarg_names}


if 'starting_model' in log:
    starting_model = log['starting_model']
else:
    starting_model = find_nearest_model(log['new_mass'], log['zval_lookup'], log['xcore'], lookup_table_dir)


# Now iterate a bunch of times

generate_snapshot_model(starting_model, run_dir, log['build_fun'], change_mass=True, **model_kwargs)


