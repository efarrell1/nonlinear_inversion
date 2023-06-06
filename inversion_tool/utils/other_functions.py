

def gyre_parametric_models(modlabels):
    moddir = '/Users/eoin/Documents/Snapshot_Seismic/models/parametric_nprofiles'
    moddirs = [moddir + '/' + modlabel for modlabel in modlabels]
    gyre_mods = [GyreMod(moddir) for moddir in moddirs]

    return gyre_mods


def get_static_models(subdir, modnames, rot=0):
    moddir = '/Users/eoin/Documents/Snapshot_Seismic/work/' + subdir + '/combined/LOGS'
    pfdirs = [moddir + '/' + x + '.data' for x in modnames]

    static_mods = [StellarStructure(pfdir, rot=rot) for pfdir in pfdirs]

    return static_mods


def all_gyre_models_from_subdir(subdir):

    moddir = '/Users/eoin/Documents/Snapshot_Seismic/work/' + subdir + '/combined/LOGS'
    gyre_dirs = sorted(glob(moddir + '/*/summary.h5'))
    gyre_dirs = ['/'.join(x.split('/')[:-1]) for x in gyre_dirs]
    print(len(gyre_dirs))

    gyre_mods = [GyreMod(moddir, read_mod=False) for moddir in gyre_dirs]
    
    # static_mods = [StellarStructure(pfdir, rot=None) for pfdir in pfdirs]

    # gyre_mods = [mod.get_all_gyre_mods() for mod in static_mods]

    # # Using a nested list comprehension
    # gyre_mods = [item for sublist in gyre_mods for item in sublist]

    return gyre_mods


def modnames_snapshot_subdir(subdir):
    directory = '/Users/eoin/Documents/Snapshot_Seismic/work/' + subdir
    file_type = '.mod'

    modnames = sorted([os.path.basename(f)[:-4] for f in os.listdir(directory) if f.endswith(file_type)])

    return modnames




def get_grid_models():

    with open(data_dir + '/grid_345_test_v3_modelnames.json', 'r') as f:
        # Load the list from the file using JSON
        model_names = json.load(f)

    name_to_params = {}
    for name in model_names:
        new_name = name.split('_m')[-1][:-4]
        new_params = new_name.split('_')
        param_vals = [re.sub("[^0-9]", "", input_string) for input_string in new_params]
        param_vals = [float(x)/1e3 for x in param_vals]
        params = ['mass', 'mcore', 'mslope', 'xc', 'xs', 'alpha', 'fa', 'fb', 'Z']
        dvals = dict(zip(params, param_vals))
        dvals['mass'] = dvals['mass'] * 10
        name_to_params[name] = dvals

    xc = 0.55
    alpha = 0.001
    mass = 3

    names = [name for name, val in name_to_params.items() if val['mass'] == mass and val['alpha'] == alpha]
    names = sorted(names)

    names = names[::2]

    names = [x.split('/')[-1][:-4] for x in names]

    return names



def create_hdf5_from_mods(static_mods, hdf5_filename, overwrite=False):

    hdf5_filename = '/Users/eoin/Documents/Snapshot_Seismic/models' + '/summaries/' + hdf5_filename

    if not overwrite and os.path.exists(hdf5_filename):
        print('File exists', hdf5_filename)
        return

    df_list = []

    for mod in static_mods:
        vals = np.column_stack((mod.nmode, mod.p, mod.dp))
        columns = ['nmode', 'p', 'dp']
        df = pd.DataFrame(vals, columns=columns)

        new_name = mod.moddir.split('_m')[-1][:-4]
        new_params = new_name.split('_')[:-1]
        param_vals = [re.sub("[^0-9]", "", input_string) for input_string in new_params]
        param_vals = [float(x)/1e3 for x in param_vals]
        params = ['mass', 'mcore', 'mslope', 'xc', 'xs', 'alpha', 'fa', 'fb', 'Z']

        for param, val in zip(params, param_vals):
            df[param] = val

        df['mass'] = 10 * df['mass']
        df['mcore'] = 10 * df['mcore']
        df['mslope'] = 10 * df['mslope']

        modname = mod.moddir.split('/')[-1]
        df['modname'] = modname
        df['subdir'] = mod.moddir.split('/')[-4]

        df['rot'] = mod.rot
        df['dp_rot'] = mod.dp_rot

        df_list.append(df)


    df = vaex.from_pandas(df=pd.concat(df_list))

    df.export_hdf5(hdf5_filename, mode='w')
        
    print('done')





def fquartic(x, a, b, c, d, e):
    return e * x ** 4 + d * x ** 3 + a * x ** 2 + b * x + c

def flinear(x, b, c):
    return  b * x + c

def fquadratic(x, a, b, c):
    return  a*x**2 + b * x + c

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






def get_params_given_period_spacing2(ps, nmin, nmax):

    slope, intercept = ps.fit_slope()

    params = {'slope': round(float(slope), 0),
              'intercept': round(float(intercept), 0),
              'P0': round(float(ps.p[0]), 5),
              'P1': round(float(ps.p[1]), 5),
              'P2': round(float(ps.p[2]), 5),
              'P3': round(float(ps.p[3]), 5),
              'P4': round(float(ps.p[4]), 5),
              'P5': round(float(ps.p[5]), 5),
              'P6': round(float(ps.p[6]), 5),
              'P6': round(float(ps.p[6]), 5),
              'P7': round(float(ps.p[7]), 5),
              'P8': round(float(ps.p[8]), 5),
              'P9': round(float(ps.p[9]), 5),
              'P10': round(float(ps.p[10]), 5),
              'P11': round(float(ps.p[11]), 5),
              'Plast': round(float(ps.p[-1]), 5),
              'dP0': round(float(ps.dp[0]), 5),
              'dP1': round(float(ps.dp[1]), 5),
              'dP2': round(float(ps.dp[2]), 5),
              'dP3': round(float(ps.dp[3]), 5),
              'dP4': round(float(ps.dp[4]), 5),
              'dP5': round(float(ps.dp[5]), 5),
              'dP6': round(float(ps.dp[6]), 5),
              'dP7': round(float(ps.dp[7]), 5),
              'dP8': round(float(ps.dp[8]), 5),
              'dP0_rot': round(float(ps.dp_rot[0]), 5),
              'dP1_rot': round(float(ps.dp_rot[1]), 5),
              'dP2_rot': round(float(ps.dp_rot[2]), 5),
              'dP3_rot': round(float(ps.dp_rot[3]), 5),
              'dP4_rot': round(float(ps.dp_rot[4]), 5),
              'dP5_rot': round(float(ps.dp_rot[5]), 5),
              'dP6_rot': round(float(ps.dp_rot[6]), 5),
              'dP7_rot': round(float(ps.dp_rot[7]), 5),
              'dP8_rot': round(float(ps.dp_rot[8]), 5)}

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

    return params


