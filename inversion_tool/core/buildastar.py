import numpy as np

from .mesa_reader import MesaData as MR
from scipy import interpolate

import os
import codecs
import glob
import itertools
from pathlib import Path
import matplotlib.pyplot as plt
import subprocess
import json
import fortranformat as ff

#  ------ takes in .mod files and modifies the abundance structure ------- #
# Add new functions / methods to modify stellar abundance structure as
# functions like mod_grad_const_hemass and add to functions dictionary
# at bottom of page



def run_buildastar(moddir, build_fun, yscale='linear', keep_original=True, i_original=0, change_mass=False, **kwargs):
    fmoddir = cfg.draw + '/' + moddir  # original folder
    new_moddir = fmoddir + '/modified'  # directory for new models
    Path(new_moddir).mkdir(parents=True, exist_ok=True)  # create new folder
    models = sorted(glob.glob(fmoddir + '/*.mod'))  # list of .mod files to modify


    # Combinations of modifications to make
    keys = kwargs.keys()
    values = kwargs.values()
    combinations = [dict(zip(kwargs.keys(), combination)) for combination
                    in itertools.product(*kwargs.values())]
    dict_save = []
    print(" --- Building models --- ")
    for mod, combo in itertools.product(models, combinations):
        print("Creating", combo)
        dict_save.append(combo)
        names = ["{:04.0f}".format(1000*combo[key]) for key in sorted(combo)]
        names.insert(0, mod[:-4])
        new_name = '_'.join(names) + '.mod'
        new_name = new_moddir + '/' + new_name.split(fmoddir)[-1]
        build_new_h_gradient(mod, new_name, build_fun, **combo)
        if change_mass:
            print('changing mass')
            change_coreratio(mod, new_name, combo['menv'])
        profileplot(mod, new_name, yscale=yscale, **kwargs)

    if keep_original:
        subprocess.call(["cp", models[0], new_moddir + '/original.mod'])
        subprocess.call(["touch", fmoddir + '/.original'])
        with open(fmoddir + '/.original', 'w') as f:
            f.write(str(i_original))

    subprocess.call(["touch", fmoddir + '/' + build_fun])
    with open(fmoddir + '/' + build_fun, 'w') as f:
        json.dump(dict_save, f)
    print("\n")
    print(" --- Creating Inlists --- ")
    print(models[0])
    if 'arb1' in MR(models[0], file_type='model').bulk_names:
        subprocess.call([cfg.newmod_script, moddir, 'arb'])
    else:
        subprocess.call([cfg.newmod_script, moddir])


def profileplot(modname, run_dir, newmodname, yscale='linear', **kwargs):
    oldmodel = MR(modname, file_type='model')
    newmodel = MR(newmodname, file_type='model')

    # does the plot
    fig, (ax) = plt.subplots(1)

    starmass = getattr(oldmodel, 'M/Msun')
    # mass = starmass * (1 - np.cumsum(oldmodel.dq))
    mass = starmass * (np.cumsum(oldmodel.dq[::-1]))[::-1]

    oldh1 = oldmodel.h1
    newh1 = newmodel.h1

    h1, = ax.plot(mass, oldh1, label='old H1', linestyle=':')
    he4, = ax.plot(mass, oldmodel.he4, label='old He4', linestyle=':')
    c12, = ax.plot(mass, oldmodel.c12, label='old C12', linestyle=':')
    o16, = ax.plot(mass, oldmodel.o16, label='old O16', linestyle=':')
    # if arb:
    #     arb, = ax.plot(mass, getattr(oldmodel, arbn), label='old ' + arbn, linestyle=':')

    starmass = getattr(newmodel, 'M/Msun')
    # mass = starmass * (1 - np.cumsum(newmodel.dq))
    mass = starmass * (np.cumsum(newmodel.dq[::-1]))[::-1]
    
    ax.plot(mass, newmodel.h1, label='new H1', linestyle='-',
            color=h1.get_color())
    ax.plot(mass, newmodel.he4, label='new He4', linestyle='-',
            color=he4.get_color())
    ax.plot(mass, newmodel.c12, label='new C12', linestyle='-',
            color=c12.get_color())
    ax.plot(mass, newmodel.o16, label='new O16', linestyle='-',
            color=o16.get_color())

    if 'arb1' in newmodel.bulk_names:
        ax.plot(mass, newmodel.arb1, label='new arb1', linestyle='-')

    ax.legend(ncol=2, fontsize=4)
    ax.set_xlabel('Mass')
    ax.set_ylabel('H1')
    ax.set_yscale(yscale)

    image_dir = run_dir + "/images"
    Path(image_dir).mkdir(parents=True, exist_ok=True)
    fname = newmodname.split('/')[-1]
    fname = fname.split('.mod')[0]
    figdir = image_dir + '/' + fname

    fig.subplots_adjust(left=0.16, right=0.99, bottom=0.14, top=0.98)

    fig.savefig(figdir, dpi=300)
    # plt.show()
    plt.close(fig)


def calc_modified_profiles(model, isotopes, f_type, **kwargs):
    """
    function to take in old list of starmass, dqs, h and he and return new list
    of h and he
    """
    starmass = getattr(model, 'M/Msun')
    dq = model.dq
    mass = starmass * (np.cumsum(dq[::-1]))[::-1]
    dm = starmass * dq

    iso_profiles = np.array([getattr(model, iso) for iso in isotopes])
    iso_profiles = {iso: profile for (iso, profile)
                    in zip(isotopes, iso_profiles)}

    kwargs = {key:value for key, value in kwargs.items() if key != 'mass'}
    newabunds = functions[f_type](iso_profiles, isotopes, mass, dm, starmass,
                                  **kwargs)
    return newabunds


def get_i_hcore(iso_profiles, mass):
    xmin = iso_profiles['h1'][-1] + 0.01
    ibegin = next((x[0] for x in enumerate(iso_profiles['h1'])
                  if x[1] < iso_profiles['h1'][0]-0.01), len(mass) - 1)
    iend = next((x[0] for x in enumerate(iso_profiles['h1'][ibegin:])
                if x[1] < xmin), len(mass) - 1)
    if iend < len(mass) - 1:
        iend += ibegin

    return ibegin, iend


def add_arb1_to_model(modname, newmodeldir):
    """
    modeldir: Directory for model
    newmodeldir: Output directory for model
    """
    model = MR(modname, file_type='model')
    tempname = '/Users/eoin/Desktop/.temp.txt'

    not_isos = ['zone', 'lnd', 'lnT', 'lnR', 'L', 'dq', 'conv_vel', 'neut', 'mlt_vc']
    isotopes = [x for x in model.bulk_names if x not in not_isos]
    
    if 'arb1' in isotopes:
        print('Arb1 already in .mod')
        return
    
    # saves header and footer data
    inputFile = open(modname, "r+")
    inFilemod = inputFile.readlines()
    inputFile.close()
    
    nhead = 1 + next((x[0] for x in enumerate(inFilemod)
                     if x[1].strip()[:3] == 'lnd'))
    nfoot = 2 + next((x[0] for x in enumerate(inFilemod[::-1])
                     if x[1].strip()[:14] == 'previous model'))
 
    header = inFilemod[:nhead]
    footer = inFilemod[-nfoot:]
    # nfmt = len(not_isos) + len(isotopes) - 2 + 1
    # fmt = ('%5i    ' + nfmt*'%0.17E     ' + '%0.17E')

    # gets new abundances
    # newabundances = calc_modified_profiles(model, isotopes, f_type, **kwargs)
    # iso_profiles = np.array([getattr(model, iso) for iso in isotopes])
    # print(np.shape(iso_profiles))
 
    # rewrites file to python exponentials
    # inFile = codecs.open(modname, "r", "utf-8")
    # outFile = codecs.open(tempname, "w", "utf-8")
    # for line in inFile:
    #     newline = line.replace('D', 'e')
    #     outFile.write(newline)
    # inFile.close()
    # outFile.close()

    # modifies numpy array containing data
    moddata = np.genfromtxt(tempname, skip_header=nhead, skip_footer=5)
    # print(np.shape(moddata))
    iarb = model.bulk_names.index('he4') + 1
    moddata = np.insert(moddata, iarb, np.zeros(np.shape(moddata)[0]), axis=1)
    # print(np.shape(moddata))

    ihe4 = header[-1].find('he4') + 3
    header[-1] = header[-1][:ihe4] + ' '*23 + 'arb1' + header[-1][ihe4:]
    header = [x.replace('mesa_49.net', 'custom_arb.net') for x in header]
    header = [x.replace('species                              49', 'species                              50') for x in header]

    # gets new abundances
    # newabundances = calc_modified_profiles(model, isotopes, f_type, **kwargs)
    # iso_profiles = np.array([getattr(model, iso) for iso in isotopes])
    # moddata = np.genfromtxt(modname, skip_header=nhead, skip_footer=5)

    # Now save in MESA .mod format
    ff.config.RECORD_SEPARATOR = ''
    iformat = ff.FortranRecordWriter('1I5')
    lineformat = ff.FortranRecordWriter('(1E27.17)')

    inputFilenew = open(newmodeldir, "w")
    inputFilenew.writelines(header)
    for row in moddata:
        i = iformat.write([row[0]])
        n = lineformat.write(row[1:])
        inputFilenew.write(i + n + '\n')
    inputFilenew.writelines(footer)
    inputFilenew.close()


def read_inout(modname, newmodeldir):
    """
    modeldir: Directory for model
    newmodeldir: Output directory for model
    """
    model = MR(modname, file_type='model')
    tempname = '/Users/eoin/Desktop/.temp.txt'

    not_isos = ['zone', 'lnd', 'lnT', 'lnR', 'L', 'dq', 'conv_vel', 'neut']
    isotopes = [x for x in model.bulk_names if x not in not_isos]
    
    # get header and footer data
    with open(modname, "r+") as f:
        inFilemod = f.readlines()
    
    nhead = 1 + next((x[0] for x in enumerate(inFilemod)
                     if x[1].strip()[:3] == 'lnd'))
    nfoot = 2 + next((x[0] for x in enumerate(inFilemod[::-1])
                     if x[1].strip()[:14] == 'previous model'))
 
    header = inFilemod[:nhead]
    footer = inFilemod[-nfoot:]

    # gets new abundances
    # newabundances = calc_modified_profiles(model, isotopes, f_type, **kwargs)
    iso_profiles = np.array([getattr(model, iso) for iso in isotopes])
    moddata = np.genfromtxt(modname, skip_header=nhead, skip_footer=5)

    # Now save in MESA .mod format
    ff.config.RECORD_SEPARATOR = ''
    iformat = ff.FortranRecordWriter('1I5')
    lineformat = ff.FortranRecordWriter('(1E27.17)')

    inputFilenew = open(newmodeldir, "w")
    inputFilenew.writelines(header)
    for row in moddata:
        i = iformat.write([row[0]])
        n = lineformat.write(row[1:])
        inputFilenew.write(i + n + '\n')
    inputFilenew.writelines(footer)
    inputFilenew.close()



def change_mass(modname, newmodeldir, new_mass):
    """
    modeldir: Directory for model
    newmodeldir: Output directory for model
    """
    model = MR(modname, file_type='model')
    tempname = '/Users/eoin/Desktop/.temp.txt'

    not_isos = ['zone', 'lnd', 'lnT', 'lnR', 'L', 'dq', 'conv_vel', 'neut']
    isotopes = [x for x in model.bulk_names if x not in not_isos]
    
    # get header and footer data
    with open(modname, "r+") as f:
        inFilemod = f.readlines()
    
    nhead = 1 + next((x[0] for x in enumerate(inFilemod)
                     if x[1].strip()[:3] == 'lnd'))
    nfoot = 2 + next((x[0] for x in enumerate(inFilemod[::-1])
                     if x[1].strip()[:14] == 'previous model'))
 
    header = inFilemod[:nhead]
    footer = inFilemod[-nfoot:]

    # gets new abundances
    # newabundances = calc_modified_profiles(model, isotopes, f_type, **kwargs)
    iso_profiles = np.array([getattr(model, iso) for iso in isotopes])
    moddata = np.genfromtxt(modname, skip_header=nhead, skip_footer=5)

    # Now save in MESA .mod format
    ff.config.RECORD_SEPARATOR = ''
    iformat = ff.FortranRecordWriter('1I5')
    lineformat = ff.FortranRecordWriter('(1E27.17)')

    lineformat_mass = ff.FortranRecordWriter('(ES28.16)')
    imass = [x.strip()[:6] for x in header].index('M/Msun')
    nums = lineformat_mass.write([new_mass])
    newline = str(header[imass][:32] + nums + '\n')
    newline = newline.replace('E', 'D')
    header[imass] = newline
    
    inputFilenew = open(newmodeldir, "w")
    inputFilenew.writelines(header)
    for row in moddata:
        i = iformat.write([row[0]])
        n = lineformat.write(row[1:])
        inputFilenew.write(i + n + '\n')
    inputFilenew.writelines(footer)
    inputFilenew.close()



def change_coreratio(modname, newmodeldir, menv):
    """
    modeldir: Directory for model
    newmodeldir: Output directory for model
    """
    model = MR(modname, file_type='model')
    tempname = '/Users/eoin/Desktop/.temp.txt'

    not_isos = ['zone', 'lnd', 'lnT', 'lnR', 'L', 'dq', 'conv_vel', 'neut']
    isotopes = [x for x in model.bulk_names if x not in not_isos]
    iso_profiles = np.array([getattr(model, iso) for iso in isotopes])
    iso_profiles = {iso: profile for (iso, profile)
                    in zip(isotopes, iso_profiles)}
    
    starmass = getattr(model, 'M/Msun')
    dq = model.dq
    mass = starmass * (np.cumsum(dq))
        
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-4), len(mass) - 1)
    # print(iend, iso_profiles['h1'][iend], mass[iend], starmass)

    
    mcore = starmass - mass[iend]
    new_mass = mcore + menv
    # print(new_mass, mcore, menv, starmass, mass[iend])
    
    # get header and footer data
    with open(newmodeldir, "r+") as f:
        inFilemod = f.readlines()
    
    nhead = 1 + next((x[0] for x in enumerate(inFilemod)
                     if x[1].strip()[:3] == 'lnd'))
    nfoot = 2 + next((x[0] for x in enumerate(inFilemod[::-1])
                     if x[1].strip()[:14] == 'previous model'))
 
    header = inFilemod[:nhead]
    footer = inFilemod[-nfoot:]

    # gets new abundances
    # newabundances = calc_modified_profiles(model, isotopes, f_type, **kwargs)
    iso_profiles = np.array([getattr(model, iso) for iso in isotopes])
    moddata = np.genfromtxt(newmodeldir, skip_header=nhead, skip_footer=5)

    # Now save in MESA .mod format
    ff.config.RECORD_SEPARATOR = ''
    iformat = ff.FortranRecordWriter('1I5')
    lineformat = ff.FortranRecordWriter('(1E27.17)')

    lineformat_mass = ff.FortranRecordWriter('(ES28.16)')
    imass = [x.strip()[:6] for x in header].index('M/Msun')
    nums = lineformat_mass.write([new_mass])
    newline = str(header[imass][:32] + nums + '\n')
    newline = newline.replace('E', 'D')
    header[imass] = newline
    
    inputFilenew = open(newmodeldir, "w")
    inputFilenew.writelines(header)
    for row in moddata:
        i = iformat.write([row[0]])
        n = lineformat.write(row[1:])
        inputFilenew.write(i + n + '\n')
    inputFilenew.writelines(footer)
    inputFilenew.close()
    # print(inputFilenew)




def build_new_h_gradient(modname, newmodeldir, f_type, mesa_v=15140, arb=False, **kwargs):
    """
    modeldir: Directory for model
    newmodeldir: Output directory for model
    """
    model = MR(modname, file_type='model')
    not_isos = ['zone', 'lnd', 'lnT', 'lnR', 'L', 'dq', 'conv_vel', 'neut', 'mlt_vc']
    isotopes = [x for x in model.bulk_names if x not in not_isos]

    inputFile = open(modname, "r+")
    inFilemod = inputFile.readlines()
    inputFile.close()
    
    nhead = 1 + next((x[0] for x in enumerate(inFilemod)
                     if x[1].strip()[:3] == 'lnd'))
    nfoot = 2 + next((x[0] for x in enumerate(inFilemod[::-1])
                     if x[1].strip()[:14] == 'previous model'))
 
    header = inFilemod[:nhead]
    footer = inFilemod[-nfoot:]

    newabundances = calc_modified_profiles(model, isotopes, f_type, **kwargs)

    # modifies numpy array containing data
    
    if nhead <=200:
        # if using older MESA versions which save as 5.67D08 which
        # isn't understood by numpy. Newer versions use 5.67E08
        tempname = '/Users/eoin/Desktop/.temp.txt'
        # rewrites file to python exponentials
        inFile = codecs.open(modname, "r", "utf-8")
        outFile = codecs.open(tempname, "w", "utf-8")
        for line in inFile:
            newline = line.replace('D', 'e')
            outFile.write(newline)
        inFile.close()
        outFile.close()
        moddata = np.genfromtxt(tempname, skip_header=nhead, skip_footer=5)
    else:
        moddata = np.genfromtxt(modname, skip_header=nhead, skip_footer=5)
    

    for iso, abund in zip(isotopes, newabundances):
        index = model.bulk_names.index(iso)
        moddata[:, index] = abund

    # Now save in MESA .mod format
    ff.config.RECORD_SEPARATOR = ''
    iformat = ff.FortranRecordWriter('1I5')
    lineformat = ff.FortranRecordWriter('(1E27.17)')

    inputFilenew = open(newmodeldir, "w")
    inputFilenew.writelines(header)
    for row in moddata:
        i = iformat.write([row[0]])
        n = lineformat.write(row[1:])
        inputFilenew.write(i + n + '\n')
    inputFilenew.writelines(footer)
    inputFilenew.close()
    


def mod_grad_const_hemass(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    # modifies gradient at constant total H and He Mass
    def get_ind_steep(itry_a, oldabund, indices=False):
        # returns list of total helium masses and corresponding indices
        f_abund = interpolate.interp1d(np.array(norm_old_mass).astype(float),
                                       np.array(oldabund[ibegin:iend]).astype(float),
                                       kind='cubic', fill_value='extrapolate')
        newendmass = mass[itry_a] + newdeltam
        itry_b = next((x[0] for x in enumerate(mass) if x[1] > newendmass), 0)
        norm_new_mass = (mass[itry_a:itry_b] - mass[itry_a]) / newdeltam
        norm_new_mass = np.array(norm_new_mass).astype(float)

        p1 = np.array((itry_a - ibegin) * [oldabund[0]])
        p2 = oldabund[:ibegin]
        p3 = f_abund(norm_new_mass)
        iextra = iend - itry_b
        if iextra > 0:
            p4 = oldabund[iend:]
            p5 = iextra * [oldabund[-1]]
        else:
            p4 = oldabund[itry_b:]
            p5 = []

        # psizes = [len(list(p)) for p in (p1, p2, p3, p4, p5)]
        newabund = np.concatenate((p1, p2, p3, p4, p5))
        print(len(newabund), len(dm))
        if indices:
            return [itry_a, itry_b, np.dot(dm, newabund)]
        else:
            return newabund

    def get_ind_shallow(itry_b, oldabund, indices=False):
        # returns list of total helium masses and corresponding indices

        f_abund = interpolate.interp1d(np.array(norm_old_mass).astype(float), np.array(oldabund[ibegin:iend]).astype(float),
                                       kind='cubic', fill_value='extrapolate')
        newendmass = mass[itry_b] - newdeltam
        itry_a = next((x[0] for x in enumerate(mass) if x[1] > newendmass), 0)
        norm_new_mass = (mass[itry_a:itry_b] - mass[itry_a]) / newdeltam
        norm_new_mass = np.array(norm_new_mass).astype(float)

        iextra = itry_a - ibegin
        if iextra > 0:
            p1 = iextra * [oldabund[0]]
            p2 = oldabund[:ibegin]
            p3 = f_abund(norm_new_mass)
            p4 = oldabund[itry_b:]
        else:
            p1 = oldabund[ibegin - itry_a:ibegin]
            p2 = f_abund(norm_new_mass)
            p3 = oldabund[itry_b:]
            p4 = []

        newabund = np.concatenate((p1, p2, p3, p4))
        if indices:
            return [itry_a, itry_b, np.dot(dm, newabund)]
        else:
            return newabund

    def get_abunds():
        if slopefactor >= 1:
            ias, ibs, h1_masses = np.array([get_ind_steep(itry, old_h1,
                                            indices=True) for itry in
                                            range(ibegin, iend)]).T
            ia, ib, m = next(((int(ia), int(ib), m)
                             for ia, ib, m in zip(ias, ibs, h1_masses)
                             if m >= total_h1_mass), (ibegin, iend,
                             total_h1_mass))
            # print(m)
            newabunds = [get_ind_steep(ia, iso_profiles[iso], indices=False)
                         for iso in isotopes]
        elif slopefactor < 1:
            ias, ibs, h1_masses = np.array([get_ind_shallow(itry, old_h1,
                                            indices=True) for itry in
                                            range(iend, len(old_h1))]).T
            ia, ib, m = next(((int(ia), int(ib), m) for ia, ib, m in zip(ias,
                             ibs, h1_masses) if m >= total_h1_mass),
                             (ibegin, iend, total_h1_mass))
            # print(m)
            newabunds = [get_ind_shallow(ib, iso_profiles[iso], indices=False)
                         for iso in isotopes]

        return newabunds

    ibegin, iend = get_i_hcore(iso_profiles, mass)
    slopefactor = kwargs['slopefactor']
    total_he4_mass = np.dot(dm, iso_profiles['he4'])
    total_h1_mass = np.dot(dm, iso_profiles['h1'])
    old_h1 = iso_profiles['h1']
    old_he4 = iso_profiles['he4']

    # original slope
    deltam_gradient = mass[iend] - mass[ibegin]
    deltah_gradient = iso_profiles['h1'][iend] - iso_profiles['h1'][ibegin]
    oldslope = -1 * deltah_gradient / deltam_gradient

    # figure out new slope
    newslope = oldslope * slopefactor
    newdeltam = -1 * deltah_gradient / newslope

    # normalised mass in gradient region
    norm_old_mass = (mass[ibegin:iend] - mass[ibegin])/deltam_gradient
    newabunds = get_abunds()
    return newabunds


def mod_h_grad(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes just the H gradient above the H depleted core
    (without maintaining constant He core mass)
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    slopefactor = kwargs['slopefactor']
    # mass of the core
    deltam_core = starmass - mass[iend]

    # new slope # change in mass and h over the gradient region
    deltam_gradient = mass[iend] - mass[ibegin]
    deltah_gradient = iso_profiles['h1'][iend] - iso_profiles['h1'][ibegin]

    # original slope
    slope = -1 * deltah_gradient / deltam_gradient
    # new slope
    newslope = slope * slopefactor
    # new change in mass of gradient region
    newdeltam = -1 * deltah_gradient / newslope

    hmass_grad = np.dot(dm[ibegin:iend], iso_profiles['h1'][ibegin:iend])
    he4mass_grad = np.dot(dm[ibegin:iend], iso_profiles['he4'][ibegin:iend])

    # original slope
    slope = -1 * deltah_gradient / deltam_gradient
    # new slope
    newslope = slope * slopefactor
    # new change in mass of gradient region
    newdeltam = -1 * deltah_gradient / newslope
    # print(slope, newslope, deltam_gradient, newdeltam)

    # not important
    changeh = (1 - slopefactor) * slope
    # new mass of top of gradient
    newendmass = mass[iend] - newdeltam

    # new mass coordinate at top of gradient
    inewmass = next((x[0] for x in enumerate(mass) if x[1] > newendmass), 0)
    # print(ibegin, inewmass)

    # normalised mass in gradient region
    norm_old_mass = (mass[ibegin:iend] - mass[ibegin])/deltam_gradient
    norm_new_mass = (mass[inewmass:iend] - mass[inewmass]) / \
    (mass[iend] - mass[inewmass])

    def get_abund_isotope(ibegin, iend, inewmass, norm_old_mass, norm_new_mass,
                          oldabund):
        """
        retusn isotopes abundance for new structure
        """
        # creates interpolation for normalised mass coordiante in gradient
        f_abund = interpolate.interp1d(norm_old_mass, oldabund[ibegin:iend],
                                       kind='cubic', fill_value='extrapolate')

        if ibegin > inewmass:
            newabund = np.concatenate((oldabund[ibegin-inewmass:ibegin],
                                      f_abund(norm_new_mass), oldabund[iend:]))
        else:
            extrai = inewmass - ibegin
            newabund = np.concatenate((extrai * [oldabund[0]],
                                      oldabund[:ibegin],
                                      f_abund(norm_new_mass), oldabund[iend:]))
        return newabund

    newabunds = [get_abund_isotope(ibegin, iend, inewmass, norm_old_mass,
                 norm_new_mass, iso_profiles[iso]) for iso in isotopes]

    return newabunds


def mod_h_grad_absolute(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient as specified
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # grad_frac is the fraction of Menv that consists of a Gradient
    grad_frac = kwargs['grad_frac']

    z = 0.02
    try:
        he_surface = kwargs['he_surface']
        h_surf = 1 - z - he_surface
    except KeyError:
        # From Choi+2016, p.3
        # z = 0.02
        yprotosol = 0.2703
        zprotosol = 0.0142
        yprim = 0.249  # y primordial
        he_surface = yprim + (yprotosol - yprim) * float(z)/(zprotosol)
        h_surf = 1 - z - he_surface

    # mass of the envelope
    menv = mass[iend]
    mflat = menv * (1 - grad_frac)
    ib_new = next((x[0] for x in enumerate(mass)
                      if x[1] > mflat), len(mass) - 1)
    masses_grad = mass[ib_new: iend + 1]

    h1_grad = [h_surf + (m - mflat)/(menv - mflat) * (iso_profiles['h1'][iend] - h_surf)
               for m in masses_grad]

    h1_flat = [h_surf for i in range(0, ib_new)]
    h1_core = iso_profiles['h1'][iend + 1:]

    # print(h1_core)
    # print(iso_profiles['h1'][iend], iso_profiles['h1'][iend + 1])
    # print(h1_flat[0], h1_flat[-1], h1_grad[0], h1_grad[-1], h1_core[0])

    h1_prof = np.concatenate((h1_flat, h1_grad, h1_core))

    metals = [x for x in isotopes if x not in ['h1', 'he4']]

    z_prof = np.array([iso_profiles[x] for x in metals])
    zsum = np.sum(z_prof, axis=0)
    # zenv = zsum[:iend + 1]
    he4_prof = 1 - h1_prof - zsum

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h1_prof
        elif isotope == 'he4':
            return he4_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    # print([len(x) for x in newabunds])

    return newabunds


def mod_draw_grad(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient however we like
    kwargs:
        - fix: either 'h1' or 'he4', depending on which one you want to keep
               constant (as the total depends on metallicity)
        - abund_surf: surface abundance of fix
        - grad_frac: fraction of Menv that contains a gradient in H/He

    Only works for Core He burning models
    """
    # ibegin, iend = get_i_hcore(iso_profiles, mass)
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-4), len(mass) - 1)
    # grad_frac is the fraction of Menv that consists of a Gradient
    grad_frac = kwargs['grad_frac']

    try:
        fix = kwargs['fix']
    except KeyError:
        fix = 'h1'

    try:
        abund_surf = kwargs['abund_surf']
    except KeyError:
        # By default fix H1 and let He vary depending on metallicity
        # and set surface H/He value based on Choi+2016
        z = 0.02
        yprotosol = 0.2703
        zprotosol = 0.0142
        yprim = 0.249  # y primordial
        he_surface = yprim + (yprotosol - yprim) * float(z)/(zprotosol)
        if fix == 'he4':
            abund_surf = he_surface
        elif fix == 'h1':
            abund_surf = 1 - z - he_surface

    # mass of the envelope
    menv = mass[iend]
    mflat = menv * (1 - grad_frac)
    ib_new = next((x[0] for x in enumerate(mass)
                  if x[1] > mflat), len(mass) - 1)
    masses_grad = mass[ib_new: iend + 1]

    abund_grad = [abund_surf + (m - mflat)/(menv - mflat) *
                  abs((iso_profiles[fix][iend] - abund_surf))
                  for m in masses_grad]

    abund_flat = [abund_surf for i in range(0, ib_new)]
    abund_core = iso_profiles[fix][iend + 1:]

    abund_prof = np.concatenate((abund_flat, abund_grad, abund_core))

    metals = [x for x in isotopes if x not in ['h1', 'he4']]

    z_prof = np.array([iso_profiles[x] for x in metals])
    zsum = np.sum(z_prof, axis=0)
    other_prof = 1 - abund_prof - zsum

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == fix:
            return abund_prof
        elif isotope in ['h1', 'he4']:
            return other_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_change_z(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient however we like
    kwargs:
        - fix: either 'h1' or 'he4', depending on which one you want to keep
               constant (as the total depends on metallicity)
        - abund_surf: surface abundance of fix

    Only works for Core He burning models
    """
    # This is the edge of he core
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-4), len(mass) - 1)

    # How to find edge of convective He core
    ihemax = np.argmax(iso_profiles['he4'])
    hemax = iso_profiles['he4'][ihemax]
    iconv = ihemax + next((x[0] for x in
                          enumerate(iso_profiles['he4'][ihemax:])
                          if x[1] < hemax - 0.01), len(mass) - 1 - ihemax)

    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']

    metals = [x for x in isotopes if x not in ['h1', 'he4']]
    z_prof = np.array([iso_profiles[x] for x in metals])
    zsum = np.sum(z_prof, axis=0)
    znew = zsum * z_fact

    # print(iend)
    h_prof_env = np.add(iso_profiles['h1'][:iend], znew[:iend])
    h_prof_core = iso_profiles['h1'][iend:]
    h_prof = np.concatenate((h_prof_env, h_prof_core))

    he_prof_env = iso_profiles['he4'][:iend]
    he_prof_shell = np.add(iso_profiles['he4'][iend:iconv], znew[iend:iconv])
    he_prof_core = iso_profiles['he4'][iconv:]
    he_prof = np.concatenate((he_prof_env, he_prof_shell, he_prof_core))

    # print(len(h_prof), len(he_prof))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'he4':
            return he_prof
        elif isotope == 'h1':
            return h_prof
        else:
            iso_conv = iso_profiles[isotope][iconv:]
            iso_rest = iso_profiles[isotope][:iconv] * z_fact
            return np.concatenate((iso_rest, iso_conv))

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_change_core_abund(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient however we like
    kwargs:
        - fix: either 'h1' or 'he4', depending on which one you want to keep
               constant (as the total depends on metallicity)
        - abund_surf: surface abundance of fix

    Only works for Core He burning models
    """
    # This is the edge of he core
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-4), len(mass) - 1)

    # # How to find edge of convective He core
    ihemax = np.argmax(iso_profiles['he4'])
    hemax = iso_profiles['he4'][ihemax]
    iconv = ihemax + next((x[0] for x in
                          enumerate(iso_profiles['he4'][ihemax:])
                          if x[1] < hemax - 0.01), len(mass) - 1 - ihemax)

    # print(len(mass), iconv, iso_profiles['he4'][iconv], iso_profiles['he4'][-1])

    # Alternative to finding edge of convective He core
    # he_reverse = iso_profiles['he4'][::-1]
    # he_center = he_reverse[0]
    # iconv_r = next((x[0] for x in enumerate(he_reverse)
    #                if x[1] > he_center + 0.01), 0)
    # iconv = len(mass) - iconv_r
    #
    # print(len(mass), iconv, iso_profiles['he4'][iconv], iso_profiles['he4'][-1])

    #  yc is central He4 abundance.
    yc = kwargs['yc']

    # moddir = '/Users/eoin/Documents/red_supergiants/models/SMC/M160'
    # folder = '/standard/LOGS/history.data'

    moddir = '/Users/eoin/Documents/snapshot_models/models/grid20/M160'
    folder = '/standard/LOGS/history.data'

    h = MR(moddir + folder)

    ib_he = next((x[0] for x in enumerate(h.center_h1)
                 if x[1] < 1e-4), len(h.center_h1) - 1)

    ie_he = ib_he + next((x[0] for x in enumerate(h.center_he4[ib_he:])
                         if x[1] < 1e-3), len(h.center_he4[ib_he:]) - 1)

    he4s = h.center_he4[ib_he:ie_he]
    c12s = h.center_c12[ib_he:ie_he]
    o16s = h.center_o16[ib_he:ie_he]

    # functions to find C12 and O16 values as a function of Yc
    fc12 = interpolate.interp1d(he4s, c12s, kind='linear',
                                fill_value='extrapolate')
    fo16 = interpolate.interp1d(he4s, o16s, kind='linear',
                                fill_value='extrapolate')

    he4_cen = iso_profiles['he4'][-1]
    he4_edge = iso_profiles['he4'][iconv]
    he4_f = (yc - he4_cen) / (he4_edge - he4_cen)
    he4_core = [x + he4_f * (he4_edge - x)
                for x in iso_profiles['he4'][iconv:]]
    he4_prof = np.concatenate((iso_profiles['he4'][:iconv], he4_core))

    c12_cen = iso_profiles['c12'][-1]
    c12_edge = iso_profiles['c12'][iconv]
    c12_f = (fc12(yc) - c12_cen) / (c12_edge - c12_cen)
    c12_core = [x + c12_f * (c12_edge - x)
                for x in iso_profiles['c12'][iconv:]]
    c12_prof = np.concatenate((iso_profiles['c12'][:iconv], c12_core))

    o16_cen = iso_profiles['o16'][-1]
    o16_edge = iso_profiles['o16'][iconv]
    o16_f = (fo16(yc) - o16_cen) / (o16_edge - o16_cen)
    o16_core = [x + o16_f * (o16_edge - x)
                for x in iso_profiles['o16'][iconv:]]
    o16_prof = np.concatenate((iso_profiles['o16'][:iconv], o16_core))

    profiles = {'he4': he4_prof, 'c12': c12_prof, 'o16': o16_prof}

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope in profiles:
            return profiles[isotope]
        else:
            return iso_profiles[isotope]

    # print(iso_profiles['he4'][-1], iso_profiles['c12'][-1], iso_profiles['o16'][-1])
    # print(he4_prof[-1], c12_prof[-1], o16_prof[-1])
    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_env_const_hemass(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    # modifies gradient at constant total H and He Mass
    # applies to core He burning stars as only works outside He core

    ibegin, iend = get_i_hcore(iso_profiles, mass)

    # grad_frac is the fraction of Menv that consists of a Gradient
    grad_frac = kwargs['grad_frac']
    # mass of the envelope
    menv = mass[iend]
    mflat = menv * (1 - grad_frac)
    ib_new = next((x[0] for x in enumerate(mass)
                   if x[1] > mflat), len(mass) - 1)
    masses_grad = mass[ib_new: iend + 1]
    masses_flat = mass[:ib_new]

    env_he4_mass = np.dot(dm[:iend], iso_profiles['he4'][:iend])
    env_h1_mass = np.dot(dm[:iend], iso_profiles['h1'][:iend])
    old_h1 = iso_profiles['h1']
    old_he4 = iso_profiles['he4']

    hsurf_trys = np.arange(0.05, 1, 0.005)

    def compute_hmass(h_surf):
        h1_grad = [h_surf + (m - mflat)/(menv - mflat) *
                   (iso_profiles['h1'][iend] - h_surf)
                   for m in masses_grad]
        total_grad = np.dot(dm[ib_new:iend + 1], h1_grad)

        h1_flat = [h_surf for i in range(0, ib_new)]
        total_flat = np.dot(dm[:ib_new], h1_flat)
        return total_grad + total_flat

    hmasses = np.array([compute_hmass(x) for x in hsurf_trys])
    hdiffs = np.abs(hmasses - env_h1_mass)
    # print(hdiffs)
    hsurf_optimal = hsurf_trys[np.argmin(hdiffs)]

    h1_grad = [hsurf_optimal + (m - mflat)/(menv - mflat) *
               (iso_profiles['h1'][iend] - hsurf_optimal)
               for m in masses_grad]
    h1_flat = [hsurf_optimal for i in range(0, ib_new)]

    h1_core = [x for x in old_h1[iend + 1:]]

    new_h1 = np.concatenate((h1_flat, h1_grad, h1_core))
    # print(len(new_h1), len(mass))
    # print(mass[0], mass[500], mass[iend])
    # print(env_h1_mass, np.dot(dm[:iend], new_h1[:iend]))

    h1_change = new_h1 - old_h1
    new_he4 = old_he4 - h1_change
    # print(len(new_he4))
    # print(env_he4_mass, np.dot(dm[:iend], new_he4[:iend]))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return new_h1
        elif isotope == 'he4':
            return new_he4
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_env_const_grad(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    # modifies H gradient for a given slope and He surface abundance.
    # applies to core He burning stars as only works outside He core

    ibegin, iend = get_i_hcore(iso_profiles, mass)

    # slope is the slope of a gradient
    slope = kwargs['slope']
    abund_surf = kwargs['h_surface']
    fix = 'h1'

    # mass of the envelope
    menv = mass[iend]
    mflat = menv - ((abund_surf - iso_profiles[fix][iend]) / slope)

    ib_new = next((x[0] for x in enumerate(mass)
                   if x[1] > mflat), len(mass) - 1)
    masses_grad = mass[ib_new: iend + 1]
    masses_flat = mass[:ib_new]

    abund_grad = [abund_surf + (m - mflat)/(menv - mflat) *
                  (iso_profiles[fix][iend] - abund_surf)
                  for m in masses_grad]

    abund_flat = [abund_surf for i in range(0, ib_new)]
    abund_core = iso_profiles[fix][iend + 1:]

    abund_prof = np.concatenate((abund_flat, abund_grad, abund_core))

    metals = [x for x in isotopes if x not in ['h1', 'he4']]

    z_prof = np.array([iso_profiles[x] for x in metals])
    zsum = np.sum(z_prof, axis=0)
    other_prof = 1 - abund_prof - zsum

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == fix:
            return abund_prof
        elif isotope in ['h1', 'he4']:
            return other_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_grad_absolute_flex_surface(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient as specified
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # grad_frac is the fraction of Menv that consists of a Gradient
    grad_frac = kwargs['grad_frac']
    h_grad = kwargs['h_grad']
    h_surf = kwargs['h_surface']
    h_grad = min(h_grad, h_surf)

    # mass of the envelope
    menv = mass[iend]
    mflat = menv * (1 - grad_frac)
    ib_new = next((x[0] for x in enumerate(mass)
                  if x[1] > mflat), len(mass) - 1)
    masses_grad = mass[ib_new: iend + 1]

    h1_grad = [h_grad + (m - mflat)/(menv - mflat) * (iso_profiles['h1'][iend] - h_grad)
               for m in masses_grad]

    h1_flat = [h_surf for i in range(0, ib_new)]
    h1_core = iso_profiles['h1'][iend + 1:]

    h1_prof = np.concatenate((h1_flat, h1_grad, h1_core))

    metals = [x for x in isotopes if x not in ['h1', 'he4']]

    z_prof = np.array([iso_profiles[x] for x in metals])
    zsum = np.sum(z_prof, axis=0)
    he4_prof = 1 - h1_prof - zsum

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h1_prof
        elif isotope == 'he4':
            return he4_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_change_z_outonly(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient however we like
    kwargs:
        - fix: either 'h1' or 'he4', depending on which one you want to keep
               constant (as the total depends on metallicity)
        - abund_surf: surface abundance of fix

    Only works for Core He burning models
    """
    # This is the edge of he core
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-4), len(mass) - 1)

    # How to find edge of convective He core
    ihemax = np.argmax(iso_profiles['he4'])
    hemax = iso_profiles['he4'][ihemax]
    iconv = ihemax + next((x[0] for x in
                          enumerate(iso_profiles['he4'][ihemax:])
                          if x[1] < hemax - 0.01), len(mass) - 1 - ihemax)

    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']
    env_frac = kwargs['env_frac']

    menv = mass[iend]
    menv_change = menv * env_frac  # modify this amount of envelope

    ichange = next((x[0] for x in enumerate(mass) if x[1] > menv_change), 0)

    metals = [x for x in isotopes if x not in ['h1', 'he4']]
    z_prof = np.array([iso_profiles[x] for x in metals])
    zsum = np.sum(z_prof, axis=0)
    znew = zsum * z_fact

    # print(iend)
    h_prof_env = np.add(iso_profiles['h1'][:ichange], znew[:ichange])
    h_prof_rest = iso_profiles['h1'][ichange:]
    h_prof = np.concatenate((h_prof_env, h_prof_rest))

    he_prof = iso_profiles['he4']

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'he4':
            return he_prof
        elif isotope == 'h1':
            return h_prof
        else:
            iso_new = iso_profiles[isotope][:ichange] * z_fact
            iso_rest = iso_profiles[isotope][ichange:]
            return np.concatenate((iso_new, iso_rest))

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_ms_cno_core(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes the metallicity in a MS star for a given fraction of the star
    beginning at the surface
    kwargs:
        - fix: either 'h1' or 'he4', depending on which one you want to keep
               constant (as the total depends on metallicity)
        - abund_surf: surface abundance of fix

    Only works for MS models
    """
    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']
    # frac is the inner fraction of the envelope to modify
    frac = kwargs['frac']

    menv_change = frac * mass[-1]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > mass[-1] - menv_change), 0)
    cno = ['c12', 'c13', 'n13', 'n14', 'n15', 'o14', 'o15', 'o16', 'o17',
           'o18']

    cno_prof = np.array([iso_profiles[x] for x in cno])
    cnosum = np.sum(cno_prof, axis=0)
    cnonew = cnosum * z_fact
    diff = cnosum - cnonew

    h_prof_out = iso_profiles['h1'][:ichange]
    h_prof_mod = np.add(iso_profiles['h1'][ichange:], diff[ichange:])
    h_prof = np.concatenate((h_prof_out, h_prof_mod))

    he_prof = iso_profiles['he4']

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h_prof
        elif isotope not in cno:
            return iso_profiles[isotope]
        else:
            iso_env = iso_profiles[isotope][:ichange]
            iso_mod = iso_profiles[isotope][ichange:] * z_fact
            return np.concatenate((iso_env, iso_mod))
    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds




def mod_flexsurf_arb(iso_profiles, isotopes, mass, dm, starmass, arbn='arb1',
                     **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient as specified
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # grad_frac is the fraction of Menv that consists of a Gradient
    grad_frac = kwargs['grad_frac']
    # h_grad is fraction of H at outer edge of gradient region
    h_grad = kwargs['h_grad']
    h_surf = kwargs['h_surface']
    arb_surf = kwargs['arb_surf']
    h_grad = min(h_grad, h_surf)

    arb_v = max(h_surf, h_grad)

    # mass of the envelope
    menv = mass[iend]
    mflat = menv * (1 - grad_frac)
    ib_new = next((x[0] for x in enumerate(mass)
                  if x[1] > mflat), len(mass) - 1)
    masses_grad = mass[ib_new: iend + 1]

    h1_grad = [h_grad + (m - mflat)/(menv - mflat) *
               (iso_profiles['h1'][iend] - h_grad)
               for m in masses_grad]

    # print(h_surf - arb_surf)

    h1_flat = [h_surf - arb_surf for i in range(0, ib_new)]
    h1_core = iso_profiles['h1'][iend + 1:]
    h1_prof = np.concatenate((h1_flat, h1_grad, h1_core))

    arb_flat = [arb_surf for i in range(0, ib_new)]
    arb_rest = [1e-20 for x in iso_profiles[arbn][ib_new:]]
    arb_prof = np.concatenate((arb_flat, arb_rest))

    metals = [x for x in isotopes if x not in ['h1', 'he4']]

    z_prof = np.array([iso_profiles[x] for x in metals])
    zsum = np.sum(z_prof, axis=0)
    he4_prof = 1 - h1_prof - arb_prof - zsum

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h1_prof
        elif isotope == 'he4':
            return he4_prof
        elif isotope == arbn:
            return arb_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_cno_hshell(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Changes the CNO abundance in the H burning shell
    """
    # Edge of He core
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-2), len(mass) - 1)

    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']
    # frac is the inner fraction of the envelope to modify
    frac = kwargs['frac']

    menv_change = frac * mass[iend]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > mass[iend] - menv_change), 0)
    cno = ['c12', 'c13', 'n13', 'n14', 'n15', 'o14', 'o15', 'o16', 'o17',
           'o18']

    cno_prof = np.array([iso_profiles[x] for x in cno])
    cnosum = np.sum(cno_prof, axis=0)
    cnonew = cnosum * z_fact
    diff = cnosum - cnonew

    h_prof_out = iso_profiles['h1'][:ichange]

    h_prof_mod = np.add(iso_profiles['h1'][ichange:iend], diff[ichange:iend])
    h_prof_core = iso_profiles['h1'][iend:]
    h_prof = np.concatenate((h_prof_out, h_prof_mod, h_prof_core))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h_prof
        elif isotope not in cno:
            return iso_profiles[isotope]
        else:
            iso_env = iso_profiles[isotope][:ichange]
            iso_mod = iso_profiles[isotope][ichange:iend] * z_fact
            iso_core = iso_profiles[isotope][iend:]
            return np.concatenate((iso_env, iso_mod, iso_core))

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_cno_envelope(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Changes the CNO abundance only in the envelope and not in the H burning
    shell
    """
    # Edge of He core
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-2), len(mass) - 1)

    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']
    # frac is the inner fraction of the envelope to leave as it is
    frac = kwargs['frac']

    menv_change = frac * mass[iend]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > mass[iend] - menv_change), 0)
    cno = ['c12', 'c13', 'n13', 'n14', 'n15', 'o14', 'o15', 'o16', 'o17',
           'o18']

    cno_prof = np.array([iso_profiles[x] for x in cno])
    cnosum = np.sum(cno_prof, axis=0)
    cnonew = cnosum * z_fact
    diff = cnosum - cnonew

    # h_prof_out = iso_profiles['h1'][:ichange]

    h_prof_mod = np.add(iso_profiles['h1'][:ichange], diff[:ichange])
    h_prof_core = iso_profiles['h1'][ichange:]
    h_prof = np.concatenate((h_prof_mod, h_prof_core))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h_prof
        elif isotope not in cno:
            return iso_profiles[isotope]
        else:
            # iso_env = iso_profiles[isotope][]
            iso_mod = iso_profiles[isotope][:ichange] * z_fact
            iso_core = iso_profiles[isotope][ichange:]
            return np.concatenate((iso_mod, iso_core))

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_metals_envelope(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Changes the CNO abundance only in the envelope and not in the H burning
    shell
    """
    # Edge of He core
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-2), len(mass) - 1)

    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']
    # frac is the inner fraction of the envelope to leave as it is
    frac = kwargs['frac']

    menv_change = frac * mass[iend]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > mass[iend] - menv_change), 0)
    metals = [x for x in isotopes if x not in ['h1', 'he4']]

    metal_prof = np.array([iso_profiles[x] for x in metals])
    metalsum = np.sum(metal_prof, axis=0)
    metalnew = metalsum * z_fact
    diff = metalsum - metalnew

    # h_prof_out = iso_profiles['h1'][:ichange]

    h_prof_mod = np.add(iso_profiles['h1'][:ichange], diff[:ichange])
    h_prof_core = iso_profiles['h1'][ichange:]
    h_prof = np.concatenate((h_prof_mod, h_prof_core))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h_prof
        elif isotope not in metals:
            return iso_profiles[isotope]
        else:
            # iso_env = iso_profiles[isotope][]
            iso_mod = iso_profiles[isotope][:ichange] * z_fact
            iso_core = iso_profiles[isotope][ichange:]
            return np.concatenate((iso_mod, iso_core))

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_cno_envelope(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Changes the CNO abundance only in the envelope and not in the H burning
    shell
    """
    # Edge of He core
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-2), len(mass) - 1)

    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']
    # frac is the inner fraction of the envelope to leave as it is
    frac = kwargs['frac']

    menv_change = frac * mass[iend]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > mass[iend] - menv_change), 0)
    cno = ['c12', 'c13', 'n13', 'n14', 'n15', 'o14', 'o15', 'o16', 'o17',
           'o18']

    cno_prof = np.array([iso_profiles[x] for x in cno])
    cnosum = np.sum(cno_prof, axis=0)
    cnonew = cnosum * z_fact
    diff = cnosum - cnonew

    # h_prof_out = iso_profiles['h1'][:ichange]

    h_prof_mod = np.add(iso_profiles['h1'][:ichange], diff[:ichange])
    h_prof_core = iso_profiles['h1'][ichange:]
    h_prof = np.concatenate((h_prof_mod, h_prof_core))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h_prof
        elif isotope not in cno:
            return iso_profiles[isotope]
        else:
            # iso_env = iso_profiles[isotope][]
            iso_mod = iso_profiles[isotope][:ichange] * z_fact
            iso_core = iso_profiles[isotope][ichange:]
            return np.concatenate((iso_mod, iso_core))

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_add_arb2(iso_profiles, isotopes, mass, dm, starmass, arbn='arb2',
                **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient as specified
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # arb_frac is the fraction of Menv that contains arb1
    arb_frac = kwargs['arb_frac']
    # h_surf = kwargs['h_surface']
    arb_surf = kwargs['arb_surf']

    # mass of the envelope
    menv = mass[iend]
    mflat = menv * arb_frac
    ib_new = next((x[0] for x in enumerate(mass)
                  if x[1] > mflat), len(mass) - 1)

    arb_flat = [arb_surf for i in range(0, ib_new)]
    arb_rest = [1e-20 for x in iso_profiles[arbn][ib_new:]]
    arb_prof = np.concatenate((arb_flat, arb_rest))

    # arb_flat_diff = np.subtract(iso_profiles['h1'][:ib_new], arb_flat)

    he4_flat = np.subtract(iso_profiles['he4'][:ib_new], arb_flat)
    he4_rest = iso_profiles['he4'][ib_new:]
    he4_prof = np.concatenate((he4_flat, he4_rest))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'he4':
            return he4_prof
        elif isotope == arbn:
            return arb_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_add_arb1(iso_profiles, isotopes, mass, dm, starmass, arbn='arb1',
                **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient as specified
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # arb_frac is the fraction of Menv that contains arb1
    arb_frac = kwargs['arb_frac']
    # h_surf = kwargs['h_surface']
    arb_surf = kwargs['arb_surf']

    # mass of the envelope
    menv = mass[iend]
    mflat = menv * arb_frac
    ib_new = next((x[0] for x in enumerate(mass)
                  if x[1] > mflat), len(mass) - 1)

    arb_flat = [arb_surf for i in range(0, ib_new)]
    arb_rest = [1e-20 for x in iso_profiles[arbn][ib_new:]]
    arb_prof = np.concatenate((arb_flat, arb_rest))

    # arb_flat_diff = np.subtract(iso_profiles['h1'][:ib_new], arb_flat)

    he4_flat = np.subtract(iso_profiles['he4'][:ib_new], arb_flat)
    he4_rest = iso_profiles['he4'][ib_new:]
    he4_prof = np.concatenate((he4_flat, he4_rest))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'he4':
            return he4_prof
        elif isotope == arbn:
            return arb_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds
# list of different ways to modify models


def mod_add_h1(iso_profiles, isotopes, mass, dm, starmass, arbn='h1',
                **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient as specified
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # arb_frac is the fraction of Menv that contains arb1
    arb_frac = kwargs['arb_frac']
    # h_surf = kwargs['h_surface']
    arb_surf = kwargs['arb_surf']

    # mass of the envelope
    menv = mass[iend]
    mflat = menv * arb_frac
    ib_new = next((x[0] for x in enumerate(mass)
                  if x[1] > mflat), len(mass) - 1)

    arb_flat = [arb_surf for i in range(0, ib_new)]
    arb_rest = [1e-20 for x in iso_profiles[arbn][ib_new:]]
    arb_prof = iso_profiles['h1'] + np.concatenate((arb_flat, arb_rest))

    # arb_flat_diff = np.subtract(iso_profiles['h1'][:ib_new], arb_flat)

    he4_flat = np.subtract(iso_profiles['he4'][:ib_new], arb_flat)
    he4_rest = iso_profiles['he4'][ib_new:]
    he4_prof = np.concatenate((he4_flat, he4_rest))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'he4':
            return he4_prof
        elif isotope == arbn:
            return arb_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_metals_envelope_tohe(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Changes the CNO abundance only in the envelope and not in the H burning
    shell
    """
    # Edge of He core
    iend = next((x[0] for x in enumerate(iso_profiles['h1'])
                if x[1] < 1e-2), len(mass) - 1)

    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']
    # frac is the inner fraction of the envelope to leave as it is
    frac = kwargs['frac']

    menv_change = frac * mass[iend]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > mass[iend] - menv_change), 0)
    metals = [x for x in isotopes if x not in ['h1', 'he4']]

    metal_prof = np.array([iso_profiles[x] for x in metals])
    metalsum = np.sum(metal_prof, axis=0)
    metalnew = metalsum * z_fact
    diff = metalsum - metalnew

    # h_prof_out = iso_profiles['h1'][:ichange]

    he_prof_mod = np.add(iso_profiles['he4'][:ichange], diff[:ichange])
    he_prof_core = iso_profiles['he4'][ichange:]
    he_prof = np.concatenate((he_prof_mod, he_prof_core))

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'he4':
            return he_prof
        elif isotope not in metals:
            return iso_profiles[isotope]
        else:
            # iso_env = iso_profiles[isotope][]
            iso_mod = iso_profiles[isotope][:ichange] * z_fact
            iso_core = iso_profiles[isotope][ichange:]
            return np.concatenate((iso_mod, iso_core))

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_grad_const_hemass_popIII(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    # modifies gradient at constant total H and He Mass
    def get_ind_steep(itry_a, oldabund, indices=False):
        # returns list of total helium masses and corresponding indices
        f_abund = interpolate.interp1d(norm_old_mass, oldabund[ibegin:iend],
                                       kind='cubic', fill_value='extrapolate')
        newendmass = mass[itry_a] + newdeltam
        itry_b = next((x[0] for x in enumerate(mass) if x[1] > newendmass), 0)
        norm_new_mass = (mass[itry_a:itry_b] - mass[itry_a]) / newdeltam

        p1 = np.array((itry_a - ibegin) * [oldabund[0]])
        p2 = oldabund[:ibegin]
        p3 = f_abund(norm_new_mass)
        iextra = iend - itry_b
        if iextra > 0:
            p4 = oldabund[iend:]
            p5 = iextra * [oldabund[-1]]
        else:
            p4 = oldabund[itry_b:]
            p5 = []

        # psizes = [len(list(p)) for p in (p1, p2, p3, p4, p5)]
        newabund = np.concatenate((p1, p2, p3, p4, p5))
        if indices:
            return [itry_a, itry_b, np.dot(dm, newabund)]
        else:
            return newabund

    def get_ind_shallow(itry_b, oldabund, indices=False):
        # returns list of total helium masses and corresponding indices
        f_abund = interpolate.interp1d(norm_old_mass, oldabund[ibegin:iend],
                                       kind='cubic', fill_value='extrapolate')
        newendmass = mass[itry_b] - newdeltam
        itry_a = next((x[0] for x in enumerate(mass) if x[1] > newendmass), 0)
        norm_new_mass = (mass[itry_a:itry_b] - mass[itry_a]) / newdeltam

        iextra = itry_a - ibegin
        if iextra > 0:
            p1 = iextra * [oldabund[0]]
            p2 = oldabund[:ibegin]
            p3 = f_abund(norm_new_mass)
            p4 = oldabund[itry_b:]
        else:
            p1 = oldabund[ibegin - itry_a:ibegin]
            p1t = oldabund[:ibegin - itry_a]
            p2 = f_abund(norm_new_mass)
            p3 = oldabund[itry_b:]
            p4 = []
            # if indices and itry_a == 2031:
            #     print(len(p1), len(p1t))

        newabund = np.concatenate((p1, p2, p3, p4))
        if indices:
            return [itry_a, itry_b, np.dot(dm, newabund)]
        else:
            return newabund

    def get_abunds():
        if slopefactor >= 1:
            ias, ibs, h1_masses = np.array([get_ind_steep(itry, old_h1,
                                            indices=True) for itry in
                                            range(ibegin, iend)]).T
            ia, ib, m = next(((int(ia), int(ib), m)
                             for ia, ib, m in zip(ias, ibs, h1_masses)
                             if m >= total_h1_mass), (ibegin, iend,
                             total_h1_mass))
            # print(m)
            newabunds = [get_ind_steep(ia, iso_profiles[iso], indices=False)
                         for iso in isotopes]
        elif slopefactor < 1:
            print(ibegin, iend)
            ias, ibs, h1_masses = np.array([get_ind_shallow(itry, old_h1,
                                            indices=True) for itry in
                                            range(iend, len(old_h1))]).T
            ia, ib, m = next(((int(ia), int(ib), m) for ia, ib, m in zip(ias,
                             ibs, h1_masses) if m >= total_h1_mass),
                             (ibegin, iend, total_h1_mass))
            print(h1_masses)
            newabunds = [get_ind_shallow(ib, iso_profiles[iso], indices=False)
                         for iso in isotopes]

        return newabunds

    ibegin, iend = get_i_hcore(iso_profiles, mass)
    print(mass[-1] - mass[ibegin], mass[-1] - mass[iend])
    # print(ibegin, iend, len(mass), mass[ibegin], mass[iend])
    # print(iso_profiles['h1'][ibegin], iso_profiles['h1'][iend], len(mass), mass[-1] - mass[ibegin], mass[-1] - mass[iend])
    # print(mass[ibegin], mass[iend])
    slopefactor = kwargs['slopefactor']
    total_he4_mass = np.dot(dm, iso_profiles['he4'])
    total_h1_mass = np.dot(dm, iso_profiles['h1'])
    old_h1 = iso_profiles['h1']
    old_he4 = iso_profiles['he4']

    # original slope
    deltam_gradient = mass[iend] - mass[ibegin]
    deltah_gradient = iso_profiles['h1'][iend] - iso_profiles['h1'][ibegin]
    oldslope = -1 * deltah_gradient / deltam_gradient

    # figure out new slope
    newslope = oldslope * slopefactor
    newdeltam = -1 * deltah_gradient / newslope
    print(newdeltam)
    itry_b = 2183

    newendmass = mass[itry_b] - newdeltam
    itry_a = next((x[0] for x in enumerate(mass) if x[1] > newendmass), 0)
    norm_new_mass = (mass[itry_a:itry_b] - mass[itry_a]) / newdeltam

    iextra = itry_a - ibegin
    print(iextra)

    # normalised mass in gradient region
    norm_old_mass = (mass[ibegin:iend] - mass[ibegin])/deltam_gradient
    newabunds = get_abunds()
    return newabunds


def mod_ms_metals_envelope(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes the metallicity in a MS star for a given fraction of the star
    beginning at the surface
    kwargs:
        - fix: either 'h1' or 'he4', depending on which one you want to keep
               constant (as the total depends on metallicity)
        - abund_surf: surface abundance of fix

    Only works for MS models
    """
    # z_fact is the factor to scale the metals by
    z_fact = kwargs['z_fact']
    # frac is the outer fraction of the envelope to modify
    frac = kwargs['frac']

    menv_change = frac * mass[-1]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > menv_change), 0)
    cno = ['c12', 'c13', 'n13', 'n14', 'n15', 'o14', 'o15', 'o16', 'o17',
           'o18']

    metals = [x for x in isotopes if x not in ['h1', 'he4']]

    cno_prof = np.array([iso_profiles[x] for x in metals])
    cnosum = np.sum(cno_prof, axis=0)
    cnonew = cnosum * z_fact
    diff = cnosum - cnonew

    h_prof_mod = np.add(iso_profiles['h1'][:ichange], diff[:ichange])
    h_prof_in = iso_profiles['h1'][ichange:]
    h_prof = np.concatenate((h_prof_mod, h_prof_in))

    he_prof = iso_profiles['he4']

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h_prof
        elif isotope not in metals:
            return iso_profiles[isotope]
        else:
            iso_env = iso_profiles[isotope][:ichange] * z_fact
            iso_mod = iso_profiles[isotope][ichange:]
            return np.concatenate((iso_env, iso_mod))
    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def add_arb_to_solar(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    adds arb1 to the .mod file
    """
    # isotopes.append('arb1')
    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'arb1':
            return np.zeros(len(mass))
        else:
            return iso_profiles[isotope]
    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds


def mod_fuel_supply(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies fuel supply in the core
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # core_frac is the mass fraction of the core
    core_frac = kwargs['core_frac']
    # arb_core is the fraction of arb1 in the core
    arb_core = kwargs['arb_core']

    # mass of the envelope
    menv_change = core_frac * mass[-1]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > mass[-1] - menv_change), 0)

    arb_rest = [1e-20 for x in iso_profiles['arb1'][:ichange]]
    arb_core = [arb_core for x in iso_profiles['arb1'][ichange:]]
    arb_prof = np.concatenate((arb_rest, arb_core))

    h1_new = iso_profiles['h1'] - arb_prof

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h1_new
        elif isotope == 'arb1':
            return arb_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds

def mod_mucore(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies fuel supply in the core
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # core_frac is the mass fraction of the core
    core_frac = kwargs['core_frac']
    # arb_core is the fraction of arb1 in the core
    arb_core = kwargs['arb_core']

    # mass of the envelope
    menv_change = core_frac * mass[-1]
    ichange = next((x[0] for x in enumerate(mass)
                    if x[1] > mass[-1] - menv_change), 0)

    arb_rest = [1e-20 for x in iso_profiles['arb1'][:ichange]]
    arb_core = [arb_core for x in iso_profiles['arb1'][ichange:]]
    arb_prof = np.concatenate((arb_rest, arb_core))

    he4_new = iso_profiles['he4'] - arb_prof

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'he4':
            return he4_new
        elif isotope == 'arb1':
            return arb_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_consthe_homogen(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    # produces homogenous star at constant total H and He Mass

    ibegin, iend = get_i_hcore(iso_profiles, mass)
    total_he4_mass = np.dot(dm, iso_profiles['he4'])
    total_h1_mass = np.dot(dm, iso_profiles['h1'])
    old_h1 = iso_profiles['h1']
    old_he4 = iso_profiles['he4']
    new_he4 = [total_he4_mass/starmass for x in old_he4]
    new_h1 = [total_h1_mass/starmass for x in old_h1]

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'he4':
            return new_he4
        elif isotope == 'h1':
            return new_h1
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]
    return newabunds



def mod_mass_gainer_heburning(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    changes just the H & He gradient above the H depleted core in an absolute
    way that sets the gradient as specified
    """
    ibegin, iend = get_i_hcore(iso_profiles, mass)
    # grad_frac is the fraction of Menv that consists of a Gradient
    frac = kwargs['frac']
    h_surf = kwargs['h_surf']

    # mass of the envelope
    menv = mass[iend]
    mflat = menv * (frac)
    ichange = next((x[0] for x in enumerate(mass)
                   if x[1] > mflat), len(mass) - 1)
    # masses_grad = mass[ib_new: iend + 1]

    h_prof_mod = [h_surf for x in iso_profiles['h1'][:ichange]]
    h_prof_in = iso_profiles['h1'][ichange:]
    h_prof = np.concatenate((h_prof_mod, h_prof_in))

    he_prof_mod = np.add(iso_profiles['he4'][:ichange], iso_profiles['h1'][:ichange] - h_surf)
    he_prof_in = iso_profiles['he4'][ichange:]
    he_prof = np.concatenate((he_prof_mod, he_prof_in))


    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        if isotope == 'h1':
            return h_prof
        elif isotope == 'he4':
            return he_prof
        else:
            return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds



def mod_coreratio(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies the core mass ratio
    """
    xmin = iso_profiles['h1'][-1] + 1e-4
    ibegin = next((x[0] for x in enumerate(iso_profiles['h1'])
                  if x[1] < iso_profiles['h1'][0]-0.01), len(mass) - 1)
    iend = next((x[0] for x in enumerate(iso_profiles['h1'][ibegin:])
                if x[1] < xmin), len(mass) - 1)
    if iend < len(mass) - 1:
        iend += ibegin
    
    # menv is the input envelope mass. 
    menv = kwargs['menv']
    menv_actual = mass[iend]
    menv_diff = menv_actual - menv
    
    normalised_mass = mass/starmass
    
    if menv_diff > 0:
        ichange = next((x[0] for x in enumerate(mass)
                       if x[1] > menv_diff), len(mass) - 1)
        mass_new = (mass[ichange:] - menv_diff)/(starmass - menv_diff)
    else:
        print("Error menv_diff < 0")


    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        # f_init = interpolate.interp1d(normalised_mass, iso_profiles[isotope], fill_value='extrapolate')
        iso = iso_profiles[isotope][ichange:]
        f_new = interpolate.interp1d(mass_new, iso, fill_value='extrapolate')
        new_prof = f_new(normalised_mass)
        
        return new_prof

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_mcore_ms(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies "convective core" mass ratio of a MS star
    it's not really the convective part, more the part that has a flat abundance
    profile in the core region.
    The rest of the abundance profiles are just pushed outwards
    """

    # mcore is the new convective core mass ratio
    mcore = kwargs['mcore']

    # Get the current convective core mass ratio
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.01

    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mcore_old = mass[icore]/mass[0]

    print(mcore, mcore_old)
    
    icore_new = next((x[0] for x in enumerate(mass) if x[1] < mcore * mass[0]))


    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        iso_prof = iso_profiles[isotope]

        if mcore > mcore_old:

            # increase size of core
            nsurf_to_remove = icore - icore_new

            prof_outer = iso_prof[:icore]
            mass_outer = mass[:icore]

            prof_core = iso_prof[icore:]
            mass_core = mass[icore:]

            mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
            mass_fit_core = (mass_core - mass_core[0])/(mass_core[-1] - mass_core[0])

            fout = interpolate.interp1d(mass_fit_outer, prof_outer)
            fcore = interpolate.interp1d(mass_fit_core, prof_core)

            mnew_outer = mass[:icore - nsurf_to_remove]
            mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

            mnew_core = mass[icore - nsurf_to_remove:-nsurf_to_remove]
            mnew_core = (mnew_core - mnew_core[0])/(mnew_core[-1] - mnew_core[0])
            
            prof_outer_new = fout(mnew_outer)
            prof_core_new = fcore(mnew_core)

            prof_inner = np.full(nsurf_to_remove, iso_prof[-1])

            prof_new = np.concatenate((prof_outer_new, prof_core_new, prof_inner))

            # print(len(prof_outer_new), len(prof_core_new), len(prof_inner))
            # print(icore, icore_new)
            # print(len(mass), len(prof_new))


        elif mcore < mcore_old:

            # decrease size of core
            nsurf_to_remove = icore_new - icore

            prof_outer = iso_prof[:icore]
            mass_outer = mass[:icore]

            prof_core = iso_prof[icore:]
            mass_core = mass[icore:]

            mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
            mass_fit_core = (mass_core - mass_core[0])/(mass_core[-1] - mass_core[0])

            fout = interpolate.interp1d(mass_fit_outer, prof_outer)
            fcore = interpolate.interp1d(mass_fit_core, prof_core)

            mnew_outer = mass[nsurf_to_remove:icore + nsurf_to_remove]
            mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

            mnew_core = mass[icore + nsurf_to_remove:]
            mnew_core = (mnew_core - mnew_core[0])/(mnew_core[-1] - mnew_core[0])
            
            prof_outer_new = fout(mnew_outer)
            prof_core_new = fcore(mnew_core)

            prof_surface = np.full(nsurf_to_remove, iso_prof[-1])

            prof_new = np.concatenate((prof_surface, prof_outer_new, prof_core_new))


            # # decrease size of core
            # nsurf_to_remove = icore_new - icore

            # prof_outer = iso_prof[:-nsurf_to_remove]
            # mass_outer = mass[:-nsurf_to_remove]
            # mass_fit = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])

            # f = interpolate.interp1d(mass_fit, prof_outer)

            # mnew = mass[nsurf_to_remove:]
            # mnew = (mnew - mnew[0])/(mnew[-1] - mnew[0])

            # prof_outer_new = f(mnew)

            # prof_surface = np.full(nsurf_to_remove, iso_prof[0])

            # prof_new = np.concatenate((prof_surface, prof_outer_new))


        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_slope_ms(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies slope of hydrogen profile of a MS star
    should maintain the same core and surface abundances and core mass
    where core is the region of flat abundance in the center
    """

    # mslope is the mass over which the slope changes
    mslope = kwargs['mslope']

    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))


    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d(mass_fit_slope, prof_slope)

        mnew_outer = np.array(mass[:ienv_new]).astype(float)
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = np.array(mass[ienv_new:icore]).astype(float)
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_straight_slope_ms(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies slope of hydrogen profile of a MS star
    should maintain the same core and surface abundances and core mass
    where core is the region of flat abundance in the center
    """

    # mslope is the mass over which the slope changes
    mslope = kwargs['mslope']

    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))


    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d([0, 1], [prof_slope[0], prof_slope[1]])

        # print(prof_slope[0], prof_slope[-1], fslope(0), fslope(0.5), fslope(1))

        mnew_outer = mass[:ienv_new]
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = mass[ienv_new:icore]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        

        mnew_outer = np.array(mnew_outer).astype(float)
        mnew_slope = np.array(mnew_slope).astype(float)
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)
        prof_slope_new = prof_slope[0] + mnew_slope * (prof_slope[-1] - prof_slope[0])

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


def mod_same(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Doesn't modify profile
    """

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """
        return iso_profiles[isotope]

    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds





def mod_3part_ms(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies slope of hydrogen profile of a MS star
    should maintain the same core and surface abundances and core mass
    where core is the region of flat abundance in the center
    """


    # mslope is the mass over which the slope changes
    if 'mcore' in kwargs:
        mcore = kwargs['mcore']
    else:
        # Get the current convective core mass ratio
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.01

        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mcore_old = mass[icore]/mass[0]

        mcore = mcore_old

    if 'mslope' in kwargs:
        mslope = kwargs['mslope']
    else:
        # Get the current mass over which slope changes
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.02

        ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

        h_core_diff = 0.002
        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mslope_old = mass[ienv] - mass[icore]
        mslope = mslope_old

    if 'xcore' in kwargs:
        xcore = kwargs['xcore']
    else:
        xcore = iso_profiles['h1'][-1]


    if 'xsurface' in kwargs:
        xsurface = kwargs['xsurface']
    else:
        xsurface = iso_profiles['h1'][0]



    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))

    
    # Compute new H1 profile here:

    icore_new = next((x[0] for x in enumerate(mass) if x[1] < mcore * mass[0]))
    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mslope) * mass[0]))

    h1_new_core = np.full(len(mass) - icore_new, xcore) 
    h1_new_surface = np.full(ienv_new, xsurface)

    fslope = interpolate.interp1d([0, 1], [xsurface, xcore])
    mnew_slope = mass[ienv_new:icore_new]
    mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
    mnew_slope = np.array(mnew_slope).astype(float)

    h1_slope_new = xsurface + mnew_slope * (xcore - xsurface)
    h1_slope_new = fslope(mnew_slope)

    h1_new = np.concatenate((h1_new_surface, h1_slope_new, h1_new_core))

    # print(len(h1_new_surface), len(h1_slope_new), len(h1_new_core))
    # print(len(h1_new), len(mass))
    # print(icore_new, ienv_new)
    # print(len(mass) - icore_new)



    # mslope = kwargs['mslope']

    # # Get the current mass over which slope changes
    # h_prof = iso_profiles['h1']
    # h_prof_rev = h_prof[::-1]
    # h_core_diff = 0.02

    # ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    # h_core_diff = 0.002
    # icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    # icore = len(h_prof) - icore_rev - 1

    # mslope_old = mass[ienv] - mass[icore]

    # ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))


    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        if isotope == 'h1':
            return h1_new
        elif isotope == 'he4':
            return iso_profiles['he4'] + (iso_profiles['h1'] - h1_new)

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d([0, 1], [prof_slope[0], prof_slope[1]])

        mnew_outer = mass[:ienv_new]
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = mass[ienv_new:icore]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        
        mnew_outer = np.array(mnew_outer).astype(float)
        mnew_slope = np.array(mnew_slope).astype(float)
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)
        prof_slope_new = prof_slope[0] + mnew_slope * (prof_slope[-1] - prof_slope[0])

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds



def mod_3part_ms_curved(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies slope of hydrogen profile of a MS star
    should maintain the same core and surface abundances and core mass
    where core is the region of flat abundance in the center
    """


    # mslope is the mass over which the slope changes
    if 'mcore' in kwargs:
        mcore = kwargs['mcore']
    else:
        # Get the current convective core mass ratio
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.01

        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mcore_old = mass[icore]/mass[0]

        mcore = mcore_old

    if 'mslope' in kwargs:
        mslope = kwargs['mslope']
    else:
        # Get the current mass over which slope changes
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.02

        ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

        h_core_diff = 0.002
        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mslope_old = mass[ienv] - mass[icore]
        mslope = mslope_old

    if 'xcore' in kwargs:
        xcore = kwargs['xcore']
    else:
        xcore = iso_profiles['h1'][-1]


    if 'xsurface' in kwargs:
        xsurface = kwargs['xsurface']
    else:
        xsurface = iso_profiles['h1'][0]

    if 'fa' in kwargs:
        fa = kwargs['fa']
        ma = mcore + fa * (mslope)
    else:
        ma = 0

    if 'fb' in kwargs:
        fb = kwargs['fb']
        mb = mcore + mslope + fb * (1 - mslope - mcore)
    else:
        mb = 0


    if 'alpha_curve' in kwargs:
        alpha_curve = kwargs['alpha_curve']
    else:
        alpha_curve = 0


    # print('fa, fb, ma, mb', fa, fb, ma, mb, mcore, mslope, mass[0])

    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))

    
    # Compute new H1 profile here:

    icore_new = next((x[0] for x in enumerate(mass) if x[1] < mcore * mass[0]))
    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mslope) * mass[0]))

    h1_new_core = np.full(len(mass) - icore_new, xcore) 
    h1_new_surface = np.full(ienv_new, xsurface)

    fslope = interpolate.interp1d([0, 1], [xsurface, xcore])
    mnew_slope = mass[ienv_new:icore_new]
    mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
    mnew_slope = np.array(mnew_slope).astype(float)

    h1_slope_new = xsurface + mnew_slope * (xcore - xsurface)
    h1_slope_new = fslope(mnew_slope)

    h1_new = np.concatenate((h1_new_surface, h1_slope_new, h1_new_core))


    if ma != 0 and mb != 0:
        ia = next((x[0] for x in enumerate(mass) if x[1] < ma * mass[0]))
        ib = next((x[0] for x in enumerate(mass) if x[1] < mb * mass[0]))

        x1 = mass[ia]/mass[0]
        x3 = mass[ib]/mass[0]
        y1 = h1_new[ia]
        y3 = h1_new[ib]

        x4 = mcore + mslope
        y4 = xsurface

        # Straight line method
        # c1 = h1_new[ib:ia]
        # fstraight = interpolate.interp1d([x1, x3], [y1, y3], kind='linear', fill_value='extrapolate')
        # y2min = np.array([fstraight(x) for x in mass[ib:ia]/mass[0]])
        # c2 = y2min

        # c_new = c1 * alpha_curve + c2 * (1 - alpha_curve)

        # h1_new = np.concatenate((h1_new[:ib], c_new, h1_new[ia:]))

        # a = np.array([y1, x1])
        # b = np.array([y4, x4])
        # c = np.array([y3, x3])

        # ba = a - b
        # bc = c - b

        # cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        # angle = np.arccos(cosine_angle)

        # # angle = np.degrees(angle)

        # print("Angle degrees:", np.degrees(angle))

        # angle_new = 0.5 * np.abs(angle)

        # print("Angle new degrees:", np.degrees(angle_new))

        # dx = alpha_curve * np.cos(angle)
        # dy = alpha_curve * np.sin(angle)

        # x2 = x4 - dx
        # y2 = y4 - dy

        fstraight = interpolate.interp1d([x1, x3], [y1, y3], kind='linear', fill_value='extrapolate')
        # y2min = fstraight(x2)
        # x2 = x4
        # y2max = y4
        # y2min = fstraight(x2)
        # y2 = y2max - alpha_curve * (y2max - y2min)
        
        # x2 = 0.5 * (x1 + x3)
        # ix2 = next((x[0] for x in enumerate(mass) if x[1] < x2 * mass[0]))
        ix2 = ib + 1
        x2 = mass[ix2] / mass[0]
        y2max = h1_new[ix2]
        y2min = fstraight(x2)
        y2 = y2max - alpha_curve * (y2max - y2min)




        # y2 = np.maximum(y2, y2min)

        # print('x1, y1', x1*mass[0], y1)
        # print('x2, y2', x2*mass[0], y2)
        # print('x3, y3', x3*mass[0], y3)
        # print()
        # print('x4, y4', x4*mass[0], y4)

        # print('x1, x2, x3, x4', x1 * mass[0], x2 * mass[0], x3 * mass[0], x4 * mass[0])
        # print('y1, y2, y3, y4', y1, y2, y3, y4)

        fcurve = interpolate.interp1d([x1, x2, x3], [y1, y2, y3], kind='quadratic', fill_value='extrapolate')

        ycurve = np.array([fcurve(x) for x in mass[ib:ia]/mass[0]])

        ycurve = np.minimum(ycurve, xsurface)

        h1_new = np.concatenate((h1_new[:ib], ycurve, h1_new[ia:]))
       

        # angle_1 = np.dot()

        # a0 = (x1 * x2**2 - x1**2 * x2 + x1 * y2 - x2 * y1)/(x1 - x2)

        # a1 = (x1**2 - x2**2 + y1 - y2)/(x1 - x2)

        # # derivative at x1 is equal to slope
        # dhdm = (h1_new[ia] - h1_new[ia - 1])/(mass[ia] - mass[ia - 1])

        # w = dhdm / (a1 - 2 * ma)

        # w = 1
        # a0 = w * a0
        # a1 = w * a1
        # a2 = -1

        # # equation is a2 * x^2 + a1 * x + a0

        # norm_mass = mass/mass[0]
        # def get_curved_part(mass):
        #     return a2 * mass ** 2 + a1 * mass + a0

        # mass_curved = norm_mass[ib:ia]
        # h_curved = get_curved_part(mass_curved)

        # h1_new[ib:ia] = h_curved


    # print(len(h1_new_surface), len(h1_slope_new), len(h1_new_core))
    # print(len(h1_new), len(mass))
    # print(icore_new, ienv_new)
    # print(len(mass) - icore_new)



    # mslope = kwargs['mslope']

    # # Get the current mass over which slope changes
    # h_prof = iso_profiles['h1']
    # h_prof_rev = h_prof[::-1]
    # h_core_diff = 0.02

    # ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    # h_core_diff = 0.002
    # icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    # icore = len(h_prof) - icore_rev - 1

    # mslope_old = mass[ienv] - mass[icore]

    # ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))


    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        if isotope == 'h1':
            return h1_new
        elif isotope == 'he4':
            return iso_profiles['he4'] + (iso_profiles['h1'] - h1_new)

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d([0, 1], [prof_slope[0], prof_slope[1]])

        mnew_outer = mass[:ienv_new]
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = mass[ienv_new:icore]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        
        mnew_outer = np.array(mnew_outer).astype(float)
        mnew_slope = np.array(mnew_slope).astype(float)
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)
        prof_slope_new = prof_slope[0] + mnew_slope * (prof_slope[-1] - prof_slope[0])

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds



def mod_3part_ms_curved_bot(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies slope of hydrogen profile of a MS star
    should maintain the same core and surface abundances and core mass
    where core is the region of flat abundance in the center

    Going to be different to mod_3part_ms_curved above as curved region
    will be at edge of core
    """


    # mslope is the mass over which the slope changes
    if 'mcore' in kwargs:
        mcore = kwargs['mcore']
    else:
        # Get the current convective core mass ratio
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.01

        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mcore_old = mass[icore]/mass[0]

        mcore = mcore_old

    if 'mslope' in kwargs:
        mslope = kwargs['mslope']
    else:
        # Get the current mass over which slope changes
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.02

        ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

        h_core_diff = 0.002
        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mslope_old = mass[ienv] - mass[icore]
        mslope = mslope_old

    if 'xcore' in kwargs:
        xcore = kwargs['xcore']
    else:
        xcore = iso_profiles['h1'][-1]


    if 'xsurface' in kwargs:
        xsurface = kwargs['xsurface']
    else:
        xsurface = iso_profiles['h1'][0]

    if 'fa' in kwargs:
        fa = kwargs['fa']
        ma = fa * mcore
    else:
        ma = 0

    if 'fb' in kwargs:
        fb = kwargs['fb']
        mb = mcore + fb * (mslope)
    else:
        mb = 0


    if 'alpha_curve' in kwargs:
        alpha_curve = kwargs['alpha_curve']
    else:
        alpha_curve = 0


    # print('fa, fb, ma, mb', fa, fb, ma, mb, mcore, mslope, mass[0])

    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))

    
    # Compute new H1 profile here:

    icore_new = next((x[0] for x in enumerate(mass) if x[1] < mcore * mass[0]))
    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mslope) * mass[0]))

    h1_new_core = np.full(len(mass) - icore_new, xcore) 
    h1_new_surface = np.full(ienv_new, xsurface)

    fslope = interpolate.interp1d([0, 1], [xsurface, xcore])
    mnew_slope = mass[ienv_new:icore_new]
    mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
    mnew_slope = np.array(mnew_slope).astype(float)

    h1_slope_new = xsurface + mnew_slope * (xcore - xsurface)
    h1_slope_new = fslope(mnew_slope)

    h1_new = np.concatenate((h1_new_surface, h1_slope_new, h1_new_core))


    if ma != 0 and mb != 0:
        ia = next((x[0] for x in enumerate(mass) if x[1] < ma * mass[0]))
        ib = next((x[0] for x in enumerate(mass) if x[1] < mb * mass[0]))

        x1 = mass[ia]/mass[0]
        x3 = mass[ib]/mass[0]
        y1 = h1_new[ia]
        y3 = h1_new[ib]

        x4 = mcore
        y4 = xcore


        fstraight = interpolate.interp1d([x1, x3], [y1, y3], kind='linear', fill_value='extrapolate')

        ix2 = ib + 1
        x2 = mass[ix2] / mass[0]
        y2max = h1_new[ix2]
        y2min = fstraight(x2)
        y2 = y2max - alpha_curve * (y2max - y2min)



        fcurve = interpolate.interp1d([x1, x2, x3], [y1, y2, y3], kind='quadratic', fill_value='extrapolate')

        ycurve = np.array([fcurve(x) for x in mass[ib:ia]/mass[0]])

        ycurve = np.minimum(ycurve, xsurface)

        h1_new = np.concatenate((h1_new[:ib], ycurve, h1_new[ia:]))
       

    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        if isotope == 'h1':
            return h1_new
        elif isotope == 'he4':
            return iso_profiles['he4'] + (iso_profiles['h1'] - h1_new)

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d([0, 1], [prof_slope[0], prof_slope[1]])

        mnew_outer = mass[:ienv_new]
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = mass[ienv_new:icore]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        
        mnew_outer = np.array(mnew_outer).astype(float)
        mnew_slope = np.array(mnew_slope).astype(float)
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)
        prof_slope_new = prof_slope[0] + mnew_slope * (prof_slope[-1] - prof_slope[0])

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds



def mod_3part_ms_icz(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies slope of hydrogen profile of a MS star
    should maintain the same core and surface abundances and core mass
    where core is the region of flat abundance in the center
    """


    # mslope is the mass over which the slope changes
    if 'mcore' in kwargs:
        mcore = kwargs['mcore']
    else:
        # Get the current convective core mass ratio
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.01

        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mcore_old = mass[icore]/mass[0]

        mcore = mcore_old

    if 'mslope' in kwargs:
        mslope = kwargs['mslope']
    else:
        # Get the current mass over which slope changes
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.02

        ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

        h_core_diff = 0.002
        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mslope_old = mass[ienv] - mass[icore]
        mslope = mslope_old

    if 'xcore' in kwargs:
        xcore = kwargs['xcore']
    else:
        xcore = iso_profiles['h1'][-1]


    if 'xsurface' in kwargs:
        xsurface = kwargs['xsurface']
    else:
        xsurface = iso_profiles['h1'][0]

    if 'mcz' in kwargs:
        mcz = kwargs['mcz']
    else:
        mcz = None

    if 'wcz' in kwargs:
        wcz = kwargs['wcz']
    else:
        wcz = None

    if 'xcz' in kwargs:
        xcz = kwargs['xcz']
    else:
        xcz = None




    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))

    
    # Compute new H1 profile here:

    icore_new = next((x[0] for x in enumerate(mass) if x[1] < mcore * mass[0]))
    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mslope) * mass[0]))

    h1_new_core = np.full(len(mass) - icore_new, xcore) 
    h1_new_surface = np.full(ienv_new, xsurface)

    # Compute slope
    if mcz is None or wcz is None or xcz is None:
        fslope = interpolate.interp1d([0, 1], [xsurface, xcore])
        mnew_slope = mass[ienv_new:icore_new]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        mnew_slope = np.array(mnew_slope).astype(float)

        h1_slope_new = xsurface + mnew_slope * (xcore - xsurface)
        h1_slope_new = fslope(mnew_slope)
    else:
        i1 = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mcz) * mass[0]))
        i2 = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mcz + wcz) * mass[0]))
        
        fslope = interpolate.interp1d([0, 1], [xsurface, xcz])
        mnew_slope = mass[ienv_new:i2]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        mnew_slope = np.array(mnew_slope).astype(float)

        # h1_slope_new1 = xsurface + mnew_slope * (xcz - xsurface)
        h1_slope_new1 = fslope(mnew_slope)        


        fslope = interpolate.interp1d([0, 1], [xcz, xcore])
        mnew_slope = mass[i1:icore_new]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        mnew_slope = np.array(mnew_slope).astype(float)

        # h1_slope_new = xsurface + mnew_slope * (xcz - xsurface)
        h1_slope_new2 = fslope(mnew_slope)  

        h1_slope_new_mid = np.full(i1 - i2, xcz)

        h1_slope_new = np.concatenate((h1_slope_new1, h1_slope_new_mid, h1_slope_new2))

    h1_new = np.concatenate((h1_new_surface, h1_slope_new, h1_new_core))



    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        if isotope == 'h1':
            return h1_new
        elif isotope == 'he4':
            return iso_profiles['he4'] + (iso_profiles['h1'] - h1_new)

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d([0, 1], [prof_slope[0], prof_slope[1]])

        mnew_outer = mass[:ienv_new]
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = mass[ienv_new:icore]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        
        mnew_outer = np.array(mnew_outer).astype(float)
        mnew_slope = np.array(mnew_slope).astype(float)
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)
        prof_slope_new = prof_slope[0] + mnew_slope * (prof_slope[-1] - prof_slope[0])

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds




def mod_ms_7param(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    Modifies slope of hydrogen profile of a MS star
    should maintain the same core and surface abundances and core mass
    where core is the region of flat abundance in the center
    """

    # mslope is the mass over which the slope changes
    if 'mcore' in kwargs:
        mcore = kwargs['mcore']
    else:
        # Get the current convective core mass ratio
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.01

        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mcore_old = mass[icore]/mass[0]

        mcore = mcore_old

    if 'mslope' in kwargs:
        mslope = kwargs['mslope']
    else:
        # Get the current mass over which slope changes
        h_prof = iso_profiles['h1']
        h_prof_rev = h_prof[::-1]
        h_core_diff = 0.02

        ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

        h_core_diff = 0.002
        icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
        icore = len(h_prof) - icore_rev - 1

        mslope_old = mass[ienv] - mass[icore]
        mslope = mslope_old

    if 'xcore' in kwargs:
        xcore = kwargs['xcore']
    else:
        xcore = iso_profiles['h1'][-1]

    if 'xsurface' in kwargs:
        xsurface = kwargs['xsurface']
    else:
        xsurface = iso_profiles['h1'][0]

    if 'fa' in kwargs:
        fa = kwargs['fa']
        ma = mcore + fa * (mslope)
    else:
        ma = 0

    if 'fb' in kwargs:
        fb = kwargs['fb']
        mb = mcore + mslope + fb * (1 - mslope - mcore)
    else:
        mb = 0

    if 'alpha' in kwargs:
        alpha = kwargs['alpha']
    else:
        alpha = 0


    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))

    
    # Compute new H1 profile here:

    icore_new = next((x[0] for x in enumerate(mass) if x[1] < mcore * mass[0]))
    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mslope) * mass[0]))

    h1_new_core = np.full(len(mass) - icore_new, xcore) 
    h1_new_surface = np.full(ienv_new, xsurface)

    fslope = interpolate.interp1d([0, 1], [xsurface, xcore])
    mnew_slope = mass[ienv_new:icore_new]
    mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
    mnew_slope = np.array(mnew_slope).astype(float)

    h1_slope_new = xsurface + mnew_slope * (xcore - xsurface)
    h1_slope_new = fslope(mnew_slope)

    h1_new = np.concatenate((h1_new_surface, h1_slope_new, h1_new_core))

    # ival = len(h1_new_surface)

    if ma != 0 and mb != 0 and alpha > 0:

        xvals = mass[::-1]
        yvals = h1_new[::-1]

        ia = next((x[0] for x in enumerate(xvals) if x[1] > ma * xvals[-1]))
        ib = next((x[0] for x in enumerate(xvals) if x[1] > mb * xvals[-1]))
        ival = next((x[0] for x in enumerate(xvals) if x[1] > (mcore + mslope) * xvals[-1]))

        xsect = np.array(xvals[ia:ib+1], dtype="float64")
        ysect = np.array(yvals[ia:ib+1], dtype="float64")

        xsect_norm = (xsect - xsect[0])/(xsect[-1] - xsect[0])

        x = np.linspace(0, 1, 10000)
        y = alpha/(x + alpha) - (alpha/(1 + alpha))*x

        theta1 = np.arctan((yvals[ib] - yvals[ival])/(xvals[ib] - xvals[ival]))
        theta2 = np.pi + np.arctan((yvals[ival] - yvals[ia])/(xvals[ival] - xvals[ia]))

        l2 = np.sqrt((yvals[ia] - yvals[ival])**2 + (xvals[ia] - xvals[ival])**2)
        l1 = np.sqrt((yvals[ib] - yvals[ival])**2 + (xvals[ib] - xvals[ival])**2)

        rot = np.array([[l1*np.cos(theta1), l2*np.cos(theta2)], [l1*np.sin(theta1), l2*np.sin(theta2)]])

        xnew, ynew = np.dot(rot, np.column_stack((x, y)).T)
        xnew = xnew + xvals[ival]
        ynew = ynew + yvals[ival]

        f_mass_ynew = interpolate.interp1d(xnew, ynew, fill_value='extrapolate')

        h1_curve = f_mass_ynew(xsect)

        h1_new[::-1][ia:ib+1] = h1_curve

 
    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        if isotope == 'h1':
            return h1_new
        elif isotope == 'he4':
            return iso_profiles['he4'] + (iso_profiles['h1'] - h1_new)

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d([0, 1], [prof_slope[0], prof_slope[1]])

        mnew_outer = mass[:ienv_new]
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = mass[ienv_new:icore]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        
        mnew_outer = np.array(mnew_outer).astype(float)
        mnew_slope = np.array(mnew_slope).astype(float)
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)
        prof_slope_new = prof_slope[0] + mnew_slope * (prof_slope[-1] - prof_slope[0])

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds



def mod_ms_7param_v2(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    mabund is the normalised mass coordinate at the edge of the 
    sloped H profile region
    """

    mabund = kwargs['mabund']
    slope_ratio = kwargs['slope_ratio']
    xcore = kwargs['xcore']
    xsurface = kwargs['xsurface']

    xcore = np.maximum(xcore, 0)

    mslope = (xsurface - xcore)/slope_ratio
    mcore = mabund - mslope


    fa = kwargs['fa']
    ma = mcore + fa * (mslope)

    fb = kwargs['fb']
    mb = mcore + mslope + fb * (1 - mslope - mcore)

    alpha = kwargs['alpha']**2

    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))

    
    # Compute new H1 profile here:

    icore_new = next((x[0] for x in enumerate(mass) if x[1] < mcore * mass[0]))
    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mslope) * mass[0]))

    h1_new_core = np.full(len(mass) - icore_new, xcore) 
    h1_new_surface = np.full(ienv_new, xsurface)

    fslope = interpolate.interp1d([0, 1], [xsurface, xcore])
    mnew_slope = mass[ienv_new:icore_new]
    mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
    mnew_slope = np.array(mnew_slope).astype(float)

    h1_slope_new = xsurface + mnew_slope * (xcore - xsurface)
    h1_slope_new = fslope(mnew_slope)

    h1_new = np.concatenate((h1_new_surface, h1_slope_new, h1_new_core))

    # ival = len(h1_new_surface)

    if ma != 0 and mb != 0 and alpha > 0:

        xvals = mass[::-1]
        yvals = h1_new[::-1]

        ia = next((x[0] for x in enumerate(xvals) if x[1] > ma * xvals[-1]))
        ib = next((x[0] for x in enumerate(xvals) if x[1] > mb * xvals[-1]))
        ival = next((x[0] for x in enumerate(xvals) if x[1] > (mcore + mslope) * xvals[-1]))

        xsect = np.array(xvals[ia:ib+1], dtype="float64")
        ysect = np.array(yvals[ia:ib+1], dtype="float64")

        xsect_norm = (xsect - xsect[0])/(xsect[-1] - xsect[0])

        x = np.linspace(0, 1, 10000)
        y = alpha/(x + alpha) - (alpha/(1 + alpha))*x

        theta1 = np.arctan((yvals[ib] - yvals[ival])/(xvals[ib] - xvals[ival]))
        theta2 = np.pi + np.arctan((yvals[ival] - yvals[ia])/(xvals[ival] - xvals[ia]))

        l2 = np.sqrt((yvals[ia] - yvals[ival])**2 + (xvals[ia] - xvals[ival])**2)
        l1 = np.sqrt((yvals[ib] - yvals[ival])**2 + (xvals[ib] - xvals[ival])**2)

        rot = np.array([[l1*np.cos(theta1), l2*np.cos(theta2)], [l1*np.sin(theta1), l2*np.sin(theta2)]])

        xnew, ynew = np.dot(rot, np.column_stack((x, y)).T)
        xnew = xnew + xvals[ival]
        ynew = ynew + yvals[ival]

        f_mass_ynew = interpolate.interp1d(xnew, ynew, fill_value='extrapolate')

        h1_curve = f_mass_ynew(xsect)

        h1_new[::-1][ia:ib+1] = h1_curve

 
    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        if isotope == 'h1':
            return h1_new
        elif isotope == 'he4':
            return iso_profiles['he4'] + (iso_profiles['h1'] - h1_new)

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d([0, 1], [prof_slope[0], prof_slope[1]])

        mnew_outer = mass[:ienv_new]
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = mass[ienv_new:icore]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        
        mnew_outer = np.array(mnew_outer).astype(float)
        mnew_slope = np.array(mnew_slope).astype(float)
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)
        prof_slope_new = prof_slope[0] + mnew_slope * (prof_slope[-1] - prof_slope[0])

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds



def mod_ms_7param_evol(iso_profiles, isotopes, mass, dm, starmass, **kwargs):
    """
    mabund is the normalised mass coordinate at the edge of the 
    sloped H profile region
    """

    mabund = kwargs['mabund']
    slope_ratio = kwargs['slope_ratio']
    xcore = kwargs['xcore']
    xsurface = kwargs['xsurface']

    xcore = np.maximum(xcore, 0)

    mslope = (xsurface - xcore)/slope_ratio
    mcore = mabund - mslope


    fa = kwargs['fa']
    ma = mcore + fa * (mslope)

    fb = kwargs['fb']
    mb = mcore + mslope + fb * (1 - mslope - mcore)

    alpha = kwargs['alpha']**2

    # Get the current mass over which slope changes
    h_prof = iso_profiles['h1']
    h_prof_rev = h_prof[::-1]
    h_core_diff = 0.02

    ienv = next((x[0] for x in enumerate(h_prof) if np.abs(x[1] - h_prof[0]) > h_core_diff))

    h_core_diff = 0.002
    icore_rev = next((x[0] for x in enumerate(h_prof_rev) if x[1] - h_prof_rev[0] > h_core_diff))
    icore = len(h_prof) - icore_rev - 1

    mslope_old = mass[ienv] - mass[icore]

    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < mass[icore] + mslope))

    
    # Compute new H1 profile here:

    icore_new = next((x[0] for x in enumerate(mass) if x[1] < mcore * mass[0]))
    ienv_new = next((x[0] for x in enumerate(mass) if x[1] < (mcore + mslope) * mass[0]))

    h1_new_core = np.full(len(mass) - icore_new, xcore) 
    h1_new_surface = np.full(ienv_new, xsurface)

    fslope = interpolate.interp1d([0, 1], [xsurface, xcore])
    mnew_slope = mass[ienv_new:icore_new]
    mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
    mnew_slope = np.array(mnew_slope).astype(float)

    h1_slope_new = xsurface + mnew_slope * (xcore - xsurface)
    h1_slope_new = fslope(mnew_slope)

    h1_new = np.concatenate((h1_new_surface, h1_slope_new, h1_new_core))

    # ival = len(h1_new_surface)

    if ma != 0 and mb != 0 and alpha > 0:

        xvals = mass[::-1]
        yvals = h1_new[::-1]

        ia = next((x[0] for x in enumerate(xvals) if x[1] > ma * xvals[-1]))
        ib = next((x[0] for x in enumerate(xvals) if x[1] > mb * xvals[-1]))
        ival = next((x[0] for x in enumerate(xvals) if x[1] > (mcore + mslope) * xvals[-1]))

        xsect = np.array(xvals[ia:ib+1], dtype="float64")
        ysect = np.array(yvals[ia:ib+1], dtype="float64")

        xsect_norm = (xsect - xsect[0])/(xsect[-1] - xsect[0])

        x = np.linspace(0, 1, 10000)
        y = alpha/(x + alpha) - (alpha/(1 + alpha))*x

        theta1 = np.arctan((yvals[ib] - yvals[ival])/(xvals[ib] - xvals[ival]))
        theta2 = np.pi + np.arctan((yvals[ival] - yvals[ia])/(xvals[ival] - xvals[ia]))

        l2 = np.sqrt((yvals[ia] - yvals[ival])**2 + (xvals[ia] - xvals[ival])**2)
        l1 = np.sqrt((yvals[ib] - yvals[ival])**2 + (xvals[ib] - xvals[ival])**2)

        rot = np.array([[l1*np.cos(theta1), l2*np.cos(theta2)], [l1*np.sin(theta1), l2*np.sin(theta2)]])

        xnew, ynew = np.dot(rot, np.column_stack((x, y)).T)
        xnew = xnew + xvals[ival]
        ynew = ynew + yvals[ival]

        f_mass_ynew = interpolate.interp1d(xnew, ynew, fill_value='extrapolate')

        h1_curve = f_mass_ynew(xsect)

        h1_new[::-1][ia:ib+1] = h1_curve

 
    def get_abund_isotope(isotope):
        """
        returns isotopes abundance for new structure
        """

        if isotope == 'h1':
            return h1_new
        elif isotope == 'he4':
            return iso_profiles['he4'] + (iso_profiles['h1'] - h1_new)

        iso_prof = iso_profiles[isotope]

        prof_outer = iso_prof[:ienv]
        mass_outer = mass[:ienv]

        prof_slope = iso_prof[ienv:icore]
        mass_slope = mass[ienv:icore]

        mass_fit_outer = (mass_outer - mass_outer[0])/(mass_outer[-1] - mass_outer[0])
        mass_fit_slope = (mass_slope - mass_slope[0])/(mass_slope[-1] - mass_slope[0])

        fout = interpolate.interp1d(mass_fit_outer, prof_outer)
        fslope = interpolate.interp1d([0, 1], [prof_slope[0], prof_slope[1]])

        mnew_outer = mass[:ienv_new]
        mnew_outer = (mnew_outer - mnew_outer[0])/(mnew_outer[-1] - mnew_outer[0])

        mnew_slope = mass[ienv_new:icore]
        mnew_slope = (mnew_slope - mnew_slope[0])/(mnew_slope[-1] - mnew_slope[0])
        
        mnew_outer = np.array(mnew_outer).astype(float)
        mnew_slope = np.array(mnew_slope).astype(float)
        
        prof_outer_new = fout(mnew_outer)
        prof_slope_new = fslope(mnew_slope)
        prof_slope_new = prof_slope[0] + mnew_slope * (prof_slope[-1] - prof_slope[0])

        prof_inner = iso_prof[icore:]

        prof_new = np.concatenate((prof_outer_new, prof_slope_new, prof_inner))

        return prof_new


    newabunds = [get_abund_isotope(iso) for iso in isotopes]

    return newabunds


# list of different ways to modify models

functions = {
    'mod_grad_const_hemass': mod_grad_const_hemass,
    'mod_h_grad': mod_h_grad,
    'mod_h_grad_absolute': mod_h_grad_absolute,
    'mod_draw_grad': mod_draw_grad,
    'mod_change_z': mod_change_z,
    'mod_change_core_abund': mod_change_core_abund,
    'mod_env_const_hemass': mod_env_const_hemass,
    'mod_env_const_grad': mod_env_const_grad,
    'mod_grad_absolute_flex_surface': mod_grad_absolute_flex_surface,
    'mod_change_z_outonly': mod_change_z_outonly,
    # 'mod_ms_z_fraction': mod_ms_z_fraction,
    'mod_flexsurf_arb': mod_flexsurf_arb,
    'mod_cno_hshell': mod_cno_hshell,
    'mod_cno_envelope': mod_cno_envelope,
    'mod_metals_envelope': mod_metals_envelope,
    'mod_add_arb2': mod_add_arb2,
    'mod_add_arb1': mod_add_arb1,
    'mod_add_h1': mod_add_h1,
    'mod_ms_cno_core': mod_ms_cno_core,
    'mod_metals_envelope_tohe': mod_metals_envelope_tohe,
    'mod_grad_const_hemass_popIII': mod_grad_const_hemass_popIII,
    'mod_ms_metals_envelope': mod_ms_metals_envelope,
    'add_arb_to_solar': add_arb_to_solar,
    'mod_fuel_supply': mod_fuel_supply,
    'mod_mucore': mod_mucore,
    'mod_consthe_homogen': mod_consthe_homogen,
    'mod_mass_gainer_heburning': mod_mass_gainer_heburning,
    'mod_coreratio': mod_coreratio,
    'mod_mcore_ms': mod_mcore_ms,
    'mod_slope_ms': mod_slope_ms,
    'mod_straight_slope_ms': mod_straight_slope_ms,
    'mod_3part_ms': mod_3part_ms,
    'mod_same': mod_same,
    'mod_3part_ms_curved': mod_3part_ms_curved,
    'mod_3part_ms_curved_bot': mod_3part_ms_curved_bot,
    'mod_3part_ms_icz': mod_3part_ms_icz,
    'mod_ms_7param': mod_ms_7param,
    'mod_ms_7param_v2': mod_ms_7param_v2,
    'mod_ms_7param_evol': mod_ms_7param_evol
}


if __name__ == '__main__':
    pass
