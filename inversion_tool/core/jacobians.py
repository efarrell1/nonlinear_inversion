
jacob_param_dict_1 = {'rot': ['slope', 'intercept', 'P0', 'P1', 'P2', 'P3'],
                       'mass': ['P0', 'P1', 'P2', 'P3', 'slope', 'intercept'],
                       'alpha': ['A', 'slope', 'intercept','P0', 'P1', 'P2', 'P3'],
                       'xcore': ['slope', 'intercept', 'P0', 'P1', 'P2', 'P3', 'P4']}

jacob_param_dict_2 = {'rot': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10'],
                       'mass': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10'],
                       'alpha': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10'],
                       'xcore': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10']}

jacob_param_dict_22 = {'rot': ['dP0', 'dP1', 'dP2', 'dP3', 'dP4', 'dP5', 'dP6', 'dP7', 'dP8', 'dP9', 'dP10'],
                       'mass': ['dP0', 'dP1', 'dP2', 'dP3', 'dP4', 'dP5', 'dP6', 'dP7', 'dP8', 'dP9', 'dP10'],
                       'alpha': ['dP0', 'dP1', 'dP2', 'dP3', 'dP4', 'dP5', 'dP6', 'dP7', 'dP8', 'dP9', 'dP10'],
                       'xcore': ['dP0', 'dP1', 'dP2', 'dP3', 'dP4', 'dP5', 'dP6', 'dP7', 'dP8', 'dP9', 'dP10']}

jacob_param_dict_3 = {'rot': ['slope', 'intercept', 'P0', 'PlastP0'],
                       'mass': ['P0', 'slope', 'intercept'],
                       'alpha': ['A', 'slope', 'intercept', 'P0'],
                       'xcore': ['slope', 'A', 'intercept', 'P0', 'PlastP0']}

jacob_param_dict_4 = {'rot': ['slope', 'intercept', 'P0'],
                       'mass': ['P0', 'slope', 'intercept'],
                       'alpha': ['A', 'slope', 'intercept'],
                       'xcore': ['slope', 'intercept', 'dP1', 'dP2', 'dP3', 'dP4', 'dP5', 'dP6', 'dP7', 'dP8']}

jacob_param_dict_5 = {'rot': ['slope', 'intercept', 'P0', 'dP0', 'avg_ddP'],
                       'mass': ['slope', 'intercept', 'P0', 'dP0', 'avg_ddP'],
                       'alpha': ['A'],
                       'xcore': ['slope', 'intercept', 'P0', 'dP0', 'avg_ddP']}

jacob_param_dict_6 = {'rot': ['slope', 'intercept', 'P0', 'avg_ddP'],
                       'mass': ['slope', 'intercept', 'P0', 'avg_ddP'],
                       'alpha': ['A'],
                       'xcore': ['slope', 'intercept', 'P0', 'avg_ddP']}

jacob_param_dict_7 = {'rot': ['slope', 'intercept', 'P0', 'avg_ddP'],
                       'mass': ['slope', 'intercept', 'P0', 'avg_ddP'],
                       'xcore': ['slope', 'intercept', 'P0', 'avg_ddP']}

jacob_param_dict_8 = {'rot': ['slope', 'intercept', 'P0'],
                       'mass': ['slope', 'intercept', 'P0'],
                       'xcore': ['slope', 'intercept', 'P0']}

jacob_param_dict_9 = {'rot': ['slope', 'intercept', 'P0', 'dP4', 'dP7'],
                       'mass': ['slope', 'intercept', 'P0', 'dP4', 'dP7'],
                       'xcore': ['slope', 'intercept', 'P0', 'dP4', 'dP7']}

jacob_param_dict_10 = {'rot': ['slope', 'intercept', 'P0'],
                       'mass': ['slope', 'intercept', 'P0'],
                       'xcore': ['slope', 'intercept', 'P0', 'T']}

jacob_param_dict_11 = {'rot': ['slope', 'intercept', 'P0'],
                       'mass': ['slope', 'intercept', 'P0'],
                       'xcore': ['slope', 'intercept', 'nmin0', 'nmin1']}

jacob_param_dict_12 = {'rot': ['slope', 'intercept', 'P0'],
                       'mass': ['slope', 'intercept', 'P0'],
                       'xcore': ['slope', 'intercept', 'nmin0', 'nmin1']}

jacob_param_dict_13 = {'rot': ['slope', 'intercept', 'ddP0'],
                       'mass': ['slope', 'intercept', 'ddP0'],
                       'xcore': ['slope', 'intercept', 'ddP0']}

jacob_param_dict_14 = {'rot': ['slope', 'intercept', 'P0'],
                       'mass': ['slope', 'intercept', 'P0'],
                       'xcore': ['slope', 'intercept', 'P0']}

# jacob_param_dict_15 = {'rot': ['slope', 'intercept', 'P0', 'P1', 'P2', 'P3', 'P4', 'P7'],
#                        'mass': ['slope', 'intercept', 'P0', 'P1', 'P2', 'P3', 'P4', 'P7'],
#                        'xcore': ['slope', 'intercept', 'P0', 'P1', 'P2', 'P3', 'P4', 'P7']}

jacob_param_dict_15 = {'rot': ['P0', 'dP0', 'Plast'],
                       'mass': ['P0', 'dP0', 'Plast'],
                       'alpha': ['A'],
                       'xcore': ['P0', 'dP0', 'Plast']}

# jacob_param_dict_15b = {'rot': ['P0', 'dP0', 'Plast', 'nmin0', 'nmin1'],
#                        'mass': ['P0', 'Plast', 'nmin0', 'nmin1'],
#                        'alpha': [],
#                        'xcore': ['P0', 'nmin0', 'nmin1']}

# jacob_param_dict_15b = {'rot': ['P0', 'dP0', 'Plast'],
#                        'mass': ['P0', 'Plast'],
#                        'alpha': [],
#                        'xcore': ['nmin0', 'nmin1', 'nmin2']}

jacob_param_dict_15b = {'rot': ['P0', 'dP0', 'Plast'],
                       'mass': ['P0', 'dP0', 'Plast'],
                       'xcore': ['P0', 'dP0', 'Plast', 'nmin0', 'nmin1', 'nmin2', 'nmin3']}


# jacob_param_dict_15 = {'rot': ['P0', 'Plast'],
#                        'mass': ['P0', 'Plast', 'nmin0', 'nmin1'],
#                        'alpha': [],
#                        'xcore': ['P0', 'Plast', 'nmin0', 'nmin1']}



jacob_param_dict_20 = {'rot': ['P0', 'Plast'],
                       'mass': ['P0', 'Plast'],
                       'alpha': [],
                       'xcore': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2']}

# jacob_param_dict_20 = {'rot': ['P0', 'Plast'],
#                        'mass': ['P0', 'Plast'],
#                        'alpha': [],
#                        'xcore': ['P0', 'Plast', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13']}


jacob_param_dict_16 = {'rot': ['P0', 'dP0', 'slope', 'Plast'],
                       'mass': ['P0', 'Plast'],
                       'alpha': ['A'],
                       'xcore': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2']}

jacob_param_dict_17 = {'rot': ['P0', 'dP0', 'slope'],
                       'mass': ['P0', 'slope'],
                       'alpha': ['A'],
                       'xcore': ['P0', 'nmin1', 'nmin2', 'nmin3']}

jacob_param_dict_18 = {'rot': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2'],
                       'mass': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2'],
                       'alpha': ['A'],
                       'xcore': ['nmin0', 'nmin1', 'nmin2'],
                       'mabund': []}


jacob_param_dict_18 = {'rot': ['P0', 'Plast', 'dn0', 'dn1'],
                       'mass': ['P0', 'Plast'],
                       'xcore': ['P0', 'dn0', 'dn1']}


jacob_param_dict_18 = {'rot': ['P0', 'Plast'],
                       'mass': ['P0', 'Plast'],
                       'xcore': ['P0', 'dn0', 'dn1', 'dn2']}


jacob_param_dict_18b = {'rot': ['P0', 'Plast'],
                        'mass': ['P0', 'Plast'],
                        'xcore': ['P0', 'dn0', 'dn1']}


jacob_param_dict_18c = {'rot': ['P0', 'Plast',  'nmin0', 'nmin1'],
                       'mass': ['P0', 'Plast', 'nmin0', 'nmin1'],
                       'alpha': [],
                       'xcore': ['P0', 'Plast', 'nmin0', 'nmin1'],
                       'mabund': []}

jacob_param_dict_18d = {'rot': ['P0', 'Plast',  'nmin0', 'nmin1', 'nmin2', 'nmin3'],
                       'mass': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2', 'nmin3'],
                       'alpha': ['A'],
                       'xcore': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2', 'nmin3'],
                       'mabund': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'nmin0', 'nmin1', 'nmin2', 'nmin3']}

jacob_param_dict_18d = {'rot': ['P0', 'Plast',  'nmin0', 'nmin1', 'nmin2', 'nmin3'],
                       'mass': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2', 'nmin3'],
                       'alpha': ['A'],
                       'xcore': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2', 'nmin3'],
                       'mabund': ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'nmin0', 'nmin1', 'nmin2', 'nmin3']}

jacob_param_dict_18d = {'rot': ['P0', 'Plast',  'nmin0', 'nmin1', 'nmin2'],
                       'mass': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2'],
                       'alpha': ['A'],
                       'xcore': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2'],
                       'mabund': ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'nmin0', 'nmin1', 'nmin2']}

jacob_param_dict_18d = {'rot': ['P0', 'Plast',  'nmin0', 'nmin1', 'nmin2'],
                       'mass': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2'],
                       'alpha': [],
                       'xcore': ['P0', 'Plast', 'nmin0', 'nmin1', 'nmin2'],
                       'mabund': []}

jacob_param_dict_18d = {'rot': ['P0', 'P16'],
                       'mass': ['P0', 'P16', 'A'],
                       'alpha': ['A'],
                       'xcore': [ 'nmin0', 'nmin1', 'nmin2'],
                       'mabund': ['P3', 'P4',  'P5', 'P6', 'P7', 'P8', 'P9', 'P16', 'nmin0', 'nmin1', 'nmin2']} # , 'nmin0', 'nmin1', 'nmin2'

jacob_param_dict_805 = {'rot': ['P0', 'P11'],
                       'mass': ['P0', 'P11'],
                       'alpha': ['A'],
                       'xcore': [ 'nmin0', 'nmin1'],
                       'mabund': ['P0', 'P1',  'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'nmin0', 'nmin1']}


jacob_param_dict_18e = {'rot': ['P3', 'P16'],
                       'mass': ['P3', 'P16'],
                       'alpha': ['A'],
                       'xcore': ['P3', 'P16', 'dn0', 'dn1' ],
                       'mabund': ['P3', 'P4',  'P5', 'P6', 'P7', 'P8', 'P9', 'P16', 'dn0', 'dn1']}



jacob_param_dict_19 = {'rot': ['P0', 'Plast'],
                       'mass': ['P0', 'nmin0', 'nmin1', 'nmin2'],
                       'mabund': ['P0'],
                       'xcore': ['P0', 'nmin0', 'nmin1', 'nmin2']}



jacob_param_dict_19b = {'rot': ['P0', 'P8'],
                        'mass': ['P0', 'P8'],
                        'alpha': ['A'],
                        'xcore': ['nmin0', 'nmin1']}


jacob_param_dict_19c = {'rot': ['P0', 'P8'],
                        'mass': ['P0', 'P8', 'jshape1', 'nmin0', 'nmin1'],
                        'alpha': ['A'],
                        'xcore': ['jshape1', 'nmin0', 'nmin1'],
                        'mabund': ['jshape1', 'nmin0', 'nmin1']}





jacob_param_dict_21 = {'mabund': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'Plast'],
                       'slope_ratio': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'Plast'],
                       'xcore': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'Plast'],
                       'rot': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'Plast'],
                       'mass': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'Plast'],
                       'alpha': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'Plast']}

jacob_param_dict_21b = {'mabund': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17'],
                       'slope_ratio': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17'],
                       'xcore': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17'],
                       'rot': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17'],
                       'mass': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17'],
                       'alpha': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17']}

jacob_param_dict_21c = {'mabund': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'],
                       'slope_ratio': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'],
                       'xcore': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'],
                       'rot': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'],
                       'mass': ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'],
                       'alpha': []}


jacob_param_dict_23 = {'mabund': ['dP0', 'dP1', 'dP2', 'dP3', 'dP4', 'dP5', 'dP6', 'dP7', 'dP8', 'dP9', 'dP10'],
                       'slope_ratio': ['dP0', 'dP1', 'dP2', 'dP3', 'dP4', 'dP5', 'dP6', 'dP7', 'dP8', 'dP9', 'dP10'],
                       'xcore': ['dP0', 'dP1', 'dP2', 'dP3', 'dP4', 'dP5', 'dP6', 'dP7', 'dP8', 'dP9', 'dP10']}


jacob_param_dict_shape = {'rot': ['P0', 'P8'],
                       'mass': ['P0', 'P8'],
                       'alpha': ['A'],
                       'xcore': [ 'jshape1', 'jshape2', 'jshape3', 'jshape4', 'jshape5', 'jshape6', 'jshape7', 'jshape8', 'jshape9', 'jshape10'],
                       'mabund': ['jshape1', 'jshape2', 'jshape3', 'jshape4', 'jshape5', 'jshape6', 'jshape7', 'jshape8', 'jshape9', 'jshape10']} # , 'nmin0', 'nmin1', 'nmin2'

jacob_param_slope = {'slope_ratio': ['nmin0', 'nmin1', 'nmin2']}

jacob_param_test = {'xcore': ['nmin0', 'nmin1', 'nmin2'],
                    'slope_ratio': ['dP2', 'dP6', 'dP12']}


jacob_param_options = {1: jacob_param_dict_1,
                       2: jacob_param_dict_2,
                       3: jacob_param_dict_3,
                       4: jacob_param_dict_4,
                       5: jacob_param_dict_5,
                       6: jacob_param_dict_6,
                       7: jacob_param_dict_7,
                       8: jacob_param_dict_8,
                       9: jacob_param_dict_9,
                       10: jacob_param_dict_10,
                       11: jacob_param_dict_11,
                       12: jacob_param_dict_12,
                       13: jacob_param_dict_13,
                       14: jacob_param_dict_14,
                       15: jacob_param_dict_15,
                       '15b': jacob_param_dict_15b,
                       16: jacob_param_dict_16,
                       17: jacob_param_dict_17,
                       18: jacob_param_dict_18,
                       '18b': jacob_param_dict_18b,
                       '18c': jacob_param_dict_18c,
                       '18d': jacob_param_dict_18d,
                       '18e': jacob_param_dict_18e,
                       19: jacob_param_dict_19,
                       '19b': jacob_param_dict_19b,
                       '19c': jacob_param_dict_19c,
                       20: jacob_param_dict_20,
                       21: jacob_param_dict_21,
                       22: jacob_param_dict_22,
                       23: jacob_param_dict_23,
                       '21b': jacob_param_dict_21b,
                       18805: jacob_param_dict_805,
                       'shape': jacob_param_dict_shape,
                       'slope': jacob_param_slope,
                       'test': jacob_param_test,
                       '21c': jacob_param_dict_21c}

