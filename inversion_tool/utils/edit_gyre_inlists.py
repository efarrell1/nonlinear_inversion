#!/Users/eoin/opt/anaconda3/bin/python

'''
Modifies GYRE inlists
'''

import sys
from rename_inlist import edit_inlist, add_to_inlist

# initial ZAMS mass and subdirectory
inlist_name = sys.argv[1]
try:
	rot_val = sys.argv[2]
except IndexError:
	rot_val = 0



dict_changes = {
	'coriolis_method': 'TAR',
	'Omega_rot': str(float(rot_val)/1000),
}

quotes = ['coriolis_method']

print('Creating GYRE Inlists for ' + inlist_name, 'rot', dict_changes['Omega_rot'])

edit_inlist(inlist_name, dict_changes, inlist_name, quotes=quotes)

