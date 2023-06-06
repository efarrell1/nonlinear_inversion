#!/Users/eoin/opt/anaconda3/bin/python

'''
Modifies paused and relaxed inlists
'''

import sys
from .rename_inlist import edit_inlist, add_to_inlist

# initial ZAMS mass and subdirectory
modname = sys.argv[1]
subd = sys.argv[2]
tpl = sys.argv[3]


print('Creating Inlists for ' + modname)

dotmod_name = "'" + modname + ".mod'"

template_inlist = tpl + '/paused_mcmc/inlist_paused'
modified_inlist = subd + '/' + modname + '/' + 'inlist_paused'
dict_changes = {
	'load_model_filename': dotmod_name
}

edit_inlist(template_inlist, dict_changes, modified_inlist)

