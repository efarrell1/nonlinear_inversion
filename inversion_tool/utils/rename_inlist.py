
def edit_inlist(template_inlist, dict_changes, modified_inlist, quotes=[]):
    '''
    template_inlist: location of template inlist
    dict_chages: dicionary containing parameters as keys and values as values
    modified_inlist: what to save modified inlist as
    quotes: optional, puts quotes around value if parameter is in this list

    -- searches through template inlist and replaces values with
    -- those from dict_changes

    Output:
    saves new, modified inlist as modified_inlist
    '''

    # reads input file and saves lines as list
    inputFile = open(template_inlist, "r+")
    inFile = inputFile.readlines()
    inputFile.close()

    infile_split = [line.split(' = ') for line in inFile]
    mesa_params = [l[0] for l in infile_split]
    mesa_params = [x.strip() for x in mesa_params]


    for mesa_param, value in dict_changes.items():
        if mesa_param in quotes:
            value = "'" + value + "'"
        try:
            i = mesa_params.index(mesa_param)
            infile_split[i][1] = str(value) + '\n'
        except ValueError:
            print(mesa_param, 'not found in inlist')

    separator = ' = '
    new_inlist = [separator.join(l) for l in infile_split]

    # opens output file and writes to file
    outFile = open(modified_inlist, "w")
    outFile.writelines(new_inlist)


def add_to_inlist(template_inlist, sections, dict_changes, modified_inlist):
    '''
    template_inlist: location of template inlist
    sections: dictionary with new parameters as keys and part of inlist that
    new parameters (star_job, controls or pgstar) correspond to
    as values
    dict_changes: dicionary containing new parameters as keys and values as
    values
    modified_inlist: what to save modified inlist as

    -- searches through template inlist and adds new keys and values
    -- those from dict_changes
    -- also adds comment saying variables have been modified

    Output:
    saves new, modified inlist as modified_inlist
    '''

    # reads input file and saves lines as list
    inputFile = open(template_inlist, "r+")
    inFile = inputFile.readlines()
    inputFile.close()

    infile_split = [line.split(' = ') for line in inFile]
    mesa_params = [l[0] for l in infile_split]

    # Adds in notice that the variables have been added to the inlist
    for s in set(sections.values()):
        try:
            section = '&' + s + '\n'
            i = mesa_params.index(section)
            infile_split.insert(i+1, ['! ----- ADDED SETTINGS ----- !\n'])
        except ValueError:
            None

    for mesa_param, value in dict_changes.items():
        try:
            section = '&' + sections[mesa_param] + '\n'
            i = mesa_params.index(section)
            infile_split.insert(i+2, [mesa_param, str(value) + '\n'])
        except ValueError:
            print(mesa_param, 'not included in inlist')
        except KeyError:
            print('no section provided for ', mesa_param)

    separator = ' = '
    new_inlist = [separator.join(l) for l in infile_split]

    # opens output file and writes to file
    outFile = open(modified_inlist, "w")
    outFile.writelines(new_inlist)
