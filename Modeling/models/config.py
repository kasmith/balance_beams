import os

######
#
# Read in some basic configuration information
#
######

cfg_fl = os.path.join(os.path.dirname(__file__), "config.txt")
cf = open(cfg_fl, 'rU')
cfg_dict = dict()
for ln in cf:
    pn, tp = ln.strip('\n').split(':')
    pn = pn.strip(' ')
    tp = tp.strip(' ')
    if tp in ['T', 'True', 'TRUE', 'true']:
        cfg_dict[pn] = True
    elif tp in ['F', 'False', 'FALSE', 'false']:
        cfg_dict[pn] = False
    else:
        try:
            tp = int(tp)
            cfg_dict[pn] = tp
        except:
            try:
                tp = float(tp)
                cfg_dict[pn] = tp
            except:
                cfg_dict[pn] = tp
cf.close()

#######
#
# Set up the parameter names that are used by each of the model types
#
#######

# Removes objects from the parameter list depending on model type


def safe_remove(lst, objs):
    if not hasattr(objs, '__iter__'):
        objs = [objs]
    for obj in objs:
        try:
            lst.remove(obj)
        except ValueError:
            pass


common_param_names = ['mass_jnd', 'dist_jnd', 'com_u',
                      'block_mass_u', 'distance_u']
clean_param_names = common_param_names + ['com_range']
shapes_param_names = common_param_names + ['tower_mass_u', 'com_range']
materials_param_names = common_param_names + \
    ['brick_density', 'brick_u', 'iron_density', 'iron_u', 'com_range']
balance_param_names = common_param_names + ['com_range_vs', 'com_range_s', 'com_range_m', 'com_range_l',
                                            'beam_mass_mean', 'beam_mass_sd']
combined_param_names = common_param_names + ['tower_mass_u', 'com_range_vs', 'com_range_m', 'beam_mass_mean',
                                             'beam_mass_sd', 'brick_density', 'brick_u', 'iron_density', 'iron_u']


def params_mod_2_rules(params, doadd=True):
    safe_remove(params, ['mass_heur_p', 'dist_heur_p', 'com_u'])
    if doadd:
        params += ['use_rule1', 'use_rule2', 'use_rule3', 'use_rule3a']
    return params


def read_def_params(flnm):
    fl = open(flnm, 'rU')
    l = next(fl).split(',')
    if not len(l) == 5:
        raise Exception('Not properly formatted parameter file')
    params = dict()
    for ln in fl:
        l = ln.strip('\n').split(',')
        nm = l[0]
        if l[1] == 'None':
            dflt = None
        else:
            dflt = float(l[1])
        params[nm] = (dflt, float(l[2]), (float(l[3]), float(l[4])))
    return params


default_params = read_def_params(os.path.join(os.path.dirname(__file__), "default_params.csv"))


def fill_params(use_param_dict, param_names, def_params=None):
    fill_dict = dict()
    for n in param_names:
        if n in use_param_dict.keys():
            fill_dict[n] = use_param_dict[n]
        else:
            if def_params is None:
                fill_dict[n] = default_params[n][0]
            else:
                fill_dict[n] = def_params[n]
    return fill_dict


def get_init_params(param_names, def_params=default_params):
    return dict([(p, def_params[p][1]) for p in param_names])


def get_param_bounds(param_names, def_params=default_params):
    return dict([(p, def_params[p][2]) for p in param_names])


def read_params(paramfl):
    header = paramfl.next()
    hspl = header.strip('\n').split(',')
    # Check to see if we need to split into individuals
    ret = dict()
    if hspl[0] == 'WID':
        pnames = hspl[1:]
        for ln in paramfl:
            l = ln.strip('\n').split(',')
            ret[l[0]] = dict()
            for pnm, p in zip(pnames, l[1:]):
                ret[l[0]][pnm] = float(p)
    else:
        pnames = hspl
        params = paramfl.next().strip('\n').split(',')
        for pnm, p in zip(pnames, params):
            if pnm != 'LLH':
                ret[pnm] = float(p)
    paramfl.close()
    return ret


__all__ = ['default_params', 'shapes_param_names', 'materials_param_names', 'balance_param_names',
           'fill_params', 'get_init_params', 'get_param_bounds', 'cfg_dict', 'combined_param_names',
           'read_params', 'params_mod_2_rules', 'safe_remove', 'clean_param_names']
