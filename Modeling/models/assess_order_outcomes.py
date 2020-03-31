from __future__ import division
import shapes
import materials
import balance
import combined
from model_parts import *
from config import read_params, default_params
import os
import json
import random
import copy
import numpy as np

USED_PARAMS_FL = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                              '..','output','complete', 'all_params.csv')

DEF_PERCEPTS_FL = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               '..', 'output', 'rationality', 'percepts.json')
DEF_MODEL_FL = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               '..', 'output', 'rationality',
                               'mod_predictions.csv')

def make_mod_types(rem=['s', 'm', 'd', 'p'], base_type=''):
    if len(rem) == 0:
        return []
    ret = []
    for t in rem:
        new_base = base_type + t
        new_rem = [r for r in rem if r != t]
        ret += [new_base]
        ret += make_mod_types(new_rem, new_base)
    return ret

MOD_TYPES = make_mod_types()

params = read_params(open(USED_PARAMS_FL,'rU'))
for k in default_params.keys():
    if k not in params.keys():
        params[k] = default_params[k][0]

def _perceive_shapes_precalc(trial, tmu, bmu, du, cr, cu):
    perc = perceive_shapes(trial, tmu, bmu, du, cr)
    if cu:
        cuadj = np.random.normal(0, cu)
    else:
        cuadj = 0.
    perc['com_adjust'] = cuadj
    return perc

def _perceive_materials_precalc(trial, bmu, bd, bu, id_, iu, du, cr, cu):
    perc = perceive_materials(trial, bmu, bd, bu, id_, iu, du, cr)
    if cu:
        cuadj = np.random.normal(0, cu)
    else:
        cuadj = 0.
    perc['com_adjust'] = cuadj
    return perc

def _perceive_pivot_precalc(trial, bmu, du, crvs, crs, crm, crl, bmm, bmsd,
                            bup, cu):
    perc = perceive_balance(trial, bmu, du, crvs, crs, crm, crl, bmm, bmsd,
                            bup, 'exponential')
    if cu:
        cuadj = np.random.normal(0, cu)
    else:
        cuadj = 0.
    perc['com_adjust'] = cuadj
    return perc

def _perceive_combined_precalc(trial, bmu, tmu, bd, bu, id_, iu, du, crvs,
                               crm, bmm, bmsd, bup, cu):
    perc = perceive_combined(trial, bmu, tmu, bd, bu, id_, iu, du, crvs, crm,
                             bmm, bmsd, bup, 'exponential')
    if cu:
        cuadj = np.random.normal(0, cu)
    else:
        cuadj = 0.
    perc['com_adjust'] = cuadj
    return perc

def make_percepts_dict(n=100, oflnm=DEF_PERCEPTS_FL):
    percepts = {'N': n}
    # Shapes experiment
    sh_dict = {}
    for trnm, tr in shapes.shapes_trials.items():
        percs = [_perceive_shapes_precalc(tr, params['tower_mass_u'],
                               params['block_mass_u'], params['distance_u'],
                               params['com_range'], params['com_u'])
                 for _ in range(n)]
        sh_dict[trnm] = percs
        print 'Shapes', trnm, 'done'
    # Materials experiment
    mat_dict = {}
    for trnm, tr in materials.materials_trials.items():
        percs = [_perceive_materials_precalc(tr, params['block_mass_u'],
                                    params['brick_density'], params['brick_u'],
                                    params['iron_density'], params['iron_u'],
                                    params['distance_u'], params['com_range'],
                                    params['com_u'])
                 for _ in range(n)]
        mat_dict[trnm] = percs
        print 'Materials', trnm, 'done'
    # Pivot experiment
    piv_dict = {}
    for trnm, tr in balance.balance_trials.items():
        percs = [_perceive_pivot_precalc(tr, params['block_mass_u'],
                                  params['distance_u'], params['com_range_vs'],
                                  params['com_range_s'], params['com_range_m'],
                                  params['com_range_l'],
                                  params['beam_mass_mean'],
                                  params['beam_mass_sd'],
                                  params['beam_use_prob'], params['com_u'])
                 for _ in range(n)]
        piv_dict[trnm] = percs
        print "Pivot", trnm, 'done'
    # Compbined experiment
    comb_dict = {}
    for trnm, tr in combined.combined_trials.items():
        percs = [_perceive_combined_precalc(tr, params['block_mass_u'],
                                   params['tower_mass_u'],
                                   params['brick_density'], params['brick_u'],
                                   params['iron_density'], params['iron_u'],
                                   params['distance_u'],
                                   params['com_range_vs'],
                                   params['com_range_m'],
                                   params['beam_mass_mean'],
                                   params['beam_mass_sd'],
                                   params['beam_use_prob'],
                                   params['com_u'])
                 for _ in range(n)]
        comb_dict[trnm] = percs
        print "Combined", trnm, 'done'
    # Stitch together
    percepts['shapes'] = sh_dict
    percepts['materials'] = mat_dict
    percepts['pivot'] = piv_dict
    percepts['combined'] = comb_dict
    with open(oflnm, 'w') as ofl:
        json.dump(percepts, ofl)
    return percepts

def _apply_physics_precalc(percept):
    com = calc_com(percept)
    bcent = percept['CoM'][0] + percept['com_adjust']
    pwid = percept['CoM'][2]
    if com > bcent + pwid:
        ret = 'L'
    elif com < bcent - pwid:
        ret = 'R'
    else:
        ret = 'B'
    return ret, com

def model_single_trial(percept, struct_type, mass_heur_p, dist_heur_p,
                       mass_jnd, dist_jnd, com_u, msprt_samples, msprt_thresh):
    # Make sure this is a legal struct_type -- all unique
    assert (len(struct_type) <= 4 and
            len(set(struct_type)) == len(struct_type) and
            all([l in ['m','s','p','d'] for l in struct_type])
            ), "All struct_type letters must be in [m,s,p,d]"
    used = dict([(s, False) for s in ['m','s','p','d']])

    # Try to solve according to struct_type
    for t in struct_type:
        if t == 's':
            used['s'] = True
            symp = apply_symmetry(percept, mass_jnd, dist_jnd)
            if check_symmetry(symp, dist_jnd) == 'B':
                return 'B', used
        elif t == 'm':
            used['m'] = True
            m = apply_mass_heur(percept, mass_jnd)[0]
            if m in ['L', 'R']:
                return m, used
        elif t == 'd':
            used['d'] = True
            d = apply_dist_heur(percept, dist_jnd)
            if d in ['L', 'R']:
                return d, used
        elif t == 'p':
            used['p'] = True
            if msprt_samples == 1:
                s = _apply_physics_precalc(percept)[0]
            else:
                s = sample_physics(percept, com_u, msprt_samples, msprt_thresh)
            if s in ['L', 'R', 'B']:
                return s, used
        else:
            raise RuntimeError("Struct type not found!")
    # If we get here, everything has failed -- just guess
    r = random.random()
    if r < (1./3.):
        return 'L', used
    elif r < (2./3.):
        return 'R', used
    else:
        return 'B', used


def model_multi_trial(N, percepts, struct_type, mass_heur_p, dist_heur_p,
                       mass_jnd, dist_jnd, com_u, msprt_samples, msprt_thresh):
    os = [model_single_trial(p, struct_type, mass_heur_p, dist_heur_p,
                             mass_jnd, dist_jnd, com_u, msprt_samples,
                             msprt_thresh) for p in percepts]
    preds = dict([o, 0.] for o in ['L', 'R', 'B'])
    used = dict([(s, False) for s in ['m','s','p','d']])
    for pr, u in os:
        preds[pr] += 1./N
        for t in ['m','s','p','d']:
            used[t] += u[t] / N
    return preds, used

def assess_all_orders(percepts, params, oflnm=DEF_MODEL_FL, mods = MOD_TYPES,
                      phys_msprt_samps=1, phys_msprt_thresh=.5):
    N = percepts['N']
    with open(oflnm, 'w') as ofl:
        ofl.write('Experiment,Trial,OrderType,ModLeft,ModBal,ModRight,')
        ofl.write('UsesSym,UsesMass,UsesDist,UsesPhys\n')
        for expnm, expdat in percepts.items():
            if expnm != "N":
                for trnm, trdat in expdat.items():
                    for mt in mods:
                        preds, used = model_multi_trial(N, trdat, mt,
                                                params['mass_heur_p'],
                                                params['dist_heur_p'],
                                                params['mass_jnd'],
                                                params['dist_jnd'],
                                                params['com_u'],
                                                phys_msprt_samps,
                                                phys_msprt_thresh)
                        ofl.write(expnm + ',' + trnm + ',' + mt + ',')
                        ofl.write(str(preds['L']) + ',' + str(preds['B']) + ',')
                        ofl.write(str(preds['R']) + ',' + str(used['s']) + ',')
                        ofl.write(str(used['m']) + ',' + str(used['d']) + ',')
                        ofl.write(str(used['p']) + '\n')
                    print 'Done with', expnm, trnm
    print 'All done'



if __name__ == '__main__':
    #percepts = make_percepts_dict()
    percepts = json.load(open(DEF_PERCEPTS_FL, 'rU'))
    #assess_all_orders(percepts, params)

    nouparams = copy.deepcopy(params)
    nouparams['distance_u'] = 0.000001
    nouparams['block_mass_u'] = 0.000001
    noumodfl = os.path.join(os.path.dirname(DEF_MODEL_FL),
                            "mod_no_localization_noise.csv")
    assess_all_orders(percepts, nouparams, noumodfl)
