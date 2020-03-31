from __future__ import division
import numpy as np
import copy
import os
import random
from .perception_shapes import perceive_shapes
from .perception_materials import perceive_materials, make_cp_space_materials
from .perception_balance import perceive_balance, make_cp_space_balance
from .perception_combined import perceive_combined
from .common_parts import ALL_MOD_TYPES, MSP_MOD_TYPES
from OptimTools import async_map

###########
#
# Parts of the model to combine
#
###########

# Run through a beam and check whether you can cancel out anything
def apply_symmetry(percept, mass_jnd = 0.001, dist_jnd = 0.001):
    com = percept['CoM'][0]
    conf = copy.deepcopy(percept['Configuration'])

    new_conf = []
    while len(conf) != 0:
        this_m, this_d = conf[0]
        rest = conf[1:]
        match_idx = -1
        for i, (n_m, n_d) in enumerate(rest):
            if abs((this_d - com) - (com - n_d)) < dist_jnd and \
                abs(this_m - n_m) < mass_jnd:
                match_idx = i

        conf = rest
        if match_idx == -1:
            new_conf.append((this_m, this_d))
        else:
            conf = rest
            del(conf[match_idx])

    return {'CoM': percept['CoM'], 'Configuration': new_conf, 'BeamWeight': percept['BeamWeight']}

def check_symmetry(percept, dist_jnd = 0.001):
    # Only apply balance if it's symmetric & balance is in the center
    if len(percept['Configuration']) == 0 and (abs(percept['CoM'][0]) < dist_jnd or percept['BeamWeight'] < .001):
        return 'B'
    else:
        return 'NONE'

# Checks whether total mass on each side of the balance is greater than mass_jnd - returns 'L' or 'R' if so
#  use_full_balance is a switch that checks whether to only look for blocks centered off the balance
#   (and not just on one side of the center)
def apply_mass_heur(percept, mass_jnd = .01, use_full_balance = True):
    com = percept['CoM'][0]
    if use_full_balance:
        hwid = percept['CoM'][1]
        com_l = com + hwid
        com_r = com - hwid
    else:
        com_l = com_r = com

    net_mass = 0
    for m, d in percept['Configuration']:
        if d > com_l:
            net_mass += m
        elif d < com_r:
            net_mass -= m

    if net_mass > mass_jnd:
        ret = "L"
    elif net_mass < -mass_jnd:
        ret = "R"
    else:
        ret = "NONE"

    return ret, net_mass


def apply_dist_heur(percept, dist_jnd = .01, use_dist_max = True, use_full_balance = True):
    com = percept['CoM'][0]
    if use_full_balance:
        hwid = percept['CoM'][1]
        com_l = com + hwid
        com_r = com - hwid
    else:
        com_l = com_r = com

    ld = 0
    rd = 0
    if use_dist_max:
        for _, d in percept['Configuration']:
            if d > com_l:
                ld = max(ld, d - com_l)
            elif d < com_r:
                rd = max(rd, -d + com_r)
    else:
        nl = 0
        nr = 0
        for _, d in percept['Configuration']:
            if d > com_l:
                ld += d
                nl += 1
            elif d < com_r:
                rd -= d
                nr += 1
        ld /= nl
        rd /= nr

    if (ld - rd) > dist_jnd:
        return "L"
    elif (rd - ld) > dist_jnd:
        return "R"
    else:
        return "NONE"


def calc_com(percept):
    strq = 0
    smass = percept['BeamWeight']
    for m, d in percept['Configuration']:
        strq += m*d
        smass += m
    return (strq / smass)

def apply_physics(percept):

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

# Runs physics a few times and picks using MSPRT
def sample_physics(percept, msprt_samples = 2, msprt_thresh = 0.9):
    samps = {'L':1, 'B': 1, 'R':1} # Flat prior of a single sample each
    for i in xrange(msprt_samples):
        samps[apply_physics(percept)[0]] += 1
        for ch in ['L','R','B']:
            belief = samps[ch] / (i+4)
            if belief > msprt_thresh:
                return ch
    return 'NONE'

########
#
# Run the full model on a single trial with many perceptions
# (must run through a 'perception' first
# Returns model_res, sym_pct, mass_heur_pct
#  model_res: a dict (with keys 'L', 'R', 'B') representing the model probabilities of each decision
#  sym_pct, mass_heur_pct: the proportion of decisions attributable to the symmetry & mass heuristics
#
########

def model_single_trial(percept, struct_type,
                       mass_jnd, dist_jnd, msprt_samples, msprt_thresh):
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
                s = apply_physics(percept)[0]
            else:
                s = sample_physics(percept, msprt_samples, msprt_thresh)
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

def model_trial_by_strat(percepts, struct_type,
                       mass_jnd, dist_jnd, msprt_samples=1,
                       msprt_thresh=.051, use_full_balance=True,
                       use_dist_max=True):
    N = len(percepts) # Figure out how many repeats of items
    os = [model_single_trial(p, struct_type,
                             mass_jnd, dist_jnd, msprt_samples,
                             msprt_thresh) for p in percepts]
    preds = dict([o, 0.] for o in ['L', 'R', 'B'])
    used = dict([(s, False) for s in ['m','s','p','d']])
    for pr, u in os:
        preds[pr] += 1./N
        for t in ['m','s','p','d']:
            used[t] += u[t] / N
    return preds, used

def model_trial(percepts, strats,
                mass_jnd, dist_jnd, msprt_samples=1,
                msprt_thresh=.051, use_full_balance=True,
                use_dist_max=True):
    return dict([(st, model_trial_by_strat(percepts, st, mass_jnd, dist_jnd,
                                           msprt_samples, msprt_thresh,
                                           use_full_balance, use_dist_max))
                 for st in strats])

########
#
# Run the full model on a single trial, then multiple trials
# Differs based on type of trial to perceive
#
########
def _precalc_com(percept, com_u):
    if com_u:
        cuadj = np.random.normal(0, com_u)
    else:
        cuadj = 0.
    percept['com_adjust'] = cuadj
    return percept

def model_shapes(trial, strats=ALL_MOD_TYPES, tower_mass_u = None, block_mass_u = None, distance_u = None, com_range = 1.,
                 mass_jnd=0.001, dist_jnd=0.001, com_u=None,
                 msprt_samples = 1, msprt_thresh = .51,
                 n_times = 500, use_full_balance = True, use_dist_max = True):
    percepts = [perceive_shapes(trial, tower_mass_u, block_mass_u, distance_u, com_range) for _ in range(n_times)]
    percepts = [_precalc_com(p, com_u) for p in percepts]
    return model_trial(percepts, strats, mass_jnd, dist_jnd,
                       msprt_samples, msprt_thresh,
                       use_full_balance, use_dist_max)


def model_all_shapes(modeltrials, strats=ALL_MOD_TYPES, tower_mass_u=None, block_mass_u=None, distance_u=None, com_range=1.,
                     mass_jnd=0.001, dist_jnd=0.001, com_u=None,
                     msprt_samples = 1, msprt_thresh = .51, n_times = 500, use_full_balance = True,
                     beam_type = None, use_dist_max = True, parallelize = True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return model_shapes(t, strats, tower_mass_u, block_mass_u, distance_u, com_range,
                                    mass_jnd, dist_jnd, com_u, msprt_samples, msprt_thresh, n_times,
                                    use_full_balance, use_dist_max)
    return mapfn(mfn, modeltrials)

# Materials

def model_materials(trial, strats=ALL_MOD_TYPES, block_mass_u=None,
                    brick_density=2.4, brick_u=None, iron_density=9., iron_u=None,
                    distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001,
                    com_u=None,  msprt_samples=1, msprt_thresh=.51, n_times = 500, use_full_balance = True,
                    type = 'smp', use_dist_max = True):
    percepts = [perceive_materials(trial, block_mass_u, brick_density, brick_u, iron_density, iron_u, distance_u, com_range) \
                for _ in range(n_times)]
    percepts = [_precalc_com(p, com_u) for p in percepts]
    return model_trial(percepts, strats, mass_jnd, dist_jnd,
                       msprt_samples, msprt_thresh,
                       use_full_balance, use_dist_max)

def model_all_materials(modeltrials, strats=ALL_MOD_TYPES, block_mass_u=None,
                        brick_density=2.4, brick_u=None, iron_density=9., iron_u=None,
                        distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001,
                        com_u=None, msprt_samples = 1, msprt_thresh = .51,
                        n_times = 500, use_full_balance = True, type = 'smp', beam_type=None, use_dist_max=True,
                        parallelize = True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return model_materials(t, strats, block_mass_u, brick_density, brick_u, iron_density, iron_u, distance_u,
                                       com_range, mass_jnd, dist_jnd, com_u,
                                       msprt_samples, msprt_thresh, n_times, use_full_balance, use_dist_max)
    return mapfn(mfn, modeltrials)

# Balance

def model_balance(trial, strats=ALL_MOD_TYPES, block_mass_u=None, distance_u=None,
                  com_range_vs=1., com_range_s=1., com_range_m=1.,
                  com_range_l=1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob = 1.,
                  mass_jnd=0.001, dist_jnd=0.001, com_u=None, msprt_samples = 1, msprt_thresh = .51,
                  n_times = 500, use_full_balance = True, beam_type = 'trunc_norm', use_dist_max=True):
    percepts = [perceive_balance(trial, block_mass_u, distance_u, com_range_vs, com_range_s,
                                 com_range_m, com_range_l, beam_mass_mean, beam_mass_sd,
                                 beam_use_prob, beam_type) for _ in range(n_times)]
    percepts = [_precalc_com(p, com_u) for p in percepts]
    return model_trial(percepts, strats, mass_jnd, dist_jnd,
                       msprt_samples, msprt_thresh,
                       use_full_balance, use_dist_max)

def model_all_balance(modeltrials, strats=ALL_MOD_TYPES, block_mass_u=None,
                      distance_u=None, com_range_vs=1., com_range_s=1., com_range_m=1.,
                      com_range_l=1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob = 1.,
                      mass_jnd=0.001, dist_jnd=0.001, com_u=None, msprt_samples = 1, msprt_thresh = .51,
                      n_times = 500, use_full_balance = True,
                      beam_type = 'trunc_norm', use_dist_max=True, parallelize=True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return model_balance(t, strats, block_mass_u, distance_u, com_range_vs, com_range_s, com_range_m,
                                     com_range_l, beam_mass_mean, beam_mass_sd,
                                     beam_use_prob, mass_jnd, dist_jnd, com_u, msprt_samples, msprt_thresh,
                                     n_times, use_full_balance, beam_type, use_dist_max)
    return mapfn(mfn, modeltrials)


# Combined

def model_combined(trial, strats=ALL_MOD_TYPES, block_mass_u=None, tower_mass_u=None,
                   brick_density=2.4, brick_u=None, iron_density=9.,
                   iron_u=None, distance_u=None, com_range_vs=1., com_range_m=1., beam_mass_mean=4., beam_mass_sd=None,
                   beam_use_prob=1., mass_jnd=0.001, dist_jnd=0.001, com_u=None,
                   msprt_samples=1, msprt_thresh=.51,
                   n_times = 500, use_full_balance = True, type = 'smp', beam_type = 'trunc_norm', use_dist_max=True):

    percepts = [perceive_combined(trial,block_mass_u,tower_mass_u,brick_density,brick_u,iron_density,iron_u,
                                  distance_u,com_range_vs,com_range_m,beam_mass_mean,beam_mass_sd,
                                  beam_use_prob, beam_type) for _ in range(n_times)]
    percepts = [_precalc_com(p, com_u) for p in percepts]
    return model_trial(percepts, strats, mass_jnd, dist_jnd,
                       msprt_samples, msprt_thresh,
                       use_full_balance, use_dist_max)

def model_all_combined(modeltrials, strats=ALL_MOD_TYPES, block_mass_u=None,
                       tower_mass_u=None, brick_density=2.4, brick_u=None,
                       iron_density=9., iron_u=None, distance_u=None, com_range_vs=1., com_range_m=1.,
                       beam_mass_mean=4., beam_mass_sd=None, beam_use_prob=1.,
                       mass_jnd=0.001, dist_jnd=0.001, com_u=None, msprt_samples = 1, msprt_thresh = .51,
                       n_times = 500, use_full_balance = True, type = 'smp',
                       beam_type = 'trunc_norm', use_dist_max=True, parallelize=True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return model_combined(t,strats,block_mass_u,tower_mass_u,brick_density,brick_u,
                                      iron_density,iron_u,distance_u,
                                      com_range_vs,com_range_m,beam_mass_mean,beam_mass_sd,beam_use_prob,mass_jnd, dist_jnd,com_u,msprt_samples, msprt_thresh,
                                      n_times,use_full_balance, beam_type, use_dist_max)
    return mapfn(mfn, modeltrials)


# Clean - for ferretti data (and anything else that only uses the basic rules (built off of materials data)

def model_clean(trial, strats=ALL_MOD_TYPES, block_mass_u=None, distance_u=None, com_range=1., mass_jnd=0.001,
                dist_jnd=0.001, com_u=None, msprt_samples=1, msprt_thresh=.51, n_times = 500,
                use_full_balance = True, type = 'smp', use_dist_max = True):
    percepts = [perceive_materials(trial, block_mass_u, None, None, None, None, distance_u, com_range) \
                for _ in range(n_times)]
    percepts = [_precalc_com(p, com_u) for p in percepts]
    return model_trial(percepts, strats, mass_jnd, dist_jnd,
                       msprt_samples, msprt_thresh,
                       use_full_balance, use_dist_max)

def model_all_clean(modeltrials, strats=ALL_MOD_TYPES, block_mass_u=None,
                        distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001,
                        com_u=None, msprt_samples = 1, msprt_thresh = .51,
                        n_times = 500, use_full_balance = True, type = 'smp', beam_type=None, use_dist_max=True,
                        parallelize = True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return model_clean(t, strats, block_mass_u, distance_u,
                                       com_range, mass_jnd, dist_jnd, com_u,
                                       msprt_samples, msprt_thresh, n_times, use_full_balance, type, use_dist_max)
    return mapfn(mfn, modeltrials)
