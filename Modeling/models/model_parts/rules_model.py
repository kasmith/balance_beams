from __future__ import division
from .perception_shapes import perceive_shapes
from .perception_materials import perceive_materials
from .perception_balance import perceive_balance
from .perception_combined import perceive_combined
from OptimTools import async_map, SPSA

RULE_TYPES = ['R1', 'R2', 'R3', 'R3a', 'R4']

def parse_rule_weights(percept):
    com = percept['CoM'][0]
    hwid = percept['CoM'][1]
    com_l = com + hwid
    com_r = com - hwid

    parsed_perc = []
    for m, d in percept['Configuration']:
        if d > com_l or d < com_r:
            parsed_perc.append((m, d - com))
    return parsed_perc

# Beam rules


def beam_rules(parsed, mass_jnd=0.001, dist_jnd=0.001):
    mass_r = 0
    mass_l = 0
    maxd_l = 0
    maxd_r = 0
    totd_l = 0
    totd_r = 0
    torque_l = 0
    torque_r = 0

    for m, d in parsed:
        if d < 0:
            mass_r += m
            maxd_r = max(-d, maxd_r)
            totd_r -= d
            torque_r -= d * m
        else:
            mass_l += m
            maxd_l = max(d, maxd_l)
            totd_l += d
            torque_l += d * m
    netmass = mass_l - mass_r

    wg = abs(netmass) > mass_jnd
    dg = abs(maxd_l - maxd_r) > dist_jnd

    # If there is no difference
    if not wg and not dg:
        rule_1 = rule_2 = rule_3 = rule_3add = rule_4 = 'B'
    # If the weight is greater but distance is not
    elif wg and not dg:
        if netmass > mass_jnd:
            rule_1 = rule_2 = rule_3 = rule_3add = rule_4 = 'L'
        else:
            rule_1 = rule_2 = rule_3 = rule_3add = rule_4 = 'R'
    # If the distance is greater but weight is not
    elif dg and not wg:
        rule_1 = 'B'
        if (maxd_l - maxd_r) > dist_jnd:
            rule_2 = rule_3 = rule_3add = rule_4 = 'L'
        else:
            rule_2 = rule_3 = rule_3add = rule_4 = 'R'
    # If there is a difference in both weight & distance
    else:
        # Rule 1 & 2 go with weight
        if netmass > mass_jnd:
            rule_1 = rule_2 = 'L'
        else:
            rule_1 = rule_2 = 'R'
        # Other rules check for side consistency
        if netmass > mass_jnd and (maxd_l - maxd_r) > dist_jnd:
            rule_3 = rule_3add = rule_4 = 'L'
        elif netmass < -mass_jnd and (maxd_r - maxd_l) > dist_jnd:
            rule_3 = rule_3add = rule_4 = 'R'
        else:
            # In a conflict, rule 3 muddles through, rule 3add checks for
            # addition, rule4 integration
            rule_3 = 'NONE'
            if mass_l + totd_l > mass_r + totd_r:
                rule_3add = 'L'
            elif mass_l + totd_l < mass_r + totd_r:
                rule_3add = 'R'
            else:
                rule_3add = 'B'
            if torque_l > torque_r:
                rule_4 = 'L'
            elif torque_l < torque_r:
                rule_4 = 'R'
            else:
                rule_4 = 'B'

    return rule_1, rule_2, rule_3, rule_3add, rule_4


def rules_trial(percepts, mass_jnd=0.001, dist_jnd=0.001):

    nsims = len(percepts)
    rtypes = RULE_TYPES
    rules = [beam_rules(parse_rule_weights(p), mass_jnd, dist_jnd)
             for p in percepts]

    routs = dict([[r, [{'L': 0., 'R': 0., 'B': 0.}, []]] for r in rtypes])

    for r1, r2, r3, r3a, r4 in rules:
        routs['R1'][0][r1] += 1.
        routs['R2'][0][r2] += 1.
        if r3 == 'NONE':
            routs['R3'][0]['L'] += 1./3.
            routs['R3'][0]['R'] += 1./3.
            routs['R3'][0]['B'] += 1./3.
        else:
            routs['R3'][0][r3] += 1.
        routs['R3a'][0][r3a] += 1.
        routs['R4'][0][r4] += 1.

    for rt in rtypes:
        routs[rt][0]['L'] /= nsims
        routs[rt][0]['R'] /= nsims
        routs[rt][0]['B'] /= nsims

    return routs

# Rules models with individual parts


def rules_shapes(trial, strats=RULE_TYPES, block_mass_u=None, tower_mass_u=None, distance_u=None, com_range=1.,
                 mass_jnd=0.001, dist_jnd=0.001, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_shapes(
        trial, tower_mass_u, block_mass_u, distance_u, com_range) for _ in range(n_times)]
    return rules_trial(percepts, mass_jnd, dist_jnd)


def rules_all_shapes(modeltrials, strats=RULE_TYPES, block_mass_u=None, tower_mass_u=None, distance_u=None, com_range=1.,
                     mass_jnd=0.001, dist_jnd=0.001, n_times=500, beam_type='trunc_norm', use_full_balance=True,
                     use_dist_max=True, msprt_samples=1, type='msp', parallelize=True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_shapes(t, block_mass_u, tower_mass_u, distance_u, com_range, mass_jnd, dist_jnd, n_times)
    return mapfn(mfn, modeltrials)


def rules_materials(trial, strats=RULE_TYPES, block_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9., iron_u=None,
                    distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_materials(trial, block_mass_u, brick_density, brick_u, iron_density, iron_u,
                                   distance_u, com_range) for _ in range(n_times)]
    return rules_trial(percepts, mass_jnd, dist_jnd)


def rules_all_materials(modeltrials, strats=RULE_TYPES, block_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9., iron_u=None,
                        distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001, n_times=500, beam_type='trunc_norm', use_full_balance=True,
                        use_dist_max=True, msprt_samples=1, type='msp', parallelize=True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_materials(t, block_mass_u, brick_density, brick_u, iron_density, iron_u, distance_u,
                                       com_range, mass_jnd, dist_jnd, n_times)
    return mapfn(mfn, modeltrials)


def rules_balance(trial, strats=RULE_TYPES, block_mass_u=None, distance_u=None, com_range_vs=1., com_range_s=1., com_range_m=1.,
                  com_range_l=1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob=1., mass_jnd=0.001,
                  dist_jnd=0.001, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_balance(trial, block_mass_u, distance_u, com_range_vs, com_range_s, com_range_m, com_range_l,
                                 beam_mass_mean, beam_mass_sd, beam_use_prob, beam_type) for _ in range(n_times)]
    return rules_trial(percepts, mass_jnd, dist_jnd)


def rules_all_balance(modeltrials, strats=RULE_TYPES, block_mass_u=None, distance_u=None, com_range_vs=1., com_range_s=1., com_range_m=1.,
                      com_range_l=1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob=1., mass_jnd=0.001,
                      dist_jnd=0.001, n_times=500, beam_type='trunc_norm', use_full_balance=True, use_dist_max=True, msprt_samples=1,
                      type='msp', parallelize=True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_balance(t, block_mass_u, distance_u, com_range_vs, com_range_s, com_range_m, com_range_l,
                                     beam_mass_mean, beam_mass_sd, beam_use_prob, mass_jnd, dist_jnd, n_times, beam_type)
    return mapfn(mfn, modeltrials)


def rules_combined(trial, strats=RULE_TYPES, block_mass_u=None, tower_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9.,
                   iron_u=None, distance_u=None, com_range_vs=1., com_range_m=1., beam_mass_mean=4., beam_mass_sd=None,
                   beam_use_prob=1., mass_jnd=0.001, dist_jnd=0.001, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_combined(trial, block_mass_u, tower_mass_u, brick_density, brick_u, iron_density, iron_u,
                                  distance_u, com_range_vs, com_range_m, beam_mass_mean, beam_mass_sd,
                                  beam_use_prob, beam_type) for _ in range(n_times)]
    return rules_trial(percepts, mass_jnd, dist_jnd)


def rules_all_combined(modeltrials, strats=RULE_TYPES, block_mass_u=None, tower_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9.,
                       iron_u=None, distance_u=None, com_range_vs=1., com_range_m=1., beam_mass_mean=4., beam_mass_sd=None,
                       beam_use_prob=1., mass_jnd=0.001, dist_jnd=0.001, n_times=500, beam_type='trunc_norm', parallelize=True,
                       use_dist_max=True, use_full_balance=True, msprt_samples=1, type='msp'):

    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_combined(t, block_mass_u, tower_mass_u, brick_density, brick_u, iron_density,
                                      iron_u, distance_u, com_range_vs, com_range_m, beam_mass_mean,
                                      beam_mass_sd, beam_use_prob, mass_jnd, dist_jnd, n_times, beam_type)
    return mapfn(mfn, modeltrials)


def rules_clean(trial, strats=RULE_TYPES, block_mass_u=None, distance_u=None, com_range=1., mass_heur_p=0., dist_heur_p=0., mass_jnd=0.001,
                dist_jnd=0.001, com_u=None,
                n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_materials(trial, block_mass_u, None, None, None, None, distance_u, com_range)
                for _ in range(n_times)]
    return rules_trial(percepts, mass_jnd, dist_jnd)


def rules_all_clean(modeltrials, strats=RULE_TYPES, block_mass_u=None, distance_u=None, com_range=1., mass_heur_p=0., dist_heur_p=0., mass_jnd=0.001,
                    dist_jnd=0.001, com_u=None,
                    n_times=500, beam_type='trunc_norm', use_full_balance=True, use_dist_max=True, msprt_samples=1,
                    type='msp', parallelize=True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t):
        return rules_clean(t, block_mass_u, distance_u, com_range, mass_heur_p, dist_heur_p, mass_jnd,
                           dist_jnd, com_u, n_times, beam_type)

    return mapfn(mfn, modeltrials)
