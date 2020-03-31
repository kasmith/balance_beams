from __future__ import division
from perception_shapes import perceive_shapes
from perception_materials import perceive_materials
from perception_balance import perceive_balance
from perception_combined import perceive_combined
from OptimTools import async_map, SPSA

def parse_rule_weights(percept):
    com = percept['CoM'][0]
    hwid = percept['CoM'][1]
    com_l = com+hwid
    com_r = com-hwid

    parsed_perc = []
    for m,d in percept['Configuration']:
        if d > com_l or d < com_r:
            parsed_perc.append((m, d - com))
    return parsed_perc

# Beam rules
def beam_rules(parsed, mass_jnd = 0.001, dist_jnd = 0.001):
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
            torque_r -= d*m
        else:
            mass_l += m
            maxd_l = max(d, maxd_l)
            totd_l += d
            torque_l += d*m
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
            # In a conflict, rule 3 muddles through, rule 3add checks for addition, rule4 integration
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


def rules_trial(percept, use_rule1 = .2, use_rule2 = .25, use_rule3 = .3333, use_rule3a = .5,
                mass_jnd = 0.001, dist_jnd = 0.001, lapse = 0.00001, return_all = False):

    if not (0 <= use_rule1 <= 1 and 0 <= use_rule2 <= 1 and 0 <= use_rule3 <= 1 and 0 <= use_rule3a <= 1):
        raise Exception('Illegal rule probabilities')

    ur1 = use_rule1
    ur2 = use_rule2 * (1-ur1)
    ur3 = use_rule3 * (1-ur1-ur2)
    ur3a = use_rule3a * (1-ur1-ur2-ur3)
    ur4 = 1-ur1-ur2-ur3-ur3a

    parsed = parse_rule_weights(percept)
    r1, r2, r3, r3a, r4 = beam_rules(parsed, mass_jnd, dist_jnd)

    if return_all:
        return r1, r2, r3, r3a, r4
    else:
        chdict = {'L':0, 'R':0, 'B':0}
        for u, ch in [(ur1, r1), (ur2, r2), (ur3, r3), (ur3a, r3a), (ur4, r4)]:
            if ch == 'NONE':
                chdict['L'] += u/3.
                chdict['R'] += u/3.
                chdict['B'] += u/3.
            else:
                chdict[ch] += u

        if abs(sum(chdict.values())-1) > .000000001:
            raise Exception('Just checking because I should not be here!')

        mpct = {'L': (chdict['L'] + (lapse / 3.)) / (1 + lapse),
                'R': (chdict['R'] + (lapse / 3.)) / (1 + lapse),
                'B': (chdict['B'] + (lapse / 3.)) / (1 + lapse)}
        return mpct

##### Rules models with individual parts

def rules_shapes(trial, block_mass_u=None, tower_mass_u=None, distance_u=None, com_range = 1.,
                 mass_jnd=0.001, dist_jnd=0.001, lapse=0.000001,use_rule1 = .2, use_rule2= .25,
                 use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_shapes(trial, tower_mass_u, block_mass_u, distance_u, com_range) for _ in range(n_times)]
    rules = {'L': 0, 'R': 0, 'B': 0}
    for p in percepts:
        rmod = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse)
        rules['L'] += rmod['L'] / n_times
        rules['R'] += rmod['R'] / n_times
        rules['B'] += rmod['B'] / n_times
    return rules, 0, 0, 0, {'L':0,'R':0,'B':0,'NONE':0}

def rules_all_shapes(modeltrials, block_mass_u=None, tower_mass_u=None, distance_u=None, com_range = 1.,
                     mass_jnd=0.001, dist_jnd=0.001, lapse=0.000001,use_rule1 = .2, use_rule2= .25,
                     use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm', use_full_balance=True,
                     use_dist_max = True, msprt_samples=1, type='msp', parallelize = True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_shapes(t, block_mass_u, tower_mass_u, distance_u, com_range, mass_jnd, dist_jnd, lapse,
                                      use_rule1, use_rule2, use_rule3, use_rule3a, n_times)
    return mapfn(mfn, modeltrials)


def rules_materials(trial, block_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9., iron_u=None,
                    distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001, lapse=0.000001,use_rule1 = .2, use_rule2= .25,
                    use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_materials(trial, block_mass_u, brick_density, brick_u, iron_density, iron_u,
                                   distance_u, com_range) for _ in range(n_times)]
    rules = {'L': 0, 'R': 0, 'B': 0}
    for p in percepts:
        rmod = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse)
        rules['L'] += rmod['L'] / n_times
        rules['R'] += rmod['R'] / n_times
        rules['B'] += rmod['B'] / n_times
    return rules, 0, 0, 0, {'L':0,'R':0,'B':0,'NONE':0}

def rules_all_materials(modeltrials, block_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9., iron_u=None,
                    distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001, lapse=0.000001,use_rule1 = .2, use_rule2= .25,
                    use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type = 'trunc_norm', use_full_balance=True,
                        use_dist_max=True, msprt_samples=1, type='msp', parallelize=True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_materials(t, block_mass_u, brick_density, brick_u, iron_density, iron_u, distance_u,
                                       com_range, mass_jnd, dist_jnd, lapse,
                                       use_rule1, use_rule2, use_rule3, use_rule3a, n_times)
    return mapfn(mfn, modeltrials)


def rules_balance(trial, block_mass_u=None, distance_u=None, com_range_vs=1., com_range_s=1., com_range_m=1.,
                  com_range_l = 1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob=1., mass_jnd = 0.001,
                  dist_jnd=0.001, lapse=0.000001, use_rule1 = .2, use_rule2= .25, use_rule3 = .3333, use_rule3a = .5,
                  n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_balance(trial, block_mass_u, distance_u, com_range_vs, com_range_s, com_range_m, com_range_l,
                                 beam_mass_mean, beam_mass_sd, beam_use_prob, beam_type) for _ in range(n_times)]
    rules = {'L': 0, 'R': 0, 'B': 0}
    for p in percepts:
        rmod = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse)
        rules['L'] += rmod['L'] / n_times
        rules['R'] += rmod['R'] / n_times
        rules['B'] += rmod['B'] / n_times
    return rules, 0, 0, 0, {'L':0,'R':0,'B':0,'NONE':0}

def rules_all_balance(modeltrials, block_mass_u=None, distance_u=None, com_range_vs=1., com_range_s=1., com_range_m=1.,
                  com_range_l = 1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob=1., mass_jnd = 0.001,
                  dist_jnd=0.001, lapse=0.000001, use_rule1 = .2, use_rule2= .25, use_rule3 = .3333, use_rule3a = .5,
                  n_times=500, beam_type='trunc_norm', use_full_balance=True, use_dist_max = True, msprt_samples=1,
                      type='msp', parallelize = True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_balance(t, block_mass_u, distance_u, com_range_vs, com_range_s, com_range_m, com_range_l,
                                     beam_mass_mean, beam_mass_sd, beam_use_prob, mass_jnd, dist_jnd, lapse, use_rule1,
                                     use_rule2, use_rule3, use_rule3a, n_times, beam_type)
    return mapfn(mfn, modeltrials)


def rules_combined(trial, block_mass_u=None, tower_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9.,
                   iron_u=None, distance_u=None, com_range_vs=1., com_range_m=1., beam_mass_mean=4., beam_mass_sd=None,
                   beam_use_prob=1., mass_jnd = 0.001, dist_jnd=0.001, lapse=0.000001, use_rule1 = .2, use_rule2= .25,
                   use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_combined(trial, block_mass_u, tower_mass_u, brick_density, brick_u, iron_density, iron_u,
                                  distance_u, com_range_vs, com_range_m, beam_mass_mean, beam_mass_sd,
                                  beam_use_prob, beam_type) for _ in range(n_times)]
    rules = {'L':0, 'R': 0, 'B': 0}
    for p in percepts:
        rmod = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse)
        rules['L'] += rmod['L'] / n_times
        rules['R'] += rmod['R'] / n_times
        rules['B'] += rmod['B'] / n_times
    return rules, 0, 0, 0, {'L':0,'R':0,'B':0,'NONE':0}

def rules_all_combined(modeltrials, block_mass_u=None, tower_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9.,
                   iron_u=None, distance_u=None, com_range_vs=1., com_range_m=1., beam_mass_mean=4., beam_mass_sd=None,
                   beam_use_prob=1., mass_jnd = 0.001, dist_jnd=0.001, lapse=0.000001, use_rule1 = .2, use_rule2= .25,
                   use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm', parallelize = True,
                   use_dist_max = True, use_full_balance=True, msprt_samples=1, type='msp'):

    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_combined(t, block_mass_u, tower_mass_u, brick_density, brick_u, iron_density,
                                      iron_u, distance_u, com_range_vs, com_range_m, beam_mass_mean,
                                      beam_mass_sd, beam_use_prob, mass_jnd, dist_jnd, lapse, use_rule1,
                                      use_rule2, use_rule3, use_rule3a, n_times, beam_type)
    return mapfn(mfn, modeltrials)


# Get a set of the rules out (reduces computation for fitting individual)

def rules_shapes_full(trial, block_mass_u=None, tower_mass_u=None, distance_u=None, com_range = 1.,
                 mass_jnd=0.001, dist_jnd=0.001, lapse=0.000001,use_rule1 = .2, use_rule2= .25,
                 use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_shapes(trial, tower_mass_u, block_mass_u, distance_u, com_range) for _ in range(n_times)]
    rules = dict([(r, {'L': 0, 'R': 0, 'B': 0}) for r in ['R1','R2','R3','R3a','R4']])
    for p in percepts:
        rs = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse,
                                      return_all=True)

        for rnm, rch in zip(['R1','R2','R3','R3a','R4'], rs):
            if rch != 'NONE':
                rules[rnm][rch] += 1. / n_times
            else:
                rules[rnm]['L'] += 1. / (3. * n_times)
                rules[rnm]['B'] += 1. / (3. * n_times)
                rules[rnm]['R'] += 1. / (3. * n_times)


    for rnm in ['R1','R2','R3','R3a','R4']:
        rules[rnm]['L'] = (rules[rnm]['L'] + (lapse / 3.)) / (1 + lapse)
        rules[rnm]['B'] = (rules[rnm]['B'] + (lapse / 3.)) / (1 + lapse)
        rules[rnm]['R'] = (rules[rnm]['R'] + (lapse / 3.)) / (1 + lapse)

    return rules, 0, 0, 0, {'L':0,'R':0,'B':0,'NONE':0}

def rules_all_shapes_full(modeltrials, block_mass_u=None, tower_mass_u=None, distance_u=None, com_range = 1.,
                     mass_jnd=0.001, dist_jnd=0.001, lapse=0.000001,use_rule1 = .2, use_rule2= .25,
                     use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm', use_full_balance=True,
                     use_dist_max = True, msprt_samples=1, type='msp', parallelize = True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_shapes_full(t, block_mass_u, tower_mass_u, distance_u, com_range, mass_jnd, dist_jnd, lapse,
                                      use_rule1, use_rule2, use_rule3, use_rule3a, n_times)
    return mapfn(mfn, modeltrials)

def rules_materials_full(trial, block_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9., iron_u=None,
                    distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001, lapse=0.000001,use_rule1 = .2, use_rule2= .25,
                    use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_materials(trial, block_mass_u, brick_density, brick_u, iron_density, iron_u,
                                   distance_u, com_range) for _ in range(n_times)]
    rules = dict([(r, {'L': 0, 'R': 0, 'B': 0}) for r in ['R1', 'R2', 'R3', 'R3a', 'R4']])
    for p in percepts:
        rs = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse,
                         return_all=True)

        for rnm, rch in zip(['R1', 'R2', 'R3', 'R3a', 'R4'], rs):
            if rch != 'NONE':
                rules[rnm][rch] += 1. / n_times
            else:
                rules[rnm]['L'] += 1. / (3. * n_times)
                rules[rnm]['B'] += 1. / (3. * n_times)
                rules[rnm]['R'] += 1. / (3. * n_times)

    for rnm in ['R1', 'R2', 'R3', 'R3a', 'R4']:
        rules[rnm]['L'] = (rules[rnm]['L'] + (lapse / 3.)) / (1 + lapse)
        rules[rnm]['B'] = (rules[rnm]['B'] + (lapse / 3.)) / (1 + lapse)
        rules[rnm]['R'] = (rules[rnm]['R'] + (lapse / 3.)) / (1 + lapse)

    return rules, 0, 0, 0, {'L': 0, 'R': 0, 'B': 0, 'NONE': 0}

def rules_all_materials_full(modeltrials, block_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9., iron_u=None,
                    distance_u=None, com_range=1., mass_jnd=0.001, dist_jnd=0.001, lapse=0.000001,use_rule1 = .2, use_rule2= .25,
                    use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type = 'trunc_norm', use_full_balance=True,
                        use_dist_max=True, msprt_samples=1, type='msp', parallelize=True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_materials_full(t, block_mass_u, brick_density, brick_u, iron_density, iron_u, distance_u,
                                       com_range, mass_jnd, dist_jnd, lapse,
                                       use_rule1, use_rule2, use_rule3, use_rule3a, n_times)
    return mapfn(mfn, modeltrials)


def rules_balance_full(trial, block_mass_u=None, distance_u=None, com_range_vs=1., com_range_s=1., com_range_m=1.,
                  com_range_l = 1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob=1., mass_jnd = 0.001,
                  dist_jnd=0.001, lapse=0.000001, use_rule1 = .2, use_rule2= .25, use_rule3 = .3333, use_rule3a = .5,
                  n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_balance(trial, block_mass_u, distance_u, com_range_vs, com_range_s, com_range_m, com_range_l,
                                 beam_mass_mean, beam_mass_sd, beam_use_prob, beam_type) for _ in range(n_times)]
    rules = dict([(r, {'L': 0, 'R': 0, 'B': 0}) for r in ['R1', 'R2', 'R3', 'R3a', 'R4']])
    for p in percepts:
        rs = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse,
                         return_all=True)

        for rnm, rch in zip(['R1', 'R2', 'R3', 'R3a', 'R4'], rs):
            if rch != 'NONE':
                rules[rnm][rch] += 1. / n_times
            else:
                rules[rnm]['L'] += 1. / (3. * n_times)
                rules[rnm]['B'] += 1. / (3. * n_times)
                rules[rnm]['R'] += 1. / (3. * n_times)

    for rnm in ['R1', 'R2', 'R3', 'R3a', 'R4']:
        rules[rnm]['L'] = (rules[rnm]['L'] + (lapse / 3.)) / (1 + lapse)
        rules[rnm]['B'] = (rules[rnm]['B'] + (lapse / 3.)) / (1 + lapse)
        rules[rnm]['R'] = (rules[rnm]['R'] + (lapse / 3.)) / (1 + lapse)

    return rules, 0, 0, 0, {'L': 0, 'R': 0, 'B': 0, 'NONE': 0}

def rules_all_balance_full(modeltrials, block_mass_u=None, distance_u=None, com_range_vs=1., com_range_s=1., com_range_m=1.,
                  com_range_l = 1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob=1., mass_jnd = 0.001,
                  dist_jnd=0.001, lapse=0.000001, use_rule1 = .2, use_rule2= .25, use_rule3 = .3333, use_rule3a = .5,
                  n_times=500, beam_type='trunc_norm', use_full_balance=True, use_dist_max = True, msprt_samples=1,
                      type='msp', parallelize = True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_balance_full(t, block_mass_u, distance_u, com_range_vs, com_range_s, com_range_m, com_range_l,
                                     beam_mass_mean, beam_mass_sd, beam_use_prob, mass_jnd, dist_jnd, lapse, use_rule1,
                                     use_rule2, use_rule3, use_rule3a, n_times, beam_type)
    return mapfn(mfn, modeltrials)

def rules_combined_full(trial, block_mass_u=None, tower_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9.,
                   iron_u=None, distance_u=None, com_range_vs=1., com_range_m=1., beam_mass_mean=4., beam_mass_sd=None,
                   beam_use_prob=1., mass_jnd = 0.001, dist_jnd=0.001, lapse=0.000001, use_rule1 = .2, use_rule2= .25,
                   use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_combined(trial, block_mass_u, tower_mass_u, brick_density, brick_u, iron_density, iron_u,
                                  distance_u, com_range_vs, com_range_m, beam_mass_mean, beam_mass_sd,
                                  beam_use_prob, beam_type) for _ in range(n_times)]
    rules = dict([(r, {'L': 0, 'R': 0, 'B': 0}) for r in ['R1', 'R2', 'R3', 'R3a', 'R4']])
    for p in percepts:
        rs = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse,
                         return_all=True)

        for rnm, rch in zip(['R1', 'R2', 'R3', 'R3a', 'R4'], rs):
            if rch != 'NONE':
                rules[rnm][rch] += 1. / n_times
            else:
                rules[rnm]['L'] += 1. / (3. * n_times)
                rules[rnm]['B'] += 1. / (3. * n_times)
                rules[rnm]['R'] += 1. / (3. * n_times)

    for rnm in ['R1', 'R2', 'R3', 'R3a', 'R4']:
        rules[rnm]['L'] = (rules[rnm]['L'] + (lapse / 3.)) / (1 + lapse)
        rules[rnm]['B'] = (rules[rnm]['B'] + (lapse / 3.)) / (1 + lapse)
        rules[rnm]['R'] = (rules[rnm]['R'] + (lapse / 3.)) / (1 + lapse)

    return rules, 0, 0, 0, {'L': 0, 'R': 0, 'B': 0, 'NONE': 0}

def rules_all_combined_full(modeltrials, block_mass_u=None, tower_mass_u=None, brick_density=2.4, brick_u=None, iron_density=9.,
                   iron_u=None, distance_u=None, com_range_vs=1., com_range_m=1., beam_mass_mean=4., beam_mass_sd=None,
                   beam_use_prob=1., mass_jnd = 0.001, dist_jnd=0.001, lapse=0.000001, use_rule1 = .2, use_rule2= .25,
                   use_rule3 = .3333, use_rule3a = .5, n_times=500, beam_type='trunc_norm', parallelize = True,
                   use_dist_max = True, use_full_balance=True, msprt_samples=1, type='msp'):

    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t): return rules_combined_full(t, block_mass_u, tower_mass_u, brick_density, brick_u, iron_density,
                                      iron_u, distance_u, com_range_vs, com_range_m, beam_mass_mean,
                                      beam_mass_sd, beam_use_prob, mass_jnd, dist_jnd, lapse, use_rule1,
                                      use_rule2, use_rule3, use_rule3a, n_times, beam_type)
    return mapfn(mfn, modeltrials)


def rules_clean(trial, block_mass_u=None, distance_u=None, com_range=1., mass_heur_p=0., dist_heur_p=0., mass_jnd=0.001,
                dist_jnd=0.001, com_u=None, lapse=0.000001, use_rule1=.2, use_rule2=.25, use_rule3=.3333, use_rule3a=.5,
                n_times=500, beam_type='trunc_norm'):
    percepts = [perceive_materials(trial, block_mass_u, None, None, None, None, distance_u, com_range) \
                for _ in range(n_times)]
    rules = {'L': 0, 'R': 0, 'B': 0}
    for p in percepts:
        rmod = rules_trial(p, use_rule1, use_rule2, use_rule3, use_rule3a, mass_jnd, dist_jnd, lapse)
        rules['L'] += rmod['L'] / n_times
        rules['R'] += rmod['R'] / n_times
        rules['B'] += rmod['B'] / n_times
    return rules, 0, 0, 0, {'L': 0, 'R': 0, 'B': 0, 'NONE': 0}



def rules_all_clean(modeltrials, block_mass_u=None, distance_u=None, com_range=1., mass_heur_p=0., dist_heur_p=0., mass_jnd=0.001,
                dist_jnd=0.001, com_u=None, lapse=0.000001, use_rule1=.2, use_rule2=.25, use_rule3=.3333, use_rule3a=.5,
                n_times=500, beam_type='trunc_norm', use_full_balance=True, use_dist_max = True, msprt_samples=1,
                      type='msp', parallelize = True):
    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    def mfn(t):
        return rules_clean(t, block_mass_u, distance_u, com_range, mass_heur_p, dist_heur_p, mass_jnd,
                           dist_jnd, com_u, lapse, use_rule1, use_rule2, use_rule3, use_rule3a, n_times, beam_type)

    return mapfn(mfn, modeltrials)
