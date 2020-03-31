from __future__ import division
import numpy as np
from .common_parts import truncNorm, getBeamWeight

def _perceive_combined_configuration(beam, block_mass_u=None, tower_mass_u = None,
                                     brick_density=2.4, brick_u=None, iron_density=9., iron_u=None, distance_u=None):
    perc_config = []

    if brick_u:
        u_brick = brick_density * np.random.lognormal(0, brick_u)
    else:
        u_brick = brick_density
    if iron_u:
        u_iron = iron_density * np.random.lognormal(0, iron_u)
    else:
        u_iron = iron_density

    for d, tp, sz, mat in beam:

        mass = sz
        if mat == 'B':
            mass *= u_brick
        if mat == 'I':
            mass *= u_iron

        if tp == 'BLOCK' and block_mass_u:
            mass *= np.random.lognormal(0, block_mass_u)
        elif tp != 'BLOCK' and tower_mass_u:
            mass *= np.random.lognormal(0, tower_mass_u)

        if distance_u:
            d += np.random.normal(0, distance_u)

        perc_config.append((mass, -d))

    return perc_config

# Negatives because left is positive in the model, but negative in Blender
def perceive_combined(trial, block_mass_u=None, tower_mass_u = None,
                     brick_density=2.4, brick_u=None, iron_density=9., iron_u=None, distance_u=None,
                     com_range_vs=1., com_range_m=1., beam_mass_mean=4., beam_mass_sd=None, beam_use_prob =1.,
                      beam_type = 'trunc_norm'):

    strut_pos = trial['StrutPos']
    strut_wid = trial['StrutSize']
    beam = trial['Beam']

    if strut_wid == .25:
        com_range = com_range_vs
    elif strut_wid == 1.:
        com_range = com_range_m
    else:
        raise Exception('Uh oh - bad strut reading')

    bw = getBeamWeight(beam_mass_mean, beam_mass_sd, beam_use_prob, beam_type)

    return {'CoM': (-strut_pos, strut_wid, com_range),
            'Configuration': _perceive_combined_configuration(beam, block_mass_u, tower_mass_u, brick_density,
                                                              brick_u, iron_density, iron_u, distance_u),
            'BeamWeight': bw}
