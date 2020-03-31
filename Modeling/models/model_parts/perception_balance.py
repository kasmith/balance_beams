from __future__ import division
import numpy as np
from .common_parts import truncNorm, getBeamWeight
from .physics_helper import PhysicsWorld

def _perceive_balance_configuration(trial, block_mass_u=None, distance_u=None):
    perc_config = []

    for d, bls in trial:
        totmass = 0
        if block_mass_u:
            totmass += bls*np.random.lognormal(0, block_mass_u)
        else:
            totmass += bls
        if distance_u:
            d += np.random.normal(0, distance_u

                                  )

        perc_config.append((totmass, -d))

    return perc_config

# Negatives because left is positive in the model, but negative in Blender
def perceive_balance(trial, block_mass_u=None, distance_u=None,
                     com_range_vs=1., com_range_s=1., com_range_m=1., com_range_l=1.,
                     beam_mass_mean=4., beam_mass_sd=None, beam_use_prob =1., beam_type = 'trunc_norm'):

    strut, beam = trial
    if strut[1] == .25:
        com_range = com_range_vs
    elif strut[1] == .5:
        com_range = com_range_s
    elif strut[1] == 1.:
        com_range = com_range_m
    elif strut[1] == 2.:
        com_range = com_range_l
    else:
        raise Exception('Uh oh - bad strut reading')

    bw = getBeamWeight(beam_mass_mean, beam_mass_sd, beam_use_prob, beam_type)

    return {'CoM': (-strut[0], strut[1], com_range),
            'Configuration': _perceive_balance_configuration(beam, block_mass_u, distance_u),
            'BeamWeight': bw}

def make_cp_space_balance(trial, block_mass_u=None, distance_u=None,
                     com_range_vs=1., com_range_s=1., com_range_m=1., com_range_l=1.,
                     beam_mass_mean=4., beam_mass_sd=None):


    strut, beam = trial
    perc_config = []

    if strut[1] == .25:
        com_range = com_range_vs
    elif strut[1] == .5:
        com_range = com_range_s
    elif strut[1] == 1.:
        com_range = com_range_m
    elif strut[1] == 2.:
        com_range = com_range_l
    else:
        raise Exception('Uh oh - bad strut reading')

    if beam_mass_sd:
        bw = truncNorm(beam_mass_mean, beam_mass_sd, min=.000001)
    else:
        bw = beam_mass_mean

    w = _physics_world(-strut[0], com_range, bw)

    for d, bls in beam:
        totmass = 0
        ## IGNORE BLOCK MASS U - NOT NEEDED BUT WOULD CAUSE COMPLICATIONS
        z = .7122
        dz = .412
        for _ in xrange(bls):
            totmass += 1
            w.add_standard_block(1, -d, z, (196,196,0))
            z += dz

        perc_config.append((totmass, -d))

    perc = {'CoM': (-strut[0], strut[1], com_range),
            'Configuration': perc_config,
            'BeamWeight': bw}

    w.perc = perc

    return w
