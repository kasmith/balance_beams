from __future__ import division
from .common_parts import standard_strut
import numpy as np

def _perceive_shapes_configuration(trial, tower_mass_u = None, block_mass_u = None, distance_u = None):
    perc_config = []
    for b, m, d in trial:
        if b == 'B':
            if block_mass_u:
                m *= np.random.lognormal(0, block_mass_u)
        else:
            if tower_mass_u:
                m *= np.random.lognormal(0, tower_mass_u)
        if distance_u:
            d += np.random.normal(0, distance_u)
        perc_config.append((m,d))
    return perc_config

def perceive_shapes(trial, tower_mass_u = None, block_mass_u = None, distance_u = None, com_range = 1.):
    return {'CoM': standard_strut(com_range),
            'Configuration': _perceive_shapes_configuration(trial, tower_mass_u, block_mass_u, distance_u),
            'BeamWeight': 4.}

