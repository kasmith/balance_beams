from __future__ import division
from .common_parts import standard_strut
import numpy as np
from .physics_helper import PhysicsWorld

def _perceive_materials_configuration(trial, block_mass_u=None, brick_density=2.4, brick_u=None,
                                      iron_density=9., iron_u=None, distance_u=None):
    perc_config = []

    if brick_u:
        u_brick = brick_density * np.random.lognormal(0, brick_u)
    else:
        u_brick = brick_density
    if iron_u:
        u_iron = iron_density * np.random.lognormal(0, iron_u)
    else:
        u_iron = iron_density

    for d, mats in trial:
        totmass = 0
        for _ in range(mats['W']):
            if block_mass_u:
                totmass += np.random.lognormal(0, block_mass_u)
            else:
                totmass += 1
        for _ in range(mats['B']):
            if block_mass_u:
                totmass += u_brick*np.random.lognormal(0, block_mass_u)
            else:
                totmass += u_brick
        for _ in range(mats['I']):
            if block_mass_u:
                totmass += u_iron*np.random.lognormal(0, block_mass_u)
            else:
                totmass += u_iron
        if distance_u:
            d += np.random.normal(0, distance_u)

        perc_config.append((totmass, d))

    return perc_config

def perceive_materials(trial, block_mass_u=None, brick_density=2.4, brick_u=None,
                                      iron_density=9., iron_u=None, distance_u=None, com_range=1.):
    return {'CoM': standard_strut(com_range),
            'Configuration': _perceive_materials_configuration(trial, block_mass_u, brick_density, brick_u,
                                                            iron_density, iron_u, distance_u),
            'BeamWeight': 4.}


# For actual physical simulation
def make_cp_space_materials(trial, block_mass_u=None, brick_density=2.4, brick_u=None,
                                      iron_density=9., iron_u=None, distance_u=None, com_range=1.):

    perc_config = []
    if brick_u:
        u_brick = brick_density * np.random.lognormal(0, brick_u)
    else:
        u_brick = brick_density
    if iron_u:
        u_iron = iron_density * np.random.lognormal(0, iron_u)
    else:
        u_iron = iron_density

    w = _physics_world(0, com_range)

    for d, mats in trial:
        if distance_u:
            d += np.random.normal(0, distance_u)
        totmass = 0

        z = .7122
        dz = .412
        for _ in xrange(mats['I']):
            if block_mass_u:
                m = u_iron*np.random.lognormal(0, block_mass_u)
            else:
                m = u_iron
            w.add_standard_block(m, d, z, (128, 128, 128))
            z += dz
            totmass += m

        for _ in xrange(mats['B']):
            if block_mass_u:
                m = u_brick * np.random.lognormal(0, block_mass_u)
            else:
                m = u_brick
            w.add_standard_block(m, d, z, (255, 0, 0))
            z += dz
            totmass += m

        for _ in xrange(mats['W']):
            if block_mass_u:
                m = np.random.lognormal(0, block_mass_u)
            else:
                m = 1
            w.add_standard_block(m, d, z, (196,196,0))
            z += dz
            totmass += m

        perc_config.append((totmass, d))

    perc = {'CoM': standard_strut(com_range),
            'Configuration': perc_config,
            'BeamWeight': 4.}

    w.perc = perc

    return w
