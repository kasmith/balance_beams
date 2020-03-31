from __future__ import division
import numpy as np
import random

# Strut: position, real width, range width
def standard_strut(com_range):
    return (0, .05, com_range)

def truncNorm(loc=0,scale=1,min=None,max=None):
    r = np.random.normal(loc,scale)
    if min is not None and r < min: return truncNorm(loc,scale,min,max)
    if max is not None and r > max: return truncNorm(loc,scale,min,max)
    return r

def getBeamWeight(mean, sd = None, p_miss = 0, type = 'trunc_norm'):

    assert type in ['trunc_norm', 'exponential', 'probabilistic']

    if type == 'trunc_norm':
        if sd:
            return truncNorm(mean, sd, min=.000001)
        else:
            return mean
    elif type == 'exponential':
        return np.random.exponential(mean)

    elif type == 'probabilistic':
        if random.random() < p_miss:
            return .000001
        else:
            return mean


def make_type_defs(types = 'smdp', calling = ''):
    chs = types + '\t'
    r = []
    for c in chs:
        if c == '\t':
            if calling != '':
                r += [calling]
        else:
            r += make_type_defs(types.replace(c,''), calling+c)
    return r

ALL_MOD_TYPES = make_type_defs()
MSP_MOD_TYPES = make_type_defs('smp')
