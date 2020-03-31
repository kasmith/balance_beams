from shapes import shapes_trials, shapes_empirical
from materials import materials_trials, materials_empirical
from balance import balance_trials, balance_empirical
from model_parts import *
from config import *
from optimization import *
from writeout import *
import time

modfn = model_all_balance

paramdict = {
    'block_mass_u' : .00001,
    'distance_u' : 0.098988625,
    'mass_heur_p' : 0.369815231,
    'mass_jnd' : 0.816178932,
    'dist_jnd' : 0.530136131,
    'com_u' : 0.429080942,
    'com_range_vs' : 0.165082296,
    'com_range_s' : 0.229277034,
    'com_range_m' : 0.278916477,
    'com_range_l' : 0.574358269,
    'beam_mass_mean' : 0.5,
    'beam_mass_sd' : 0.5,
    'lapse' : 0.153380576,
    'n_times' : 500,
    'use_full_balance' : True
}

'''
for bmm in [0.25, 0.5, 1., 2., 4., 8.]:
    for bms in [0.001, 0.1, 0.25, 0.5, 1., 2.]:
        paramdict['beam_mass_mean'] = bmm
        paramdict['beam_mass_sd'] = bms
        print bmm, bms
        print get_all_llh(modfn, balance_trials, balance_empirical, **paramdict)
'''
t = time.time()
print (get_all_llh(modfn, balance_trials, balance_empirical, **paramdict))
print ('Done:', (time.time() - t))
print ('')

t = time.time()
print (get_all_llh(modfn, balance_trials, balance_empirical, parallelize=False, **paramdict))
print ('Done:', (time.time() - t))
print ('')
