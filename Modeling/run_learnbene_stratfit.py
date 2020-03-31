from __future__ import division, print_function
from models import *
from OptimTools import async_map
import os
import json
import argparse
from scipy.optimize import minimize, minimize_scalar

strategy_types = {
    'comb': MSP_MOD_TYPES,
    'inc_dist': ALL_MOD_TYPES,
    'base_strats': ['smp', 'sp'],
    'no_phys': make_type_defs('sm'),
    'no_sym': make_type_defs('mp'),
    'no_weight': make_type_defs('sp'),
    'just_sp': ['sp'],
    'just_smp': ['smp'],
    'rules': []
}

for stp in MSP_MOD_TYPES:
    if stp not in ['sp', 'smp']:
        strategy_types['add1_'+stp] = ['sp', 'smp', stp]

parser = argparse.ArgumentParser("Fit only strategies from perceptual parameters")
parser.add_argument('-i', '--initialization', required=True,
                    type=argparse.FileType('rU'),
                    help="File with initialization parameters")
parser.add_argument('-o', '--output', help="Base name of output files",
                    default=None)
parser.add_argument('-d', '--nodebug', help="Turns off debug mode",
                    action="store_true")
parser.add_argument('-s', '--strategy', help="Name of strategy collection",
                    default="base_strats", choices=strategy_types.keys())

if __name__ == '__main__':
    args = parser.parse_args()
    print ("Arguments:")
    print (args)

    fitpars = ['com_u', 'com_range', 'mass_jnd',
               'dist_jnd', 'distance_u']
    model = model_clean
    rule_model = rules_clean
    params = clean_param_names
    trials = learnbene_trials
    empdat = learnbene_empirical

    # Cut out unneeded parameters
    safe_remove(params, 'beam_mass_sd')
    safe_remove(params, 'block_mass_u')

    strat_nm = args.strategy
    onm = args.output

    if args.nodebug:
        nsims = 1000
        maxiter = 500
        ncore = 4
    else:
        nsims = 10
        maxiter = 10
        ncore = 1
        trials = {k: trials[k] for k in trials.keys()[:4]}

    strategies = strategy_types[strat_nm]
    output_folder = "learn_bene"
    if onm is None:
        ofl_base = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "output", output_folder, strat_nm+"_stratonly")
        if not args.nodebug:
            ofl_base += "_DEBUG"

    init = json.load(args.initialization)
    iperc = init['perceptual']
    perc_params = fill_params(iperc, params)

    istrat = init['strategies']

    print("Starting run:", strategies)

    wids = empdat.keys()
    fitstrats = dict()

    # Make the set of perceptual simulations
    def _get_trial_model(trial_name):
        trial = trials[trial_name]
        return model(trial, strategies, n_times=nsims, **perc_params)

    if ncore == 1:
        tms = [_get_trial_model(tr) for tr in trials.keys()]
    else:
        tms = async_map(_get_trial_model, trials.keys(), ncpu=ncore)

    trial_outcomes = dict(zip(trials.keys(), tms))

    print("Percepts / simulations fit")

    # Optimize the strategies individually
    def _fit_ind_strats(wid):
        dat = empdat[wid]
        def _get_llh(p):
            if any([i > 10000 for i in p]):
                return 99999999999
            mix_dict = transform_strat_strengths(p, strategies)
            tr_names = list(set(dat.keys()) & set(trial_outcomes.keys()))
            return -sum([get_ind_trial_llh_from_mix(trial_outcomes[tr], mix_dict, dat[tr])
                         for tr in tr_names])
        o = minimize(_get_llh, revert_strat_strengths(istrat, strategies))
        print("Fit", wid)
        return o.x, o.fun

    if ncore == 1:
        wfits = [_fit_ind_strats(w) for w in wids]
    else:
        wfits = async_map(_fit_ind_strats, wids)

    # Fill this out
    tot_llh = 0
    ind_strats = dict()
    for w, f in zip(wids, wfits):
        rawps, llh = f
        ind_strats[w] = {
            'llh': llh,
            'strat_params': transform_strat_strengths(rawps, strategies)
        }
        tot_llh += llh

    # Write out the parameters
    param_dict = {
        "perceptual": perc_params,
        "ind_strategies": ind_strats,
        "LLH": tot_llh
    }

    json.dump(param_dict, open(ofl_base + "_params.json", 'w'))

    # Write out the data
    onm = ofl_base + "_learnbene.csv"
    write_joint_strat_outcomes(onm, model, empdat, perc_params, ind_strats,
                               trials, n_sims=nsims)
    print("Done")
