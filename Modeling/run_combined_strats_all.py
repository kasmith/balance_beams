from models import *
import os
import json
import argparse

# Settings (debug & strategy selection)

strategy_types = {
    'comb': MSP_MOD_TYPES,
    'inc_dist': ALL_MOD_TYPES,
    'base_strats': ['smp', 'sp'],
    'no_phys': make_type_defs('sm'),
    'no_sym': make_type_defs('mp'),
    'no_weight': make_type_defs('sp'),
    'just_sp': ['sp'],
    'just_smp': ['smp'],
    'just_phys': ['p'],
    'swap_dist': ['sdp', 'sp'],
    'add_dist': ['smp', 'sp', 'sdp', 'sdmp', 'smdp'],
    'rules': []
}

for stp in ALL_MOD_TYPES:
    if stp not in ['sp', 'smp']:
        strategy_types['add1_'+stp] = ['sp', 'smp', stp]


parser = argparse.ArgumentParser(description="Run basic strategy models")
parser.add_argument('-s', '--strategy', help="Name of strategy collection",
                    default="comb", choices=strategy_types.keys())
parser.add_argument('-o', '--output', help="Base name of output files",
                    default=None)
parser.add_argument('-d', '--nodebug', help="Turns off debug mode",
                    action="store_true")
parser.add_argument('-i', '--initialization', default=None,
                    type=argparse.FileType('rU'),
                    help="File with initialization parameters")
parser.add_argument('-f', '--fit_type', default="all",
                    help="How to fit (by all, individual, joint)",
                    choices=["all", "joint_strats", "joint_percept",
                             "individual", "none", "joint_halfstrat"])
parser.add_argument('--single_strat', action="store_true",
                    help="Whether to fix individuals to use a single strategy")
parser.add_argument('--random_split', action="store_true",
                    help="Randomly split trials in two instead of using order (for joint_halfstrat)")
parser.add_argument('--use_learnbene', action='store_true',
                    help="Fit the learn_benefit data rather than the regular")
parser.add_argument('--use_geomat', action='store_true',
                    help="Fit the geomat data rather than the regular")

if __name__ == '__main__':
    args = parser.parse_args()
    print ("Arguments:")
    print (args)

    assert not (args.use_learnbene and args.use_geomat),\
        "Cannot set learnbene and geomat!!!"

    if args.use_learnbene:
        fitpars = ['com_u', 'com_range', 'mass_jnd',
                   'dist_jnd', 'distance_u']
        models = {'learnbene': model_clean}
        rule_models = {'learnbene': rules_clean}
        params = {'learnbene': clean_param_names}
        trials = {'learnbene': learnbene_trials}
        empdat = {'learnbene': learnbene_empirical}
    elif args.use_geomat:
        assert args.strategy != 'rules', 'Rules not allowed for geomat'
        fitpars = ['com_u', 'com_range', 'mass_jnd',
                   'dist_jnd', 'distance_u']
        models = {'geomat': model_geomat}
        params = {'geomat': clean_param_names}
        trials = {'geomat': geomat_trials}
        empdat = {'geomat': geomat_empirical}

    else:
        # Set up the normal stuff
        fitpars = ['com_u', 'com_range', 'beam_mass_mean', 'com_range_s',
                   'com_range_vs', 'com_range_m', 'com_range_l', 'mass_jnd',
                   'dist_jnd', 'distance_u', 'tower_mass_u']

        # Define the data to use
        models = {'sh': model_shapes,
                  'mat': model_materials,
                  'bal': model_balance}
        rule_models = {'sh': rules_shapes,
                       'mat': rules_materials,
                       'bal': rules_balance}
        params = {'sh': shapes_param_names,
                  'mat': materials_param_names,
                  'bal': balance_param_names}
        trials = {'sh': shapes_trials,
                  'mat': materials_trials,
                  'bal': balance_trials}
        empdat = {'sh': shapes_empirical,
                  'mat': materials_empirical,
                  'bal': balance_empirical}
        orders = {'sh': shapes_emp_order,
                  'mat': materials_emp_order,
                  'bal': balance_emp_order}


    # Cut out unneeded parameters
    for pnames in params.values():
        safe_remove(pnames, 'beam_mass_sd')

    strat_nm = args.strategy
    onm = args.output

    if strat_nm == 'rules':
        params = {
            'sh': params_mod_2_rules(params['sh'], False),
            'mat': params_mod_2_rules(params['mat'], False),
            'bal': params_mod_2_rules(params['bal'], False)
        }
        models = rule_models # Overwrite the models to use

    if args.nodebug:
        if args.fit_type in ['all', 'joint_strats']:
            nsims = 500
        else:
            nsims = 250
        maxiter = 500
        printiter = 50
        ncore = 30
    else:
        nsims = 10
        maxiter = 10
        printiter = 1
        ncore = 1
        for m in trials.keys():
            trials[m] = {k: trials[m][k] for k in list(trials[m].keys())[:10]}

    if strat_nm == 'rules':
        strategies = ['R1', 'R2', 'R3', 'R3a', 'R4']
    else:
        strategies = strategy_types[strat_nm]

    if args.use_learnbene:
        output_folder = "learn_bene"
    elif args.use_geomat:
        output_folder = "geomat"
    else:
        output_folder = "comb_strats"

    if onm is None:
        if args.fit_type == 'all':
            ofl_base = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    "output", output_folder, strat_nm+"_all")
        elif args.fit_type == 'joint_strats':
            ofl_base = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    "output", output_folder, strat_nm+"_joint_strat")
        elif args.fit_type == 'joint_percept':
            ofl_base = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    "output", output_folder, strat_nm+"_joint_percept")
        elif args.fit_type == "individual":
            ofl_base = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    "output", output_folder, strat_nm+"_individual")
        elif args.fit_type == 'joint_halfstrat':
            ofl_base = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    "output", output_folder, strat_nm+"_halfstrat")
        if not args.nodebug:
            ofl_base += "_DEBUG"
    else:
        ofl_base = onm
    ofl_params = ofl_base + "_params.json"

    if args.initialization is not None:
        init = json.load(args.initialization)
        iperc = init['perceptual']
        istrat = init['strategies']
    else:
        iperc = None
        istrat = None

    print ('starting run', strategies)

    if args.fit_type == "all":
        perc_params, strat_mix, llh = optimize_complete(models, trials, empdat, fitpars,
                                                        params, strategies, iperc,
                                                        istrat, max_iters=maxiter,
                                                        print_iters=printiter,
                                                        n_reps=nsims, n_cores=ncore,
                                                        savestate=ofl_base+"_TMP.json")

        with open(ofl_params, 'w') as paramfl:
            json.dump({'perceptual': perc_params,
                       'strategies': strat_mix,
                       'LLH': llh},
                      paramfl)

        for mod in models.keys():
            pnames = params[mod]
            perc_ps = fill_params(dict([(p, perc_params[p])
                                        for p in perc_params.keys()
                                        if p in pnames]), pnames)
            onm = ofl_base + "_" + mod + ".csv"
            write_all_outcomes(onm, models[mod], empdat[mod], perc_ps, strat_mix,
                               trials[mod], n_sims=nsims)

    elif args.fit_type == "joint_percept":
        perc_params, ind_strats, llh = optimize_joint_perc_complete(models, trials, empdat, fitpars,
                                                        params, strategies, iperc,
                                                        istrat, max_iters=maxiter,
                                                        print_iters=printiter,
                                                        n_reps=nsims, n_cores=ncore,
                                                        savestate=ofl_base+"_TMP.json",
                                                        enforce_single_strat=args.single_strat)
        with open(ofl_params, 'w') as paramfl:
            json.dump({'perceptual': perc_params,
                       'ind_strategies': ind_strats,
                       'LLH': llh},
                      paramfl)

        for mod in models.keys():
            pnames = params[mod]
            perc_ps = fill_params(dict([(p, perc_params[p])
                                        for p in perc_params.keys()
                                        if p in pnames]), pnames)
            onm = ofl_base + "_" + mod + ".csv"
            write_joint_strat_outcomes(onm, models[mod], empdat[mod], perc_ps,
                                       ind_strats, trials[mod], n_sims=nsims)

    elif args.fit_type == "joint_halfstrat":
        # Shuffle the order for randomization tests if requested
        if args.random_split:
            new_ords = dict()
            for m in orders.keys():
                new_ords[m] = dict()
                for w in orders[m].keys():
                    trnms, idxs = zip(*orders[m][w].items())
                    idxs = list(idxs)
                    random.shuffle(idxs)
                    new_ords[m][w] = dict(zip(trnms, idxs))
            orders = new_ords
        perc_params, ind_strats, llh = optimize_joint_strats_byhalf(models, trials, empdat, fitpars,
                                                        params, orders, strategies, iperc,
                                                        istrat, max_iters=maxiter,
                                                        print_iters=printiter,
                                                        n_reps=nsims, n_cores=ncore,
                                                        savestate=ofl_base+"_TMP.json")
        with open(ofl_params, 'w') as paramfl:
            json.dump({'perceptual': perc_params,
                       'ind_strategies': ind_strats,
                       'LLH': llh},
                      paramfl)
        for mod in models.keys():
            pnames = params[mod]
            perc_ps = fill_params(dict([(p, perc_params[p])
                                        for p in perc_params.keys()
                                        if p in pnames]), pnames)
            onm = ofl_base + "_" + mod + ".csv"
            write_half_strat_outcomes(onm, models[mod], empdat[mod], perc_ps,
                                      ind_strats, orders[mod], trials[mod],
                                      n_sims=nsims)


    elif args.fit_type == "joint_strats":
        ind_perc, strat_mix, llh = optimize_joint_strats_complete(models, trials, empdat, fitpars,
                                                        params, strategies, iperc,
                                                        istrat, max_iters=maxiter,
                                                        print_iters=printiter,
                                                        n_reps=nsims, n_cores=ncore,
                                                        savestate=ofl_base+"_TMP.json")
        with open(ofl_params, 'w') as paramfl:
            json.dump({'perceptual': ind_perc,
                       'strategies': strat_mix,
                       'LLH': llh},
                      paramfl)

        for mod in models.keys():
            wids = empdat[mod].keys()
            c_strats = dict([(wid, strat_mix) for wid in wids])
            onm = ofl_base + "_" + mod + ".csv"
            write_ind_outcomes(onm, models[mod], empdat[mod], ind_perc[mod],
                               c_strats, trials[mod], n_sims=nsims)

    elif args.fit_type == "individual":
        ind_perc, ind_strat, llh = optimize_ind(models, trials, empdat, fitpars,
                                                        params, strategies, iperc,
                                                        istrat, max_iters=maxiter,
                                                        print_iters=printiter,
                                                        n_reps=nsims, n_cores=ncore,
                                                        savestate=ofl_base+"_TMP.json",
                                                        enforce_single_strat=args.single_strat)
        with open(ofl_params, 'w') as paramfl:
            json.dump({'perceptual': ind_perc,
                       'strategies': ind_strat,
                       'LLH': llh},
                      paramfl)

        for mod in models.keys():
            wids = empdat[mod].keys()
            onm = ofl_base + "_" + mod + ".csv"
            write_ind_outcomes(onm, models[mod], empdat[mod], ind_perc[mod],
                               ind_strat, trials[mod], n_sims=nsims)
