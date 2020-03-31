from models import *
from OptimTools import async_map
import argparse, os, copy, random
import numpy as np
import json


shapes_package = [model_all_shapes, shapes_trials, shapes_empirical, shapes_param_names, rules_all_shapes_full]
materials_package = [model_all_materials, materials_trials, materials_empirical, materials_param_names, rules_all_materials_full]
balance_package = [model_all_balance, balance_trials, balance_empirical, balance_param_names, rules_all_balance_full]
geomat_package = [model_all_materials, geomat_trials, geomat_empirical, materials_param_names, rules_all_materials_full]
combined_package = [model_all_combined, combined_trials, combined_empirical, combined_param_names, rules_all_combined_full]


#modtypes = ALL_MOD_TYPES
## Limit modtypes to make fits reasonable -- just reduced models and re-ordering full model
modtypes = ['msp', 'mdsp', 'mp', 'sp', 'ms', 'm', 's', 'p',
            'mps', 'spm', 'smp', 'pms', 'psm']



def split_trials(all_trials, write_fit_flnm, i=0):
    n_trials = len(all_trials)
    trs = copy.copy(all_trials.keys())
    random.shuffle(trs)
    fittrs = copy.copy(trs[:int(n_trials / 2)])
    cvtrs = copy.copy(trs[int(n_trials / 2):])

    with open(write_fit_flnm,'a') as ofl:
        ofl.write(str(i))
        for t in fittrs:
            ofl.write(',' + t)
        ofl.write('\n')

    return fittrs, cvtrs

def safe_remove(lst, objs):
    if not hasattr(objs, '__iter__'):
        objs = [objs]
    for obj in objs:
        try:
            lst.remove(obj)
        except ValueError:
            pass

def append_to(flnm, to_write):
    with open(flnm,'a') as wfl:
        wfl.write(to_write)

def do_add_ind_cv(model_func, trials, empdat, parameters, initialization = None, out_dir = '.', msprt_samples = 1,
                  n_repeats = 2, n_samples = 10, max_iters = 5, init_grid = 2, print_iter = 5, beam_type = 'exponential'):

    WIDS = empdat.keys()

    # Set up files
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    flnm_cv_trials = os.path.join(out_dir, "fit_trials.csv")
    flnm_summary = os.path.join(out_dir, 'cv_summary.txt')
    flnm_baseline = os.path.join(out_dir, 'baseline_params.csv')
    flnm_llhs = os.path.join(out_dir, 'cv_llh_grid.csv')
    flnms_params = dict([(p,os.path.join(out_dir, 'param_' + p + '.csv')) for p in parameters])

    wfl = open(flnm_cv_trials, 'w')
    wfl.close()

    with open(flnm_baseline, 'w') as wfl:
        wfl.write('Order')
        for p in parameters:
            wfl.write(',' + p)
        wfl.write(',Fit_LLH,CV_LLH\n')

    with open(flnm_llhs, 'w') as wfl:
        wfl.write('Order,Baseline')
        for p in parameters:
            wfl.write(',' + p)
        wfl.write('\n')

    for fp in flnms_params.values():
        with open(fp, 'w') as wfl:
            wfl.write('Order')
            for w in WIDS:
                wfl.write(',' + w)
            wfl.write(',Fit_LLH,CV_LLH\n')

    with open(flnm_summary,'w') as wfl:
        wfl.write('File setup complete\n\n---------------------\n\n')

    # Repeat the cross-validation a number of times
    for rep_n in range(n_repeats):

        append_to(flnm_summary, "CV set " + str(rep_n) + '\n\n')

        # Pull apart the trials
        fittrs, cvtrs = split_trials(trials, flnm_cv_trials, rep_n)
        fit_trials = dict([(t, trials[t]) for t in fittrs])
        cv_trials = dict([(t, trials[t]) for t in cvtrs])

        # Fit the baseline model
        base_params, base_fit_llh = optimize_all(model_func, fit_trials, empdat, parameters, parameters, initialization,
                                                 None, msprt_samples, max_iters, print_iter, n_samples, True, 'smp', beam_type)

        base_fit_llh *= -1 # For turning to negative

        base_cv_llh = get_all_llh(model_func, cv_trials, empdat, n_times = n_samples, use_full_balance = True,
                                  msprt_samples = msprt_samples, parallelize=True, type='smp', beam_type=beam_type,
                                  **base_params)

        with open(flnm_baseline, 'a') as wfl:
            wfl.write(str(rep_n))
            for p in parameters:
                wfl.write(',' + str(base_params[p]))
            wfl.write(',' + str(base_fit_llh) + ',' + str(base_cv_llh) + '\n')

        append_to(flnm_llhs, str(rep_n) + ',' + str(base_cv_llh))
        append_to(flnm_summary, "Baseline:\nFit LLH: " + str(base_fit_llh) + "\nCV LLH: " + str(base_cv_llh) + "\n\n")

        # Go through and fit each of the individual parameter variation models
        for p in parameters:
            popts, fitllhs = optimize_single_ind_param(model_func, fit_trials, empdat, p, parameters, base_params,
                                                        msprt_samples=msprt_samples, max_iters=max_iters,
                                                        i_range = init_grid, print_iters=print_iter, n_reps=n_samples,
                                                        use_full_balance=True, modtype='smp', beam_type=beam_type)

            tot_fit_llh = -sum(fitllhs.values())

            def ind_cv_llh(w):
                params = copy.deepcopy(base_params)
                params[p] = popts[w]
                edat = empdat[w]
                return get_ind_llh(model_func, cv_trials, edat, fit_all_trials = False, n_times = n_samples,
                                   use_full_balance = False, msprt_samples=msprt_samples, parallelize = False,
                                   type = 'smp', beam_type=beam_type, **params)

            cvllhs = async_map(ind_cv_llh, WIDS)
            tot_cv_llh = sum(cvllhs)

            with open(flnms_params[p], 'a') as wfl:
                wfl.write(str(rep_n))
                for w in WIDS:
                    wfl.write(',' + str(popts[w]))
                wfl.write(',' + str(tot_fit_llh) + ',' + str(tot_cv_llh) + '\n')

            append_to(flnm_llhs, ',' + str(tot_cv_llh))
            append_to(flnm_summary, p + ':\nFit LLH: ' + str(tot_fit_llh) + '\nCV LLH: ' + str(tot_cv_llh) + '\n\n')

        append_to(flnm_llhs,'\n')
        append_to(flnm_summary,'CV set ' + str(rep_n) + ' done\n\n------------------\n\n')


def do_rules_extension_cv(model_func, rules_model_func, trials, empdat, fit_param, param_names, initialization_mod,
                          initialization_rules, out_dir = '.', n_repeats = 2, n_samples = 10, max_iters = 5,
                          print_iter = 5, beam_type = 'exponential'):

    WIDS = empdat.keys()

    rules_param_names = params_mod_2_rules(copy.copy(param_names), False)

    mod_params = {}
    for p in param_names:
        if p != fit_param:
            mod_params[p] = initialization_mod[p]
    rule_params = {}
    for p in rules_param_names:
        rule_params[p] = initialization_rules[p]

    # Set up files
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    flnm_cv_trials = os.path.join(out_dir, "fit_trials.csv")
    flnm_summary = os.path.join(out_dir, 'cv_summary.txt')
    flnm_fits = os.path.join(out_dir, 'subj_fit_grid.csv')
    flnm_base_params = os.path.join(out_dir, 'baseline_mod_params.json')
    flnm_rule_params = os.path.join(out_dir, 'baseline_rule_params.json')

    wfl = open(flnm_cv_trials, 'w')
    wfl.close()

    with open(flnm_fits, 'w') as wfl:
        wfl.write('WID,Order,ModParam,ModFitLLH,ModCVLLH,RuleChoice,RuleFitLLH,RuleCVLLH\n')

    with open(flnm_base_params, 'w') as bp:
        json.dump(initialization_mod, bp)
    with open(flnm_rule_params, 'w') as bp:
        json.dump(initialization_rules, bp)

    with open(flnm_summary,'w') as wfl:
        wfl.write('File setup complete\n\n---------------------\n\n')

    # Fit the predictions for all trials by rule (caching)
    #def _get_rule_outcome(trnm):
    #    tr = trials[trnm]
    #    return rules_model_func(tr, n_times=n_samples, beam_type=beam_type,parallelize=False, **rule_params)
    #
    #rule_outcomes_raw = map(_get_rule_outcome, trials.keys())
    #rule_outcomes = ((tr, oc) for tr, oc in zip(trials.keys(), rule_outcomes_raw))

    rule_outcomes_raw = rules_model_func(trials.values(), n_times=n_samples, beam_type=beam_type,parallelize=False, **rule_params)
    rule_outcomes = dict([tr, ro[0]] for tr, ro in zip(trials.keys(), rule_outcomes_raw))

    def _get_rule_llh_from_emp(wid, trnames, rule_names=['R1','R2','R3','R3a','R4']):
        emp = empdat[wid]
        rule_llhs = dict([(rnm, 0) for rnm in rule_names])
        for tr in trnames:
            e = emp[tr]
            ocm = rule_outcomes[tr]
            for rnm in rule_names:
                rule_llhs[rnm] += np.log(ocm[rnm][e])
        return rule_llhs


    # Repeat the cross-validation a number of times
    for rep_n in range(n_repeats):

        append_to(flnm_summary, "CV set " + str(rep_n) + '\n\n')

        # Pull apart the trials
        fittrs, cvtrs = split_trials(trials, flnm_cv_trials, rep_n)
        fit_trials = dict([(t, trials[t]) for t in fittrs])
        cv_trials = dict([(t, trials[t]) for t in cvtrs])

        # Fit the individual model parameter
        mod_fit_params, mod_fit_llhs = optimize_single_ind_param(model_func, fit_trials, empdat, fit_param, param_names,
                                                         mod_params, max_iters=max_iters, n_reps=n_samples,
                                                         beam_type=beam_type,print_iters=print_iter, parallelize=False)
        tot_fit_llhs = sum(mod_fit_llhs.values())
        append_to(flnm_summary, 'Model fit; total LLH: ' + str(tot_fit_llhs) + '\n')

        # Cross-validate the individual model parameters
        mod_cv_llhs = {}
        tot_cv_llh = 0
        for w in WIDS:
            fpar = mod_fit_params[w]
            udict = copy.copy(mod_params)
            udict[fit_param] = fpar
            llh = get_ind_llh(model_func, cv_trials, empdat[w], fit_all_trials=False, n_times=n_samples,
                              beam_type=beam_type, **udict)
            mod_cv_llhs[w] = llh
            tot_cv_llh += llh

        append_to(flnm_summary, 'Model cross-validated; total LLH: ' + str(tot_cv_llh) + '\n')

        # Fit the individual rules
        def _get_fit_rule(wid):
            llhs = _get_rule_llh_from_emp(wid, fittrs)
            best_llh = -99999999999999999
            best_rule = 'NA'
            for r,llh in llhs.items():
                if llh > best_llh:
                    best_llh = llh
                    best_rule = r
            return best_rule, best_llh

        fit_rules = {}
        fit_rule_llhs = {}
        tot_rule_fit_llhs = 0
        for w in WIDS:
            r, llh = _get_fit_rule(w)
            fit_rules[w] = r
            fit_rule_llhs[w] = llh
            tot_rule_fit_llhs += llh
        append_to(flnm_summary, 'Rules fit; total LLH: ' + str(tot_rule_fit_llhs) + '\n')

        # Cross_validate the rules
        cv_rule_llhs = {}
        tot_rule_cv_llhs = 0
        for w in WIDS:
            r = fit_rules[w]
            llh = _get_rule_llh_from_emp(w, cvtrs)[r]
            cv_rule_llhs[w] = llh
            tot_rule_cv_llhs += llh
        append_to(flnm_summary, 'Rules cross-validated; total LLH: ' + str(tot_rule_cv_llhs) + '\n')

        # Write to the file
        with open(flnm_fits, 'a') as wfl:
            for w in WIDS:
                wfl.write(w + ',' + str(rep_n) + ',' + str(mod_fit_params[w]) + ',' + str(-mod_fit_llhs[w][0]) + ',')
                wfl.write(str(mod_cv_llhs[w]) + ',' + str(fit_rules[w]) + ',' + str(fit_rule_llhs[w]) + ',')
                wfl.write(str(cv_rule_llhs[w]) + '\n')
        append_to(flnm_summary, 'CV set ' + str(rep_n) + ' done\n\n------------------\n\n')


def do_mod_type_cv(model_func, trials, empdat, parameters, initialization = None, out_dir = '.', msprt_samples = 1,
                  n_repeats = 2, n_samples = 10, max_iters = 5, print_iter = 5):
    WIDS = empdat.keys()

    # Set up parameters by model type
    def getparamsbymt(mod):
        pnames = copy.copy(parameters)
        if 'p' not in mod:
            safe_remove(pnames, ['com_u', 'beam_mass_mean'])
        if 's' not in mod:
            safe_remove(pnames, 'dist_jnd')
        if 'm' not in mod:
            safe_remove(pnames, 'mass_heur_p')
        if 'd' not in mod:
            safe_remove(pnames, 'dist_heur_p')
        if mod in ['p', 'd', 'pd', 'dp']:
            safe_remove(pnames, 'mass_jnd')
        if mod == 's':
            safe_remove(pnames, ['com_range', 'com_range_vs', 'com_range_s', 'com_range_m', 'com_range_l'])
        return pnames

    param_dict = dict([(m, getparamsbymt(m)) for m in modtypes])

    # Set up files
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    flnm_cv_trials = os.path.join(out_dir, "fit_trials.csv")
    flnm_summary = os.path.join(out_dir, 'cv_summary.txt')
    flnm_llhs = os.path.join(out_dir, 'cv_llh_grid.csv')
    flnms_mods = dict([(m, os.path.join(out_dir, 'mod_' + m + '.csv')) for m in modtypes])

    wfl = open(flnm_cv_trials, 'w')
    wfl.close()

    with open(flnm_llhs, 'w') as wfl:
        wfl.write('Order')
        for m in modtypes:
            wfl.write(',' + m)
        wfl.write('\n')

    for m, fm in flnms_mods.items():
        with open(fm, 'w') as wfl:
            wfl.write('Order')
            for p in param_dict[m]:
                wfl.write(',' + p)
            wfl.write(',Fit_LLH,CV_LLH\n')

    with open(flnm_summary, 'w') as wfl:
        wfl.write('File setup complete\n\n---------------------\n\n')

    # Repeat the cross-validation a number of times
    for rep_n in range(n_repeats):
        append_to(flnm_summary, "CV set " + str(rep_n) + '\n\n')
        append_to(flnm_llhs, str(rep_n))

        # Pull apart the trials
        fittrs, cvtrs = split_trials(trials, flnm_cv_trials, rep_n)
        fit_trials = dict([(t, trials[t]) for t in fittrs])
        cv_trials = dict([(t, trials[t]) for t in cvtrs])

        for m in modtypes:

            parameters = param_dict[m]

            popts, fitllhs = optimize_all(model_func, fit_trials, empdat, parameters, parameters, initialization,
                                                       msprt_samples=msprt_samples, max_iters=max_iters,
                                                       print_iters=print_iter, n_reps=n_samples,
                                                       use_full_balance=True, modtype=m, use_mass_ratio=False)

            tot_fit_llh = -fitllhs

            tot_cv_llh = get_all_llh(model_func, cv_trials, empdat, msprt_samples=msprt_samples, n_times=n_samples,
                                 use_full_balance=True, type=m, use_mass_ratio=False, **popts)

            with open(flnms_mods[m], 'a') as wfl:
                wfl.write(str(rep_n))
                for p in parameters:
                    wfl.write(',' + str(popts[p]))
                wfl.write(',' + str(tot_fit_llh) + ',' + str(tot_cv_llh) + '\n')

            append_to(flnm_llhs, ',' + str(tot_cv_llh))
            append_to(flnm_summary, m + ':\nFit LLH: ' + str(tot_fit_llh) + '\nCV LLH: ' + str(tot_cv_llh) + '\n\n')

        append_to(flnm_llhs, '\n')
        append_to(flnm_summary, 'CV set ' + str(rep_n) + ' done\n\n------------------\n\n')

def do_mod_type_cv_complete(model_func_dict, trials_dict, empdat_dict, fitpars, parameters_dict, initialization = None,
                            beam_type = 'exponential', out_dir = '.', msprt_samples = 1,
                            n_repeats = 2, n_samples = 10, max_iters = 5, print_iter = 5, cv_trials = None):
    mods = model_func_dict.keys()
    assert mods == trials_dict.keys() and mods == empdat_dict.keys() and mods == parameters_dict.keys()

    # Set up parameters by model type
    def getparamsbymt(mod_type, parameters):
        pnames = copy.copy(parameters)
        if 'p' not in mod_type:
            safe_remove(pnames, ['com_u', 'beam_mass_mean'])
        if 's' not in mod_type:
            safe_remove(pnames, 'dist_jnd')
        if 'm' not in mod_type:
            safe_remove(pnames, 'mass_heur_p')
        if 'd' not in mod_type:
            safe_remove(pnames, 'dist_heur_p')
        if mod_type in ['p', 'd', 'pd', 'dp']:
            safe_remove(pnames, 'mass_jnd')
        if mod_type == 's':
            safe_remove(pnames, ['com_range', 'com_range_vs', 'com_range_s', 'com_range_m', 'com_range_l'])
        return pnames

    param_dict = dict()
    fitpars_dict = dict()
    for mt in modtypes:
        param_dict[mt] = dict()
        for m in mods:
            pd = parameters_dict[m]
            param_dict[mt][m] = getparamsbymt(mt,pd)
        fitpars_dict[mt] = list(set([pname for plist in parameters_dict.values() for pname in plist if pname in fitpars]))


    # Set up files
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    flnm_cv_trials = os.path.join(out_dir, "fit_trials.csv")
    flnm_summary = os.path.join(out_dir, 'cv_summary.txt')
    flnm_llhs = os.path.join(out_dir, 'cv_llh_grid.csv')
    flnms_mods = dict([(m, os.path.join(out_dir, 'mod_' + m + '.csv')) for m in modtypes])

    if cv_trials:
        cv_trial_dict = dict()
        with open(flnm_cv_trials,'w') as wfl:
            for ln in cv_trials:
                l = ln.strip('\n').split(',')
                idx = int(l[0])
                mod = l[1]
                trs = l[2:]
                if idx not in cv_trial_dict.keys():
                    cv_trial_dict[idx] = dict()
                cv_trial_dict[idx][mod] = trs
                if idx < n_repeats:
                    wfl.write(ln)
    else:
        wfl = open(flnm_cv_trials, 'w')
        wfl.close()


    with open(flnm_llhs, 'w') as wfl:
        wfl.write('Order')
        for m in modtypes:
            wfl.write(',' + m)
        wfl.write('\n')

    for m, fm in flnms_mods.items():
        with open(fm, 'w') as wfl:
            wfl.write('Order')
            for p in fitpars_dict[m]:
                wfl.write(',' + p)
            wfl.write(',Fit_LLH,CV_LLH\n')

    with open(flnm_summary, 'w') as wfl:
        wfl.write('File setup complete\n\n---------------------\n\n')

    # Repeat the cross-validation a number of times
    for rep_n in range(n_repeats):
        append_to(flnm_summary, "CV set " + str(rep_n) + '\n\n')
        append_to(flnm_llhs, str(rep_n))

        # Pull apart the trials

        fit_trials_bymod = dict()
        cv_trials_bymod = dict()
        for m in mods:
            trials = trials_dict[m]
            if cv_trials and rep_n in cv_trial_dict.keys():
                fittrs = cv_trial_dict[rep_n][m]
                cvtrs = [tr for tr in trials if tr not in fittrs]
            else:
                fittrs, cvtrs = split_trials(trials, flnm_cv_trials, str(rep_n) + ',' + m)
            fit_trials_bymod[m] = dict([(t, trials[t]) for t in fittrs])
            cv_trials_bymod[m] = dict([(t, trials[t]) for t in cvtrs])

        for mt in modtypes:

            popts, fitllhs = optimize_complete(model_func_dict, fit_trials_bymod, empdat_dict, fitpars_dict[mt],
                                               param_dict[mt], initialization, msprt_samples=msprt_samples,
                                               max_iters=max_iters, n_reps=n_samples,
                                               modtype=mt, beam_type=beam_type)

            tot_fit_llh = -fitllhs

            tot_cv_llh = 0
            for m in mods:
                ps = fill_params(dict([(p, popts[p]) for p in popts.keys() if p in param_dict[mt][m]]), param_dict[mt][m])
                tot_cv_llh += get_all_llh(model_func_dict[m], cv_trials_bymod[m], empdat_dict[m],
                                          msprt_samples=msprt_samples, n_times=n_samples,
                                          type=mt, beam_type=beam_type, **ps)

            with open(flnms_mods[mt], 'a') as wfl:
                wfl.write(str(rep_n))
                print fitpars_dict[mt]
                for p in fitpars_dict[mt]:
                    wfl.write(',' + str(popts[p]))
                wfl.write(',' + str(tot_fit_llh) + ',' + str(tot_cv_llh) + '\n')

            append_to(flnm_llhs, ',' + str(tot_cv_llh))
            append_to(flnm_summary, mt + ':\nFit LLH: ' + str(tot_fit_llh) + '\nCV LLH: ' + str(tot_cv_llh) + '\n\n')

        append_to(flnm_llhs, '\n')
        append_to(flnm_summary, 'CV set ' + str(rep_n) + ' done\n\n------------------\n\n')

def do_comp_rules(fitpars, initialization = None, beam_type = 'exponential', out_dir = '.', msprt_samples = 1,
                  n_repeats = 2, n_samples = 10, max_iters = 5, print_iter = 5, cv_trials = None, parallelize = True,
                  remove_pars = None):
    model_func_dict = {'sh': model_all_shapes, 'mat': model_all_materials, 'bal': model_all_balance}
    rules_func_dict = {'sh': rules_all_shapes_full, 'mat': rules_all_materials_full, 'bal': rules_all_balance_full}
    trials_dict = {'sh': shapes_trials, 'mat': materials_trials, 'bal': balance_trials}
    empdat_dict = {'sh': shapes_empirical, 'mat': materials_empirical, 'bal': balance_empirical}
    parameters_dict = {'sh': shapes_param_names, 'mat': materials_param_names, 'bal': balance_param_names}

    mods = model_func_dict.keys()

    WIDS = []
    wid_lookup = dict()
    for m in mods:
        WIDS += empdat_dict[m].keys()
        for w in empdat_dict[m].keys():
            wid_lookup[w] = m
        if remove_pars is not None:
            for p in remove_pars:
                safe_remove(parameters_dict[m],p)
        if beam_type != 'trunc_norm':
            safe_remove(parameters_dict[m], 'beam_mass_sd')
        if beam_type == 'probabilistic':
            parameters_dict[m].append('beam_use_prob')

    parameters_rule_dict = dict([(m, params_mod_2_rules(copy.copy(ps), False)) for m, ps in parameters_dict.items()])
    all_params_base = list(set([pname for plist in parameters_dict.values() for pname in plist]))
    all_params_rule = list(set([pname for plist in parameters_rule_dict.values() for pname in plist]))
    fitpars_base = [p for p in fitpars if p in all_params_base]
    fitpars_rule = [p for p in fitpars if p in all_params_rule]
    fit_agg_base = [p for p in fitpars_base if p != 'mass_heur_p']
    fixed_base = [p for p in all_params_base if p not in fitpars_base]
    fixed_rules = [p for p in all_params_rule if p not in fitpars_rule]

    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    # Set up files
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    flnm_cv_trials = os.path.join(out_dir, "fit_trials.csv")
    flnm_summary = os.path.join(out_dir, 'cv_summary.txt')
    flnm_llhs = os.path.join(out_dir, 'cv_llh_grid.csv')
    flnm_rule_chs = os.path.join(out_dir, 'rule_choices.csv')
    flnm_base_params = os.path.join(out_dir, 'base_params.csv')
    flnm_rule_params = os.path.join(out_dir, 'rule_params.csv')
    flnm_weight_params = os.path.join(out_dir, 'weightp_params.csv')

    if cv_trials:
        cv_trial_dict = dict()
        with open(flnm_cv_trials, 'w') as wfl:
            for ln in cv_trials:
                l = ln.strip('\n').split(',')
                idx = int(l[0])
                mod = l[1]
                trs = l[2:]
                if idx not in cv_trial_dict.keys():
                    cv_trial_dict[idx] = dict()
                cv_trial_dict[idx][mod] = trs
                if idx < n_repeats:
                    wfl.write(ln)
    else:
        wfl = open(flnm_cv_trials, 'w')
        wfl.close()

    with open(flnm_llhs,'w') as fl_llh, open(flnm_rule_chs,'w') as fl_rule, open(flnm_weight_params, 'w') as fl_wt:
        fl_llh.write('Order,Model')
        fl_rule.write('Order')
        fl_wt.write('Order')
        for w in WIDS:
            fl_llh.write(',' + w)
            fl_rule.write(',' + w)
            fl_wt.write(',' + w)
        fl_llh.write(',Total_CV,Total_Fit\n')
        fl_rule.write('\n')
        fl_wt.write('\n')

    with open(flnm_base_params,'w') as fl_base, open(flnm_rule_params,'w') as fl_rule:
        fl_base.write('Order')
        fl_rule.write('Order')
        for p in fitpars_base:
            if p != 'mass_heur_p':
                fl_base.write(',' + p)
        for p in fitpars_rule:
            if p != 'mass_heur_p':
                fl_rule.write(',' + p)
        fl_base.write('\n')
        fl_rule.write('\n')

    with open(flnm_summary, 'w') as wfl:
        wfl.write('File setup complete\n\n---------------------\n\n')

    # Repeat cross-val a number of times
    for rep_n in range(n_repeats):
        append_to(flnm_summary, "CV set " + str(rep_n) + '\n\n')

        # Pull apart the trials
        fit_trials_bymod = dict()
        cv_trials_bymod = dict()
        for m in mods:
            trials = trials_dict[m]
            if cv_trials and rep_n in cv_trial_dict.keys():
                fittrs = cv_trial_dict[rep_n][m]
                cvtrs = [tr for tr in trials if tr not in fittrs]
            else:
                fittrs, cvtrs = split_trials(trials, flnm_cv_trials, str(rep_n) + ',' + m)
            fit_trials_bymod[m] = dict([(t, trials[t]) for t in fittrs])
            cv_trials_bymod[m] = dict([(t, trials[t]) for t in cvtrs])


        # Run the individual fits for the baseline model
        optdict = optimize_one_joint_complete(model_func_dict, fit_trials_bymod, empdat_dict, "mass_heur_p",
                                              fit_agg_base, parameters_dict, initialization, fixed_base, msprt_samples,
                                              max_iters, print_iter, 2, n_samples, True, 'msp', beam_type, True, 10,
                                              parallelize)

        # Cross-validate
        def cv_base(w):
            m = wid_lookup[w]
            mod = model_func_dict[m]
            trs = cv_trials_bymod[m]
            edat = empdat_dict[m][w]
            odp = optdict[w][0]
            useps = [p for p in odp.keys() if p in parameters_dict[m]]
            ps = dict([(p, odp[p]) for p in useps])
            #fill_params(ps)
            return get_ind_llh(mod, trs, edat, beam_type=beam_type,type='msp',parallelize=False,
                               msprt_samples=msprt_samples, n_times=n_samples, **ps)

        cvs = mapfn(cv_base, WIDS)
        cv_base_dict = dict([(w, c) for w, c in zip(WIDS, cvs)])

        firstwrite = True
        tot_fit_llh = 0
        tot_cv_llh = 0
        with open(flnm_llhs, 'a') as fl_llh, open(flnm_base_params,'a') as fl_base, \
                open(flnm_weight_params, 'a') as fl_wt:
            fl_llh.write(str(rep_n) + ',Sim')
            fl_wt.write(str(rep_n))
            for w in WIDS:
                params, fit_llh = optdict[w]
                cv_llh = cv_base_dict[w]
                if firstwrite:
                    fl_base.write(str(rep_n))
                    for p in fitpars_base:
                        if p != 'mass_heur_p':
                            fl_base.write(',' + str(params[p]))
                    fl_base.write('\n')
                    firstwrite = False
                tot_fit_llh += fit_llh
                tot_cv_llh += cv_llh
                fl_llh.write(',' + str(cv_llh))
                fl_wt.write(',' + str(params['mass_heur_p']))
            fl_llh.write(',' + str(tot_cv_llh) + ',' + str(tot_fit_llh) + '\n')
            fl_wt.write('\n')

        append_to(flnm_summary, "Base fit LLH: " + str(tot_fit_llh) + "\nBase CV LLH: " + str(tot_cv_llh) + '\n\n')


        # Run the individual rule applications
        rule_ret = optimize_for_individual_rules(fit_trials_bymod['sh'], fit_trials_bymod['mat'], fit_trials_bymod['bal'],
                                                 empdat_dict['sh'], empdat_dict['mat'], empdat_dict['bal'],
                                                 fitpars_rule, parameters_rule_dict, initialization, fixed_rules,
                                                 max_iters, print_iter, n_samples, beam_type, ['1','2','3', '3a','4'])
        rule_params, rule_assignments, rule_fit_llh = rule_ret

        with open(flnm_rule_params, 'a') as fl_rp:
            fl_rp.write(str(rep_n))
            for p in fitpars_rule:
                fl_rp.write(',' + str(rule_params[p]))
            fl_rp.write('\n')

        with open(flnm_rule_chs, 'a') as fl_rch:
            fl_rch.write(str(rep_n))
            for w in WIDS:
                fl_rch.write(',' + str(rule_assignments[w]))
            fl_rch.write('\n')

        append_to(flnm_summary, "Rule fit LLH: " + str(rule_fit_llh) + '\n')

        # Calculate the predictions on the cross-val trials
        def get_rule_pred(m):
            mod = rules_func_dict[m]
            trs = cv_trials_bymod[m]
            useps = [p for p in rule_params.keys() if p in parameters_rule_dict[m]]
            ps = dict([(p, rule_params[p]) for p in useps])
            utrs = [trs[t] for t in trs.keys()]
            rets = mod(utrs, beam_type=beam_type, parallelize=parallelize, n_times=n_samples, **ps)

            rules_ret = dict()
            for rule in ['1','2','3','3a','4']:
                rules_ret[rule] = dict()
                for ret, trnm in zip(rets, trs.keys()):
                    ch = dict()
                    for c in ['L','B','R']:
                        ch[c] = ret[0]['R'+rule][c]
                    rules_ret[rule][trnm] = ch

            return rules_ret

        rule_preds = dict([(m, get_rule_pred(m)) for m in mods])

        # Figure out individual LLHs
        rule_cv_llh = 0
        with open(flnm_llhs, 'a') as fl_llh:
            fl_llh.write(str(rep_n) + ',Rule')
            for w in WIDS:
                m = wid_lookup[w]
                rule = rule_assignments[w]
                edat = empdat_dict[m][w]
                modpreds = rule_preds[m][rule]
                ind_cv_llh = 0
                for tr in cv_trials_bymod[m].keys():
                    if tr in edat.keys():
                        ind_cv_llh += np.log(modpreds[tr][edat[tr]])
                rule_cv_llh += ind_cv_llh
                fl_llh.write(',' + str(ind_cv_llh))
            fl_llh.write(',' + str(rule_cv_llh) + ',' + str(rule_fit_llh) + '\n')

        append_to(flnm_summary, "Rule CV LLH: " + str(rule_cv_llh) + '\n\n')

        append_to(flnm_summary, 'CV set ' + str(rep_n) + ' done\n\n------------------\n\n')




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Run model cross-validation for balance beams models.")
    parser.add_argument('-o', '--output', help="Destination folder for model output", required=True)
    parser.add_argument('-t', '--type', help="The type of beams model to fit",
                        default='shapes', choices=['shapes','materials','balance','complete','combined'])
    parser.add_argument('-i', '--init_parameters', help="Filename for parameter initialization (optional)",
                        default=None, type=argparse.FileType('rU'))
    parser.add_argument('-c', '--mod_class', help="Whether to fit individual parameter changes or model type changes",
                        default='add_ind', choices=['add_ind','model_type','comp_rules','rules_extension'])
    parser.add_argument('--debug', help="Sets debug mode (much faster / fewer runs for testing)", action='store_true')
    parser.add_argument('--samples', help="Sets the number of MSPRT samples", type=int, default=1)
    parser.add_argument('--beam_type', help = "Whether to calculate beam mass using truncated normal, exponential, or probabilistic",
                        default='exponential', choices=['trunc_norm', 'exponential', 'probabilistic'])
    parser.add_argument('--remove_parameters', help="A list of parameter names to NOT fit",
                        nargs='*', default=None)
    parser.add_argument('--fixed_parameters',
                        help="A list of parameters that can be initialized from -i but won't be fit",
                        nargs='*', default=None)
    parser.add_argument('-r', '--replications', help="Number of times to replicate the cross-validation",
                        type=int, default=None)
    parser.add_argument('--cv_trials', help="Point to a file of cross-validation trials to replicate across runs",
                        default=None, type=argparse.FileType('rU'))
    parser.add_argument('--single_fit', help='Single parameter to fit for rules_extension',
                        default='mass_heur_p')
    parser.add_argument('--init_rules', help="Filename for rules initialization (for rules_extension)",
                        default=None, type=argparse.FileType('rU'))

    args = parser.parse_args()

    print 'Arguments:'
    print args

    if args.type == 'shapes' or args.type == 'complete':
        package = shapes_package
    elif args.type == 'materials':
        package = materials_package
    elif args.type == 'balance':
        package = balance_package
    elif args.type == 'combined':
        package = combined_package
    else:
        raise Exception('Improper package type')

    mod_fn, trs, emp, pnames, rule_fn = package
    pnames = copy.copy(pnames)

    if args.init_parameters is None:
        initpars = None
    else:
        initpars = read_params(args.init_parameters)

    if args.debug:
        nsamps = 10
        maxiter = 5
        igrid = 2
        nreps = 2
        piter = 5
        para = False
    else:
        nsamps = 500
        maxiter = 100
        igrid = 10
        nreps = 10
        piter = 100
        para = True

    if args.replications is not None:
        nreps = args.replications

    if args.mod_class == 'add_ind':
        do_add_ind_cv(mod_fn, trs, emp, pnames, initpars, args.output, args.samples, nreps, nsamps,
                      maxiter, igrid, piter, args.beam_type)
    elif args.mod_class == 'model_type':
        if args.type == 'complete':
            models = {'sh': model_all_shapes, 'mat': model_all_materials, 'bal': model_all_balance}
            trials = {'sh': shapes_trials, 'mat': materials_trials, 'bal': balance_trials}
            empdat = {'sh': shapes_empirical, 'mat': materials_empirical, 'bal': balance_empirical}
            params = {'sh': shapes_param_names, 'mat': materials_param_names, 'bal': balance_param_names}

            if args.beam_type != 'trunc_norm':
                for pnames in params.values():
                    safe_remove(pnames, 'beam_mass_sd')

            if args.beam_type == 'probabilistic':
                for pnames in params.values():
                    pnames.append('beam_use_prob')

            fitpars = list(set([pname for plist in params.values() for pname in plist]))
            if args.remove_parameters is not None:
                for p in args.remove_parameters:
                    safe_remove(fitpars, p)

            do_mod_type_cv_complete(models, trials, empdat, fitpars, params, initpars, args.beam_type,
                                    args.output, args.samples, nreps, nsamps, maxiter, piter, args.cv_trials)
        else:
            do_mod_type_cv(mod_fn, trs, emp, pnames, initpars, args.output, args.samples, nreps, nsamps,
                           maxiter, piter)
    elif args.mod_class == 'comp_rules':
        fitpars = list(set(shapes_param_names + materials_param_names + balance_param_names))
        if args.beam_type != 'trunc_norm':
            safe_remove(fitpars, 'beam_mass_sd')
        if args.beam_type == 'probabilistic':
            fitpars.append('beam_use_prob')
        if args.remove_parameters is not None:
            for p in args.remove_parameters:
                safe_remove(fitpars, p)
        do_comp_rules(fitpars, initpars, args.beam_type, args.output, args.samples, nreps, nsamps, maxiter,
                      piter, args.cv_trials, para, args.remove_parameters)

    elif args.mod_class == 'rules_extension':
        if args.beam_type != 'trunc_norm':
            safe_remove(pnames, 'beam_mass_sd')
        if args.beam_type == 'probabilistic':
            pnames.append('beam_use_prob')
        if args.remove_parameters is not None:
            for p in args.remove_parameters:
                safe_remove(pnames, p)
        if args.init_rules is None:
            initrules = None
        else:
            initrules = read_params(args.init_rules)

        do_rules_extension_cv(mod_fn, rule_fn, trs, emp, args.single_fit, pnames, initpars, initrules,
                              args.output, nreps, nsamps, maxiter, piter, args.beam_type)

    else:
        raise Exception('Not valid cv type')