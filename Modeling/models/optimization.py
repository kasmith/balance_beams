from __future__ import division
from .model_parts import get_all_llh, get_ind_llh, \
    rules_shapes, rules_materials, rules_balance, \
    get_all_trial_llh_from_mix, pack_empirical_data, ALL_MOD_TYPES, \
    get_ind_trial_llh_from_mix
from OptimTools import SPSA, async_map
import copy
import random
import json
import os
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from .config import *
from .writeout import *

"""
__all__ = ['optimize_all', 'optimize_joint_strats_complete',
           'optimize_joint_perc_complete',
           'cross_validate_individual_params',
           'optimize_single_ind_param', 'optimize_complete',
           'optimize_one_joint_complete',
           'optimize_for_individual_rules', 'optimize_mix_for_all',
           'optimize_ind']
"""
__all__ = ['optimize_complete', 'optimize_joint_strats_complete',
           'optimize_joint_perc_complete', 'optimize_joint_strats_byhalf',
           'optimize_ind', 'transform_strat_strengths',
           'revert_strat_strengths']


##########
#
# Optimize models
#
##########

# Turns a vector of reals into probabilities summing to 1
def transform_strat_strengths(strat_strengths, strat_names):
    exp_str = [np.exp(s) for s in strat_strengths]
    tot_str = sum(exp_str)
    tot_str = max(tot_str, 2.) # Keeps guess rate < 1/3rd
    tot_str += 1 # Account for guessing
    r = dict([sn, ss/tot_str] for sn, ss in zip(strat_names, exp_str))
    r['guess'] = 1./tot_str
    return r

def revert_strat_strengths(strat_prob_dict, strat_names):
    for n in strat_names:
        assert n in strat_prob_dict.keys(), (n + " not found in probabilities")
    sn_guess = strat_names + ['guess']
    tot_p = sum([v for k,v in strat_prob_dict.items() if k in sn_guess])
    renorm_probs = []
    for n in sn_guess:
        renorm_probs.append(strat_prob_dict[n] / tot_p)
    sum_str = 1. / renorm_probs[-1]
    return [np.log(p*sum_str) for p in renorm_probs[:-1]]

# As above but deletes / fills strategies to fit strat names
def force_revert_strat_strengths(strat_prob_dict, strat_names):
    spd_names = list(strat_prob_dict.keys())
    spd_names.remove('guess')
    spd_strs = revert_strat_strengths(strat_prob_dict, spd_names)
    str_dict = dict([(spd_names[i], spd_strs[i]) for i in range(len(spd_names))])
    ret = []
    for sn in strat_names:
        ret.append(str_dict.get(sn, 1.))
    return ret

class _strat_holder(object):
    def __init__(self, strat):
        self.strat = strat

def optimize_mix_for_all(tmdict, strategies, empirical_data, init=None):
    if init is None:
        init = np.zeros(len(strategies))
    else:
        assert len(init) == len(strategies), \
            "Initialization & strat lengths don't mix"
    def _get_neg_llh(p):
        if any([i > 10000 for i in p]):
            return 99999999999
        mix_dict = transform_strat_strengths(p, strategies)
        return -get_all_trial_llh_from_mix(tmdict, mix_dict,
                                           empirical_data)
    o = minimize(_get_neg_llh, init)
    return o.x, o.fun

def optimize_all(model_function, model_trials, empirical_data, fit_params_perc,
                 perc_param_names, strategies=ALL_MOD_TYPES,
                 initial_params_perc=None, init_strat_params = None,
                 fixed_params=None, msprt_samples=1, max_iters=500,
                 print_iters=50, n_reps=500, use_full_balance=True,
                 beam_type='exponential', use_dist_max=True):
    if initial_params_perc is None:
        initial_params_perc = get_init_params(fit_params_perc)
    else:
        def_inits = get_init_params(fit_parameters_perc)
        for p in fit_params_perc:
            if p not in initial_params_perc:
                initial_params_perc[p] = def_inits[p]
    if fixed_params is None:
        fill_dict = dict()
        for prm in default_params.keys():
            fill_dict[prm] = default_params[prm][0]
    else:
        fill_dict = dict()
        for p in parameter_names:
            if p in fixed_params:
                fill_dict[p] = initial_params[p]
            else:
                fill_dict[p] = default_params[p]

    if init_strat_params is not None:
        init_strat_params = force_revert_strat_strengths(init_strat_params,
                                                         strategies)

    all_bounds = get_param_bounds(fit_params_perc)
    firstps = [initial_params_perc[p] for p in fit_params_perc]
    bnds = [all_bounds[p] for p in fit_params_perc]

    emp = pack_empirical_data(empirical_data, model_trials.keys())
    sh = _strat_holder(init_strat_params)
    def hfnc(params):
        # Here's the perceptual parameters
        use_param = dict([(fp, prm)
                          for fp, prm in zip(fit_params_perc, params)])
        p = fill_params(use_param, perc_param_names, def_params=fill_dict)
        # Get the model predictions by strategy
        trial_models = make_trial_models(model_function, model_trials, p,
                                         n_times=n_reps)
        # Optimize for the mix
        opt_mix, opt_neg_llh = optimize_mix_for_all(trial_models, strategies,
                                                    emp, sh.strat)
        sh.strat = opt_mix
        return opt_neg_llh

    optparams, llh, _ = SPSA(
        hfnc, firstps, 1e-4, .1, bounds=bnds, maxiter=max_iters,
        print_iters=print_iters)
    optdict_perc = dict([k, p] for k, p in zip(fit_params_perc, optparams))
    trial_models = make_trial_models(model_function, model_trials, optdict_perc)
    opt_mix_raw, opt_neg_llh = optimize_mix_for_all(trial_models, strategies,
                                                    emp, init_strat_params)
    opt_mix = transform_strat_strengths(opt_mix_raw, strategies)
    return optdict_perc, opt_mix, llh


def optimize_complete(model_function_dict, model_trials_dict,
                      empirical_data_dict, fit_params_perc,
                      perc_param_names_dict, strategies=ALL_MOD_TYPES,
                      initial_params_perc=None, init_strat_params = None,
                      fixed_params=None, msprt_samples=1, max_iters=500,
                      print_iters=50, n_reps=500, use_full_balance=True,
                      beam_type='exponential', use_dist_max=True,
                      n_cores=1, savestate=None):
    models = model_function_dict.keys()
    assert (model_trials_dict.keys() == models and
           empirical_data_dict.keys() == models and
           perc_param_names_dict.keys() == models), \
           "Need all models to be identical"

    all_parameter_names = list(
        set([pname for plist in perc_param_names_dict.values()
             for pname in plist]))

    if initial_params_perc is None:
        initial_params_perc = get_init_params(fit_params_perc)
    else:
        def_inits = get_init_params(fit_params_perc)
        for p in fit_params_perc:
            if p not in initial_params_perc:
                initial_params_perc[p] = def_inits[p]
    if fixed_params is None:
        fill_dict = dict()
        for prm in default_params.keys():
            fill_dict[prm] = default_params[prm][0]
    else:
        fill_dict = dict()
        for p in all_parameter_names:
            if p in fixed_params:
                fill_dict[p] = initial_params[p]
            else:
                fill_dict[p] = default_params[p]

    if init_strat_params is not None:
        init_strat_params = force_revert_strat_strengths(init_strat_params,
                                                         strategies)

    all_bounds = get_param_bounds(fit_params_perc)
    firstps = [initial_params_perc[p] for p in fit_params_perc]
    bnds = [all_bounds[p] for p in fit_params_perc]

    flattened_emp = {}
    for mod in models:
        flattened_emp.update(pack_empirical_data(empirical_data_dict[mod],
                                                 model_trials_dict[mod].keys()))

    # First gets the predictions from each model, then combines across models
    def _get_trial_models(perc_params):
        # Set up for multiprocessing
        #  Make stim sets of (name, trial, function, parameters)
        ind_stims = []
        for mod in models:
            parameter_names = perc_param_names_dict[mod]
            model_function = model_function_dict[mod]
            model_trials = model_trials_dict[mod]
            use_param = dict([(fp, prm) for fp, prm in zip(
                fit_params_perc, perc_params) if fp in parameter_names])
            p = fill_params(use_param, parameter_names, def_params=fill_dict)
            ind_stims += [(trnm, tr, model_function, p)
                          for trnm, tr in model_trials.items()]

        def _get_itrial(itp):
            nm, tr, fn, param = itp
            return (nm, fn(tr, strategies, n_times=n_reps, **param))

        if n_cores == 1:
            return dict(map(_get_itrial, ind_stims))
        else:
            return dict(async_map(_get_itrial, ind_stims, ncpu=n_cores))

    # Hold the optimization from last round to initialize new round
    sh = _strat_holder(init_strat_params)
    def hfnc(params):
        tmdict = _get_trial_models(params)
        opt_mix, opt_neg_llh = optimize_mix_for_all(tmdict, strategies,
                                                    flattened_emp,
                                                    sh.strat)
        sh.strat = opt_mix
        return opt_neg_llh

    optparams, llh, _ = SPSA(hfnc, firstps, 1e-4, .1, bounds=bnds,
                             maxiter=max_iters, print_iters=print_iters,
                             savestate=savestate)
    optdict_perc = dict([k, p] for k, p in zip(fit_params_perc, optparams))
    tm_all = _get_trial_models(optparams)
    #tm_flat = {}
    #for tms in tm_all.values():
    #    tm_flat.update(tms)
    opt_mix_raw, opt_neg_llh = optimize_mix_for_all(tm_all, strategies,
                                                    flattened_emp,
                                                    init_strat_params)
    opt_mix = transform_strat_strengths(opt_mix_raw, strategies)
    return optdict_perc, opt_mix, llh



def optimize_joint_perc_complete(model_function_dict, model_trials_dict,
                      empirical_data_dict, fit_params_perc,
                      perc_param_names_dict, strategies=ALL_MOD_TYPES,
                      initial_params_perc=None, init_strat_params = None,
                      fixed_params=None, msprt_samples=1, max_iters=500,
                      print_iters=50, n_reps=500, use_full_balance=True,
                      beam_type='exponential', use_dist_max=True,
                      n_cores=1, savestate=None, enforce_single_strat=False):

    models = model_function_dict.keys()
    assert (model_trials_dict.keys() == models and
           empirical_data_dict.keys() == models and
           perc_param_names_dict.keys() == models), \
           "Need all models to be identical"

    WIDS = dict()
    strat_par_dict = dict()
    for mod in models:
        ws = empirical_data_dict[mod].keys()
        WIDS[mod] = ws
        if init_strat_params is None:
            strat_par_dict.update(dict([(w, np.zeros(len(strategies)))
                                    for w in ws]))
        else:
            strat_par_dict.update(dict([(w,
                                         force_revert_strat_strengths(init_strat_params,
                                                                      strategies))
                                    for w in ws]))


    all_parameter_names = list(
        set([pname for plist in perc_param_names_dict.values()
             for pname in plist]))

    if initial_params_perc is None:
        initial_params_perc = get_init_params(fit_params_perc)
    else:
        def_inits = get_init_params(fit_params_perc)
        for p in fit_params_perc:
            if p not in initial_params_perc:
                initial_params_perc[p] = def_inits[p]
    if fixed_params is None:
        fill_dict = dict()
        for prm in default_params.keys():
            fill_dict[prm] = default_params[prm][0]
    else:
        fill_dict = dict()
        for p in all_parameter_names:
            if p in fixed_params:
                fill_dict[p] = initial_params[p]
            else:
                fill_dict[p] = default_params[p]

    all_bounds = get_param_bounds(fit_params_perc)
    firstps = [initial_params_perc[p] for p in fit_params_perc]
    bnds = [all_bounds[p] for p in fit_params_perc]

    # First gets the predictions from each model, then combines across models
    def _get_trial_models(perc_params):
        # Set up for multiprocessing
        #  Make stim sets of (model, name, trial, function, parameters)
        ind_stims = []
        for mod in models:
            parameter_names = perc_param_names_dict[mod]
            model_function = model_function_dict[mod]
            model_trials = model_trials_dict[mod]
            use_param = dict([(fp, prm) for fp, prm in zip(
                fit_params_perc, perc_params) if fp in parameter_names])
            p = fill_params(use_param, parameter_names, def_params=fill_dict)
            ind_stims += [(mod, trnm, tr, model_function, p)
                          for trnm, tr in model_trials.items()]

        def _get_itrial(itp):
            mod, nm, tr, fn, param = itp
            return (mod, nm, fn(tr, strategies, n_times=n_reps, **param))

        if n_cores == 1:
            tms = map(_get_itrial, ind_stims)
        else:
            tms = async_map(_get_itrial, ind_stims, ncpu=n_cores)

        ret = dict([(m, {}) for m in models])
        for mod, nm, o in tms:
            ret[mod][nm] = o
        return ret

    # Hold the optimization from last round to initialize new round
    def _get_strats_by_perc(params):
        tmdict = _get_trial_models(params)
        # Make stim sets of (wid, model, data)
        io_stims = []
        for m in models:
            for wid, dat in empirical_data_dict[m].items():
                io_stims.append((wid, m, dat))
        def _do_ind_opt(ios):
            wid, mod, dat = ios
            trs = tmdict[mod]
            init_par = strat_par_dict[wid]
            # Get the llh from a flattened strategy set
            def _get_llh(p):
                if any([i > 10000 for i in p]):
                    return 99999999999
                mix_dict = transform_strat_strengths(p, strategies)
                tr_names = list(set(dat.keys()) & set(trs.keys()))
                return -sum([get_ind_trial_llh_from_mix(trs[tr],
                                                        mix_dict,
                                                        dat[tr])
                             for tr in tr_names])
            if enforce_single_strat:
                # For a single strategy we must fit the guess rate for each parameter
                n_pars = len(strategies)
                def _opt_by_param(i):
                    def _hfnc(p):
                        raw_par = [-99999999 for _ in range(n_pars)]
                        raw_par[i] = p
                        return _get_llh(raw_par)
                    return minimize_scalar(_hfnc, bracket=(0, 3))
                all_opts = [_opt_by_param(i) for i in range(n_pars)]
                all_llhs = [o.fun for o in all_opts]
                best_idx = all_llhs.index(min(all_llhs))
                o = all_opts[best_idx]
                remade_param = [-99999999 for _ in range(n_pars)]
                remade_param[best_idx] = o.x
                return (wid, {'model': mod,
                              'strat_params': remade_param,
                              'llh': o.fun})
            else:
                # To allow a mix of strategies we can optimize get_llh
                o = minimize(_get_llh, init_par)
                return (wid, {'model': mod,
                              'strat_params': o.x,
                              'llh': o.fun})

        # Run the model
        if n_cores == 1:
            ind_parse = dict(map(_do_ind_opt, io_stims))
        else:
            ind_parse = dict(async_map(_do_ind_opt, io_stims))

        for wid in strat_par_dict.keys():
            strat_par_dict[wid] = ind_parse[wid]['strat_params']
        return ind_parse

    def hfnc(params):
        iparse = _get_strats_by_perc(params)
        tot_llh = sum([x['llh'] for x in iparse.values()])
        return tot_llh

    optparams, llh, _ = SPSA(hfnc, firstps, 1e-4, .1, bounds=bnds,
                             maxiter=max_iters, print_iters=print_iters,
                             savestate=savestate)
    optdict_perc = dict([k, p] for k, p in zip(fit_params_perc, optparams))
    istrats = _get_strats_by_perc(optparams)
    for wid in istrats:
        istrats[wid]['strat_params'] = \
            transform_strat_strengths(istrats[wid]['strat_params'],
                                      strategies)
    return optdict_perc, istrats, llh


def _joint_strat_recomb_params(key_list, params):
    ret = dict()
    for kl, p in zip(key_list, params):
        wid, paramname, mod = kl
        if mod not in ret.keys():
            ret[mod] = dict()
        if wid not in ret[mod].keys():
            ret[mod][wid] = dict()
        ret[mod][wid][paramname] = p
    return ret

def optimize_joint_strats_complete(model_function_dict, model_trials_dict,
                      empirical_data_dict, fit_params_perc,
                      perc_param_names_dict, strategies=ALL_MOD_TYPES,
                      initial_params_perc=None, init_strat_params = None,
                      fixed_params=None, msprt_samples=1, max_iters=500,
                      print_iters=50, n_reps=500, use_full_balance=True,
                      beam_type='exponential', use_dist_max=True,
                      n_cores=1, savestate=None):

    models = model_function_dict.keys()
    assert (model_trials_dict.keys() == models and
           empirical_data_dict.keys() == models and
           perc_param_names_dict.keys() == models), \
           "Need all models to be identical"

    all_parameter_names = list(
        set([pname for plist in perc_param_names_dict.values()
             for pname in plist]))

    if initial_params_perc is None:
        initial_params_perc = get_init_params(fit_params_perc)
    else:
        def_inits = get_init_params(fit_params_perc)
        for p in fit_params_perc:
            if p not in initial_params_perc:
                initial_params_perc[p] = def_inits[p]
    if fixed_params is None:
        fill_dict = dict()
        for prm in default_params.keys():
            fill_dict[prm] = default_params[prm][0]
    else:
        fill_dict = dict()
        for p in all_parameter_names:
            if p in fixed_params:
                fill_dict[p] = initial_params[p]
            else:
                fill_dict[p] = default_params[p]

    if init_strat_params is not None:
        init_strat_params = force_revert_strat_strengths(init_strat_params,
                                                         strategies)
    sh = _strat_holder(init_strat_params)

    # Make WIDs and associated parameterization
    WIDS = dict()
    perc_par_init = []
    perc_par_keys = []
    flat_emp_dat = dict()
    for mod in models:
        ws = empirical_data_dict[mod].keys()
        parameter_names = perc_param_names_dict[mod]
        WIDS[mod] = ws
        for wid in ws:
            flat_emp_dat[wid] = empirical_data_dict[mod][wid]
            for pnm in parameter_names:
                if pnm in initial_params_perc.keys():
                    perc_par_init.append(initial_params_perc[pnm])
                    perc_par_keys.append((wid, pnm, mod))

    all_bounds = get_param_bounds(fit_params_perc)
    bnds = [all_bounds[p[1]] for p in perc_par_keys]

    # First gets the predictions from each model, then combines across models
    def _get_trial_models(dat):
        perc_params, mod, _ = dat
        # Set up for multiprocessing
        #  Make stim sets of (name, trial, function, parameters)
        parameter_names = perc_param_names_dict[mod]
        model_function = model_function_dict[mod]
        model_trials = model_trials_dict[mod]
        use_param = dict([(fp, perc_params[fp]) for fp in fit_params_perc
                          if fp in parameter_names])
        p = fill_params(use_param, parameter_names, def_params=fill_dict)

        def _get_itrial(tr):
            return model_function(tr, strategies, n_times=n_reps, **p)
        return dict([(trnm, _get_itrial(tr))
                     for trnm, tr in model_trials.items()])

    def _opt_over_strat(params):
        pardict = _joint_strat_recomb_params(perc_par_keys, params)
        trmod_input = []
        for mod in models:
            pdat = pardict[mod]
            for wid, dat in pdat.items():
                trmod_input.append((dat, mod, wid))
        if n_cores == 1:
            trmod_out = map(_get_trial_models, trmod_input)
        else:
            trmod_out = async_map(_get_trial_models, trmod_input, ncpu=n_cores)
        # Store the outcomes for each strategy
        ocm_strat_dict = dict()
        for tin, tout in zip(trmod_input, trmod_out):
            dat, mod, wid = tin
            ocm_strat_dict[wid] = tout
        # Optimize the strategy Parameters
        flat_wids = []
        for mod in models:
            flat_wids += WIDS[mod]
        def _neg_llh_by_strat(stratpar):
            if any([i > 10000 for i in stratpar]):
                return 99999999999
            mix_dict = transform_strat_strengths(stratpar, strategies)
            llh = 0
            for w in flat_wids:
                smods = ocm_strat_dict[w]
                emp = flat_emp_dat[w]
                utrs = list(set(emp.keys()) & set(smods.keys()))
                for tr in utrs:
                    llh -= get_ind_trial_llh_from_mix(smods[tr], mix_dict, emp[tr])
            return llh
        init_sp = sh.strat
        if init_sp is None:
            init_sp = np.zeros(len(strategies))
        o = minimize(_neg_llh_by_strat, init_sp)
        sh.strat = o.x
        return o.fun, o.x

    def hfnc(params):
        return _opt_over_strat(params)[0]

    optparams, llh, _ = SPSA(hfnc, perc_par_init, 1e-4, .1, bounds=bnds,
                             maxiter=max_iters, print_iters=print_iters,
                             savestate=savestate)
    opt_perc_params = _joint_strat_recomb_params(perc_par_keys, optparams)
    sllh, opt_mix_raw = _opt_over_strat(optparams)
    opt_mix = transform_strat_strengths(opt_mix_raw, strategies)
    return opt_perc_params, opt_mix, llh


def optimize_ind(model_function_dict, model_trials_dict,
                      empirical_data_dict, fit_params_perc,
                      perc_param_names_dict, strategies=ALL_MOD_TYPES,
                      initial_params_perc=None, init_strat_params = None,
                      fixed_params=None, msprt_samples=1, max_iters=500,
                      print_iters=50, n_reps=500, use_full_balance=True,
                      beam_type='exponential', use_dist_max=True,
                      n_cores=1, savestate=None, enforce_single_strat=False):

    models = model_function_dict.keys()
    assert (model_trials_dict.keys() == models and
           empirical_data_dict.keys() == models and
           perc_param_names_dict.keys() == models), \
           "Need all models to be identical"

    WIDS = dict()
    strat_par_dict = dict()
    mod_by_wid = dict()
    flat_wids = []
    for mod in models:
        ws = empirical_data_dict[mod].keys()
        WIDS[mod] = ws
        for w in ws:
            mod_by_wid[w] = mod
        flat_wids += ws

    if init_strat_params is not None:
        init_strat_params = force_revert_strat_strengths(init_strat_params,
                                                         strategies)
    else:
        init_strat_params = np.zeros(len(strategies))

    all_parameter_names = list(
        set([pname for plist in perc_param_names_dict.values()
             for pname in plist]))

    if initial_params_perc is None:
        initial_params_perc = get_init_params(fit_params_perc)
    else:
        def_inits = get_init_params(fit_params_perc)
        for p in fit_params_perc:
            if p not in initial_params_perc:
                initial_params_perc[p] = def_inits[p]
    if fixed_params is None:
        fill_dict = dict()
        for prm in default_params.keys():
            fill_dict[prm] = default_params[prm][0]
    else:
        fill_dict = dict()
        for p in all_parameter_names:
            if p in fixed_params:
                fill_dict[p] = initial_params[p]
            else:
                fill_dict[p] = default_params[p]

    def get_opt_from_wid(wid):
        mod = mod_by_wid[wid]
        parameter_names = perc_param_names_dict[mod]
        model_function = model_function_dict[mod]
        model_trials = model_trials_dict[mod]
        empirical = empirical_data_dict[mod]

        these_fit_params = [p for p in fit_params_perc if p in parameter_names]

        all_bounds = get_param_bounds(parameter_names)
        firstps = [initial_params_perc[p] for p in these_fit_params]
        bnds = [all_bounds[p] for p in these_fit_params]
        emp = empirical[wid]
        tr_use = list(set(emp.keys()) & set(model_trials.keys()))
        sh = _strat_holder(copy.copy(init_strat_params))


        def _fit_from_strat(perc_params):
            ind_stims = []
            use_param = dict([(fp, prm) for fp, prm in zip(
                these_fit_params, perc_params) if fp in parameter_names])
            p = fill_params(use_param, parameter_names, def_params=fill_dict)
            ind_stims += [(trnm, model_trials[trnm], model_function, p)
                          for trnm in tr_use]

            def _get_itrial(itp):
                nm, tr, fn, param = itp
                return (nm, fn(tr, strategies, n_times=n_reps, **param))

            tm_dict = dict(map(_get_itrial, ind_stims))

            def _neg_llh_from_strat(stratpar):
                if any([i > 10000 for i in stratpar]):
                    #return 99999999999, [-99999999 for _ in range(len(strategies))]
                    return 99999999999
                mix_dict = transform_strat_strengths(stratpar, strategies)
                ret = [get_ind_trial_llh_from_mix(tm_dict[tr],
                                                        mix_dict,
                                                        emp[tr])
                             for tr in tr_use]
                return -sum(ret)
            if enforce_single_strat:
                # Fix participants to a single strategy
                n_pars = len(strategies)
                def _opt_by_ind(i):
                    def _hfnc(p):
                        raw_par = [-99999999 for _ in range(n_pars)]
                        raw_par[i] = p
                        return _neg_llh_from_strat(raw_par)
                    return minimize_scalar(_hfnc, bracket=(0, 3))
                all_opts = [_opt_by_ind(i) for i in range(n_pars)]
                all_llhs = [o.fun for o in all_opts]
                best_idx = all_llhs.index(min(all_llhs))
                o = all_opts[best_idx]
                remade_param = [-99999999 for _ in range(n_pars)]
                remade_param[best_idx] = o.x
                sh.strat = remade_param
                return o.fun, remade_param
            else:
                # Allow inidividual mixtures of strategies
                o = minimize(_neg_llh_from_strat, sh.strat)
                sh.strat = o.x
                return o.fun, o.x

        def hfnc(perc_params):
            return _fit_from_strat(perc_params)[0]

        optparams, llh, _ = SPSA(hfnc, firstps, 1e-4, .1, bounds=bnds,
                                 maxiter=max_iters, print_iters=print_iters,
                                 savestate=savestate[:-5]+"_"+wid+".json")
        optdict_perc = dict([k, p] for k, p in zip(these_fit_params, optparams))
        ollh, optstrat = _fit_from_strat(optparams)
        opt_mix = transform_strat_strengths(optstrat, strategies)

        return optdict_perc, opt_mix, llh

    if n_cores == 1:
        indfits = map(get_opt_from_wid, flat_wids)
    else:
        indfits = async_map(get_opt_from_wid, flat_wids, n_cores)

    perc_dict = dict([(m, dict()) for m in models])
    strat_dict = dict()
    tot_llh = 0
    for w, ifit in zip(flat_wids, indfits):
        mod = mod_by_wid[w]
        perc_dict[mod][w] = ifit[0]
        strat_dict[w] = ifit[1]
        tot_llh += ifit[2]

    return perc_dict, strat_dict, tot_llh

def optimize_joint_strats_byhalf(model_function_dict, model_trials_dict,
                      empirical_data_dict, fit_params_perc,
                      perc_param_names_dict, orders_dict,
                      strategies=ALL_MOD_TYPES,
                      initial_params_perc=None, init_strat_params = None,
                      fixed_params=None, msprt_samples=1, max_iters=500,
                      print_iters=50, n_reps=500, use_full_balance=True,
                      beam_type='exponential', use_dist_max=True,
                      n_cores=1, savestate=None):

    models = model_function_dict.keys()
    assert (model_trials_dict.keys() == models and
           empirical_data_dict.keys() == models and
           perc_param_names_dict.keys() == models and
           orders_dict.keys() == models), \
           "Need all models to be identical"

    WIDS = dict()
    strat_par_dict_0 = dict()
    emp_splits = dict()
    for mod in models:
        ws = empirical_data_dict[mod].keys()
        WIDS[mod] = ws
        # Split trials into first-half / second-half
        emp_ords = orders_dict[mod]
        trial_names = model_trials_dict[mod].keys()
        tr_cutoff = int(len(trial_names)/2)
        this_emp_split = {}
        for w in ws:
            ord_d = emp_ords[w]
            this_emp_split[w] = [[], []]
            for trnm, o in ord_d.items():
                if o < tr_cutoff:
                    this_emp_split[w][0].append(trnm)
                else:
                    this_emp_split[w][1].append(trnm)
        emp_splits[mod] = this_emp_split
        if init_strat_params is None:
            strat_par_dict_0.update(dict([(w, np.zeros(len(strategies)))
                                    for w in ws]))
        else:
            strat_par_dict_0.update(dict([(w,
                                         force_revert_strat_strengths(init_strat_params,
                                                                      strategies))
                                    for w in ws]))
    # Initialize to the same values
    strat_par_dict_1 = copy.deepcopy(strat_par_dict_0)

    all_parameter_names = list(
        set([pname for plist in perc_param_names_dict.values()
             for pname in plist]))

    if initial_params_perc is None:
        initial_params_perc = get_init_params(fit_params_perc)
    else:
        def_inits = get_init_params(fit_params_perc)
        for p in fit_params_perc:
            if p not in initial_params_perc:
                initial_params_perc[p] = def_inits[p]
    if fixed_params is None:
        fill_dict = dict()
        for prm in default_params.keys():
            fill_dict[prm] = default_params[prm][0]
    else:
        fill_dict = dict()
        for p in all_parameter_names:
            if p in fixed_params:
                fill_dict[p] = initial_params[p]
            else:
                fill_dict[p] = default_params[p]

    all_bounds = get_param_bounds(fit_params_perc)
    firstps = [initial_params_perc[p] for p in fit_params_perc]
    bnds = [all_bounds[p] for p in fit_params_perc]

    # First gets the predictions from each model, then combines across models
    def _get_trial_models(perc_params):
        # Set up for multiprocessing
        #  Make stim sets of (model, name, trial, function, parameters)
        ind_stims = []
        for mod in models:
            parameter_names = perc_param_names_dict[mod]
            model_function = model_function_dict[mod]
            model_trials = model_trials_dict[mod]
            use_param = dict([(fp, prm) for fp, prm in zip(
                fit_params_perc, perc_params) if fp in parameter_names])
            p = fill_params(use_param, parameter_names, def_params=fill_dict)
            ind_stims += [(mod, trnm, tr, model_function, p)
                          for trnm, tr in model_trials.items()]

        def _get_itrial(itp):
            mod, nm, tr, fn, param = itp
            return (mod, nm, fn(tr, strategies, n_times=n_reps, **param))

        if n_cores == 1:
            tms = map(_get_itrial, ind_stims)
        else:
            tms = async_map(_get_itrial, ind_stims, ncpu=n_cores)

        ret = dict([(m, {}) for m in models])
        for mod, nm, o in tms:
            ret[mod][nm] = o
        return ret

    # Hold the optimization from last round to initialize new round
    def _get_strats_by_perc(params):
        tmdict = _get_trial_models(params)
        # Make stim sets of (wid, model, data)
        io_stims_0 = []
        io_stims_1 = []
        for m in models:
            for wid, dat in empirical_data_dict[m].items():
                ords = emp_splits[m][wid]
                trs_0 = dict([(tr, dat[tr]) for tr in ords[0]])
                trs_1 = dict([(tr, dat[tr]) for tr in ords[1]])
                io_stims_0.append((wid, m, trs_0, 0))
                io_stims_1.append((wid, m, trs_1, 1))
        def _do_ind_opt(ios):
            wid, mod, dat, half = ios
            trs = tmdict[mod]
            init_par = [strat_par_dict_0[wid], strat_par_dict_1[wid]][half]
            # Get the llh from a flattened strategy set
            def _get_llh(p):
                if any([i > 10000 for i in p]):
                    return 99999999999
                mix_dict = transform_strat_strengths(p, strategies)
                tr_names = list(set(dat.keys()) & set(trs.keys()))
                return -sum([get_ind_trial_llh_from_mix(trs[tr],
                                                        mix_dict,
                                                        dat[tr])
                             for tr in tr_names])

            # To allow a mix of strategies we can optimize get_llh
            o = minimize(_get_llh, init_par)
            return (wid, {'model': mod,
                          'strat_params': o.x,
                          'llh': o.fun})

        # Run the model
        if n_cores == 1:
            ind_parse_0 = dict(map(_do_ind_opt, io_stims_0))
            ind_parse_1 = dict(map(_do_ind_opt, io_stims_1))
        else:
            ind_parse_0 = dict(async_map(_do_ind_opt, io_stims_0))
            ind_parse_1 = dict(async_map(_do_ind_opt, io_stims_1))

        ret_dict = dict()
        for wid in strat_par_dict_0.keys():
            strat_par_dict_0[wid] = ind_parse_0[wid]['strat_params']
            strat_par_dict_1[wid] = ind_parse_1[wid]['strat_params']
            ret_dict[wid] = {
                'model': ind_parse_0[wid]['model'],
                'strat_params_first': ind_parse_0[wid]['strat_params'],
                'strat_params_second': ind_parse_1[wid]['strat_params'],
                'llh': ind_parse_0[wid]['llh'] + ind_parse_1[wid]['llh']
            }
        return ret_dict

    def hfnc(params):
        iparse = _get_strats_by_perc(params)
        tot_llh = sum([x['llh'] for x in iparse.values()])
        return tot_llh

    optparams, llh, _ = SPSA(hfnc, firstps, 1e-4, .1, bounds=bnds,
                             maxiter=max_iters, print_iters=print_iters,
                             savestate=savestate)
    optdict_perc = dict([k, p] for k, p in zip(fit_params_perc, optparams))
    istrats = _get_strats_by_perc(optparams)
    for wid in istrats:
        istrats[wid]['strat_params_first'] = \
            transform_strat_strengths(istrats[wid]['strat_params_first'],
                                      strategies)
        istrats[wid]['strat_params_second'] = \
            transform_strat_strengths(istrats[wid]['strat_params_second'],
                                      strategies)
    return optdict_perc, istrats, llh
