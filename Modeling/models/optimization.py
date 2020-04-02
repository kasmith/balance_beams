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
                fit_params_perc, perc_params) if fp in parameter_names])
            p = fill_params(use_param, parameter_names, def_params=fill_dict)
            ind_stims += [(trnm, model_trials[trnm], model_function, p)
                          for trnm in tr_use]

            def _get_itrial(itp):
                nm, tr, fn, param = itp
                return (nm, fn(tr, strategies, n_times=n_reps, **param))

            tm_dict = dict(map(_get_itrial, ind_stims))

            def _neg_llh_from_strat(stratpar):
                if any([i > 10000 for i in stratpar]):
                    return 99999999999, [-99999999 for _ in range(len(strategies))]
                mix_dict = transform_strat_strengths(stratpar, strategies)
                return -sum([get_ind_trial_llh_from_mix(tm_dict[tr],
                                                        mix_dict,
                                                        emp[tr])
                             for tr in tr_use])
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










""" OLD CODE: DELETE SOON
def optimize_joint(model_function, model_trials, empirical_data, fit_aggregate, fit_individual, parameter_names,
                   initial_params=None, fixed_params=None, msprt_samples=1, max_iters=500, print_iters=50,
                   repeats=5, n_reps=500, use_full_balance=True, modtype='smp',
                   beam_type='trunc_norm', use_dist_max=True):
    WIDS = empirical_data.keys()
    # Load in initialization parameters
    if initial_params is None:
        agg_inits = get_init_params(fit_aggregate)
        ind_inits = dict([(w, get_init_params(fit_individual)) for w in WIDS])
    else:
        def_inits = get_init_params(fit_aggregate + fit_individual)
        # Check if we're copying initial parameters from single fitting
        if initial_params.keys()[0] in parameter_names:
            agg_inits = dict()
            ind_inits = dict()
            for p in fit_aggregate:
                if p in initial_params.keys():
                    agg_inits[p] = initial_params[p]
                else:
                    agg_inits[p] = def_inits[p]
            for w in WIDS:
                ind_inits[w] = dict()
                for p in fit_individual:
                    if p in initial_params.keys():
                        ind_inits[w][p] = initial_params[p]
                    else:
                        ind_inits[w][p] = def_inits[p]
        else:
            agg_store = dict([(p, []) for p in fit_aggregate])
            ind_inits = dict()
            for w in WIDS:
                ind_inits[w] = dict()
                for p in fit_aggregate:
                    if p in initial_params[w].keys():
                        agg_store[p].append(initial_params[w][p])
                for p in fit_individual:
                    if p in initial_params.keys():
                        ind_inits[w][p] = initial_params[p]
                    else:
                        ind_inits[w][p] = def_inits[p]
            agg_inits = dict()
            for p in fit_aggregate:
                if len(agg_store[p]) > 0:
                    agg_inits[p] = np.mean(agg_store[p])
                else:
                    agg_inits[p] = def_inits[p]

    if fixed_params is None:
        fill_dict = dict()
        for w in WIDS:
            fill_dict[w] = dict()
            for prm in default_params.keys():
                fill_dict[prm] = default_params[prm][0]
    else:
        fill_dict = dict()
        for w in WIDS:
            fill_dict[w] = dict()
            for p in parameter_names:
                if p in fixed_params:
                    if p in agg_inits:
                        fill_dict[w][p] = agg_inits[p]
                    elif p in ind_inits[w]:
                        fill_dict[w][p] = ind_inits[w][p]
                    else:
                        if initial_params.keys()[0] in parameter_names:
                            fill_dict[w][p] = initial_params[p]
                        else:
                            fill_dict[w][p] = initial_params[w][p]
                else:
                    fill_dict[w][p] = default_params[p]

    # Helper function for allowing an easy split of the likelihood function
    # into aggregate / individual
    def split_llh(params, params_to_fit, allps, indps, ind_emp, fill):
        ps = copy.deepcopy(allps)
        ps.update(copy.deepcopy(indps))
        ps.update(dict([(pnm, p) for pnm, p in zip(params_to_fit, params)]))
        ps = fill_params(ps, parameter_names, fill)
        return -get_ind_llh(model_function, model_trials, ind_emp, fit_all_trials=False, n_times=n_reps,
                            use_full_balance=use_full_balance, parallelize=False, msprt_samples=msprt_samples,
                            type=modtype, beam_type=beam_type, use_dist_max=use_dist_max, **ps)

    all_bound_raw = get_param_bounds(fit_aggregate)
    all_bnds = [all_bound_raw[p] for p in fit_aggregate]
    ind_bound_raw = get_param_bounds(fit_individual)
    ind_bnds = [ind_bound_raw[p] for p in fit_individual]

    for i in range(repeats):
        print '\nIteration', i, 'starting'
        print 'Initial agg_fits:', agg_inits, '\n'

        def help_agg(p):
            def _agg_helper(w):
                return split_llh(p, fit_aggregate, agg_inits, ind_inits[w], empirical_data[w], fill_dict[w])
            return sum(async_map(_agg_helper, WIDS))

        flat_agg_params = [agg_inits[p] for p in fit_aggregate]
        new_agg, midllh, _ = SPSA(help_agg, flat_agg_params, 1e-4, .1, bounds=all_bnds, maxiter=max_iters,
                                  print_iters=print_iters)
        agg_inits = dict([(pnm, p) for pnm, p in zip(fit_aggregate, new_agg)])
        print 'Iteration', i, 'aggregate fitting done; LLH:', midllh

        def help_ind(w):
            def _ind_helper(p):
                return split_llh(p, fit_individual, agg_inits, ind_inits[w], empirical_data[w], fill_dict[w])
            flat_ind_params = [ind_inits[w][p] for p in fit_individual]
            ps, llh, _ = SPSA(_ind_helper, flat_ind_params, 1e-4, .1, bounds=ind_bnds, maxiter=max_iters,
                              print_iters=print_iters)
            indp_dict = dict([(pnm, p) for pnm, p in zip(fit_individual, ps)])
            print 'Iteration', i, 'worker', w, 'done; LLH:', llh
            return indp_dict

        ind_outcomes = async_map(help_ind, WIDS)
        for w, nps in zip(WIDS, ind_outcomes):
            ind_inits[w] = nps

        print '\nIteration', i, 'complete\n'

    llhdict = dict([(w, -split_llh(dict(), [], agg_inits, ind_inits[w],
                                   empirical_data[w], fill_dict[w])) for w in WIDS])

    # Update individual parameters to attach aggregate & llh
    ret = dict()
    for w in WIDS:
        params = copy.deepcopy(ind_inits[w])
        params.update(copy.deepcopy(agg_inits))
        if fixed_params is not None:
            for p in fixed_params:
                if p in agg_inits:
                    params[p] = agg_inits[p]
                elif p in ind_inits[w]:
                    params[p] = ind_inits[w][p]
                else:
                    if initial_params.keys()[0] in parameter_names:
                        params[p] = initial_params[p]
                    else:
                        params[p] = initial_params[w][p]
        ret[w] = (params, llhdict[w])

    return ret


def optimize_one_joint_complete(model_function_dict, model_trials_dict, empirical_data_dict, fit_individual, fit_agg,
                                parameter_names_dict, initial_params=None, fixed_params=None, msprt_samples=1,
                                max_iters=500, print_iters=50, repeats=5, n_reps=500, use_full_balance=True,
                                modtype='smp', beam_type='trunc_norm', use_dist_max=True, initialization_range=10,
                                parallelize=True):

    models = model_function_dict.keys()
    assert model_trials_dict.keys() == models and empirical_data_dict.keys() == models and \
        parameter_names_dict.keys() == models, "Need all models to be identical"

    if parallelize:
        mapfn = async_map
    else:
        mapfn = map

    WIDS = list(set([w for e in empirical_data_dict.values()
                     for w in e.keys()]))
    mod_by_wid = dict([(w, m)
                       for m in models for w in empirical_data_dict[m].keys()])
    all_parameter_names = list(
        set([pname for plist in parameter_names_dict.values() for pname in plist]))

    agg_parameter_names = [
        pnm for pnm in all_parameter_names if pnm in fit_agg]
    # agg_parameter_names.remove(fit_individual)

    # Load in initialization parameters
    if initial_params is None:
        agg_inits = get_init_params(agg_parameter_names)
        ind_inits = dict(
            [(w, get_init_params([fit_individual])[fit_individual]) for w in WIDS])
    else:
        def_inits = get_init_params(agg_parameter_names + [fit_individual])
        # Assume we're copying initial parameters from single fittin
        if initial_params.keys()[0] not in all_parameter_names:
            raise Exception('Cannot load from individual fits')

        agg_inits = dict()
        for p in agg_parameter_names:
            if p in initial_params.keys():
                agg_inits[p] = initial_params[p]
            else:
                agg_inits[p] = def_inits[p]

        if fit_individual in initial_params.keys():
            ind_st = initial_params[fit_individual]
        else:
            ind_st = def_inits[fit_individual]
        ind_inits = dict([(w, ind_st) for w in WIDS])

    if fixed_params is None:
        fill_dict = dict()
        for prm in default_params.keys():
            fill_dict[prm] = default_params[prm][0]
    else:
        fill_dict = dict()
        for p in all_parameter_names:
            if p in fixed_params:
                if p == fit_individual:
                    fill_dict[p] = ind_inits.values()[0]
                elif p in agg_inits.keys():
                    fill_dict[p] = agg_inits[p]
                elif p in initial_params.keys():
                    fill_dict[p] = initial_params[p]
                else:
                    fill_dict[p] = default_params[p][0]
            else:
                fill_dict[p] = default_params[p][0]

    all_bound_raw = get_param_bounds(agg_parameter_names)
    all_bnds = [all_bound_raw[p] for p in agg_parameter_names]
    ind_bound_raw = get_param_bounds([fit_individual])
    ind_bnds = ind_bound_raw[fit_individual]

    # Helper function for allowing an easy split of the likelihood function
    # into aggregate / individual
    def agg_llh(params, ind_ps):
        pdict = dict([(pnm, p) for pnm, p in zip(agg_parameter_names, params)])
        pdict = fill_params(pdict, all_parameter_names, fill_dict)
        llh = 0
        for m in models:
            mod_fn = model_function_dict[m]
            emp_dat = empirical_data_dict[m]
            trs = model_trials_dict[m]
            pnames = parameter_names_dict[m]

            these_params = dict([(pnm, pdict[pnm]) for pnm in pnames])

            these_wids = emp_dat.keys()

            def ifit(w):
                ind_pdict = copy.copy(these_params)
                ind_pdict[fit_individual] = ind_ps[w]
                ind_emp = emp_dat[w]
                return -get_ind_llh(mod_fn, trs, ind_emp, fit_all_trials=False, n_times=n_reps,
                                    use_full_balance=use_full_balance, parallelize=False, msprt_samples=msprt_samples,
                                    type=modtype, beam_type=beam_type, use_dist_max=use_dist_max, **ind_pdict)
            llhs = mapfn(ifit, these_wids)
            llh += sum(llhs)
        return llh

    def ind_llh(i_param, w, agg_ps):
        mod = mod_by_wid[w]
        mod_fn = model_function_dict[mod]
        ind_emp = empirical_data_dict[mod][w]
        trs = model_trials_dict[mod]
        pnames = parameter_names_dict[mod]
        pdict = fill_params(agg_ps, all_parameter_names, fill_dict)

        these_params = dict([(pnm, pdict[pnm])
                             for pnm in pnames if pnm != fit_individual])
        these_params[fit_individual] = i_param
        return -get_ind_llh(mod_fn, trs, ind_emp, fit_all_trials=False, n_times=n_reps,
                            use_full_balance=use_full_balance, parallelize=False, msprt_samples=msprt_samples,
                            type=modtype, beam_type=beam_type, use_dist_max=use_dist_max, **these_params)

    def grid_init(w, ind_inits, agg_ps):

        stp = (ind_bnds[1] - ind_bnds[0]) / (initialization_range + 1)
        min_diff = (ind_bnds[1] - ind_bnds[0]) / 100
        starts = [ind_bnds[0] + min_diff, ind_bnds[1] - min_diff] + [ind_bnds[0] +
                                                                    (i + 1) * stp for i in range(initialization_range)] + [ind_inits[w]]
        starts_val = map(lambda p: ind_llh(p, w, agg_ps), starts)
        sidx = [i for i, j in enumerate(starts_val) if j == min(starts_val)][0]
        return starts[sidx]

    for i in range(repeats):
        print '\nIteration', i, 'starting'
        print 'Initial agg_fits:', agg_inits, '\n'

        flat_agg_params = [agg_inits[p] for p in agg_parameter_names]
        new_agg, midllh, _ = SPSA(lambda p: agg_llh(p, ind_inits), flat_agg_params, 1e-4, .1, bounds=all_bnds, maxiter=max_iters,
                                  print_iters=print_iters)
        agg_inits = dict([(pnm, p)
                          for pnm, p in zip(agg_parameter_names, new_agg)])
        print 'Iteration', i, 'aggregate fitting done; LLH:', midllh

        # On the first iteration, find a good place to initialize individual
        # parameters
        if i == 0:
            ind_starts = async_map(lambda w: grid_init(
                w, ind_inits, agg_inits), WIDS)
            ind_inits = dict([(w, s) for w, s in zip(WIDS, ind_starts)])
            print 'Initialized individuals:', ind_inits, '\n'

        def help_ind(w):
            ps, llh, _ = SPSA(lambda p: ind_llh(p[0], w, agg_inits), [ind_inits[w]], 1e-4, .1, bounds=[ind_bnds],
                              maxiter=max_iters, print_iters=print_iters)
            print 'Iteration', i, 'worker', w, 'done; LLH:', llh
            return ps[0]

        ind_outcomes = mapfn(help_ind, WIDS)
        for w, nps in zip(WIDS, ind_outcomes):
            ind_inits[w] = nps

        print '\nIteration', i, 'complete\n'

    llhdict = dict([(w, -ind_llh(ind_inits[w], w, agg_inits)) for w in WIDS])

    # Update individual parameters to attach aggregate & llh
    ret = dict()
    for w in WIDS:
        params = copy.deepcopy(agg_inits)
        params[fit_individual] = ind_inits[w]
        if fixed_params is not None:
            for p in fixed_params:
                if p in agg_inits:
                    params[p] = agg_inits[p]
                elif p == fit_individual:
                    params[p] = ind_inits[w]
                elif p in initial_params:
                    params[p] = initial_params[p]
        ret[w] = (params, llhdict[w])

    return ret


def optimize_single_ind_param(model_function, model_trials, empirical_data, fit_param, parameter_names,
                              initial_params, msprt_samples=1, max_iters=100, i_range=10, print_iters=50,
                              n_reps=500, use_full_balance=True, modtype='smp',
                              beam_type='trunc_norm', use_dist_max=True, parallelize=True):
    WIDS = empirical_data.keys()
    # Check if we're copying initial parameters from single fitting
    assert initial_params.keys()[
        0] in parameter_names, "Must copy from aggregate fitting"
    aggpars = copy.deepcopy(initial_params)

    if fit_param in initial_params.keys():
        init = aggpars[fit_param]
        del aggpars[fit_param]

    else:
        init = get_init_params([fit_param])[fit_param]

    bnds = get_param_bounds([fit_param])[fit_param]

    def _fit_single(wid):
        print 'starting', wid
        empdat = empirical_data[wid]

        def hfnc(par):
            ip = dict([(fit_param, par)])
            p = fill_params(ip, parameter_names, aggpars)
            return -get_ind_llh(model_function, model_trials, empdat, fit_all_trials=False,
                                n_times=n_reps, use_full_balance=use_full_balance, msprt_samples=msprt_samples,
                                parallelize=False, type=modtype, beam_type=beam_type, use_dist_max=use_dist_max,
                                **p)

        # Calculate the value at various points within the bounds
        stp = (bnds[1] - bnds[0]) / (i_range + 1)
        starts = [bnds[0] + (i + 1) * stp for i in range(i_range)] + [init]
        starts_val = map(hfnc, starts)

        sidx = [i for i, j in enumerate(starts_val) if j == min(starts_val)][0]

        optparams, llh, _ = SPSA(hfnc, [
                                 starts[sidx]], 1e-4, .1, bounds=[bnds], maxiter=max_iters, print_iters=print_iters)
        return optparams, llh

    if parallelize:
        fits = async_map(_fit_single, WIDS)
    else:
        fits = map(_fit_single, WIDS)
    paramdict = dict()
    llhdict = dict()
    for w, (ps, llh) in zip(WIDS, fits):
        paramdict[w] = ps[0]
        llhdict[w] = llh
    return paramdict, llhdict


def _get_ind_cv_llhs(model_function, cv_trials, empirical_data, param_dict, param_names, msprt_samples=1,
                     n_times=500, use_full_balance=True, parallelize=True, type='smp', use_mass_ratio=False):
    WIDS = empirical_data.keys()
    llhdict = dict()
    for w in WIDS:
        emp = empirical_data[w]
        params = fill_params(param_dict[w], param_names)
        llhdict[w] = get_ind_llh(model_function, cv_trials, emp, n_times=n_times, msprt_samples=msprt_samples,
                                 use_full_balance=use_full_balance, parallelize=parallelize, type=type,
                                 use_mass_ratio=use_mass_ratio, **params)
    return llhdict


def cross_validate_individual_params(model_function, model_trials, empirical_data, cv_parameters, parameter_names,
                                     output_dir='.', llh_delta=5, initial_params=None, fixed_params=None,
                                     msprt_samples=1, max_iters=500, repeats=5, n_reps=500,
                                     use_full_balance=True, modtype='full', use_mass_ratio=False, do_drop=True):

    # Create the output directory as needed
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Pick out the trials to cross-validate
    n_trials = len(model_trials)
    trs = copy.copy(model_trials.keys())
    random.shuffle(trs)
    fittrs = copy.copy(trs[:int(n_trials / 2)])
    cvtrs = copy.copy(trs[int(n_trials / 2):])
    fit_trials = dict([(t, model_trials[t]) for t in fittrs])
    cv_trials = dict([(t, model_trials[t]) for t in cvtrs])
    del trs
    cvfl = open(os.path.join(output_dir, 'CV_Trials.json'), 'w')
    json.dump({'FitTrials': fittrs, 'CVTrials': cvtrs}, cvfl)
    cvfl.close()

    # Set up the file that summarizes the llh traces
    summary_file_name = os.path.join(output_dir, 'CV_Summary.csv')
    trace_file_name = os.path.join(output_dir, 'CV_Trace.txt')
    sfl = open(summary_file_name, 'w')
    sfl.write('FileName,')
    for p in cv_parameters:
        sfl.write(p + ',')
    sfl.write('Fit_LLH,CV_LLH\n')
    sfl.close()
    tfl = open(trace_file_name, 'w')
    tfl.write('Starting cross-validation -- ')
    if do_drop:
        tfl.write('dropping from all')
    else:
        tfl.write('starting from none')
    tfl.write('\nUsing the following parameters:\n\t' +
              str(cv_parameters) + '\n\n')
    tfl.write('Allowable LLH difference: ' + str(llh_delta) + '\n\n')
    tfl.close()

    print "Beginning cross-validation with the following parameters:"
    print cv_parameters
    print ""

    # Fit the baseline model (everything uses individual parameters)
    if do_drop:
        optdict, llhdict = optimize_ind(model_function, fit_trials, empirical_data, cv_parameters, parameter_names,
                                        initial_params, fixed_params, msprt_samples, max_iters, None, n_reps,
                                        use_full_balance, modtype, use_mass_ratio)
        write_ind_params(os.path.join(
            output_dir, 'baseline.csv'), optdict, llhdict)
        fitllh = -sum(llhdict.values())
        cvllh_dict = _get_ind_cv_llhs(model_function, cv_trials, empirical_data, optdict, parameter_names, msprt_samples,
                                      n_reps, use_full_balance, True, modtype, use_mass_ratio)
        cvllh = sum(cvllh_dict.values())
    else:
        ioptdict, llh = optimize_all(model_function, fit_trials, empirical_data, cv_parameters, parameter_names,
                                     initial_params, fixed_params, msprt_samples, max_iters, None, n_reps,
                                     use_full_balance, modtype, use_mass_ratio)
        fitllh = -llh
        optdict = dict([(w, copy.copy(ioptdict))
                        for w in empirical_data.keys()])
        cvllh = get_all_llh(model_function, cv_trials, empirical_data, msprt_samples=msprt_samples, n_times=n_reps,
                            use_full_balance=use_full_balance, parallelize=True, type=modtype,
                            use_mass_ratio=use_mass_ratio, **ioptdict)

    print "Fit LLH:", fitllh
    print "CV LLH:", cvllh
    print ""

    sfl = open(summary_file_name, 'a')
    sfl.write('baseline,')
    for p in cv_parameters:
        sfl.write('T,')
    sfl.write(str(fitllh) + ',' + str(cvllh) + '\n')
    sfl.close()
    tfl = open(trace_file_name, 'a')
    tfl.write('Baseline fits\nFit LLH: ' + str(fitllh) +
              '\nCV LLH: ' + str(cvllh) + '\n\n')
    tfl.close()

    # Set up triggers for cross-validation
    if do_drop:
        ind_param = copy.copy(cv_parameters)
        agg_param = []
        iter_param = ind_param
    else:
        ind_param = []
        agg_param = copy.copy(cv_parameters)
        iter_param = agg_param
    old_llh = cvllh
    old_name = ""
    old_pardict = optdict
    running = True

    # Run the cross-validation
    while running:
        new_llhs = dict()
        new_pardicts = dict()
        for p in iter_param:
            # Shift a parameter from individual to aggregate
            this_ip = copy.copy(ind_param)
            this_ap = copy.copy(agg_param)
            if do_drop:
                this_ip.remove(p)
                this_ap.append(p)
                this_name = old_name + "NO_" + p
            else:
                this_ip.append(p)
                this_ap.remove(p)
                this_name = old_name + "ADD_" + p

            # Fit this data
            worker_opt = optimize_joint(model_function, fit_trials, empirical_data, this_ap, this_ip,
                                        parameter_names, old_pardict, fixed_params, msprt_samples, max_iters, None,
                                        repeats, n_reps, use_full_balance, modtype, use_mass_ratio)
            this_optdict = dict([(w, worker_opt[w][0])
                                 for w in worker_opt.keys()])
            this_llhdict = dict([(w, worker_opt[w][1])
                                 for w in worker_opt.keys()])
            this_cvllhdict = _get_ind_cv_llhs(model_function, cv_trials, empirical_data, this_optdict, parameter_names,
                                              msprt_samples, n_reps, use_full_balance, True, modtype, use_mass_ratio)
            this_fitllh = sum(this_llhdict.values())
            this_cvllh = sum(this_cvllhdict.values())

            # And store
            new_llhs[p] = this_cvllh
            new_pardicts[p] = this_optdict

            # Write it out
            print "Fits of", this_name
            print "Fit LLH:", this_fitllh
            print "CV LLH:", this_cvllh
            print ""
            write_ind_params(os.path.join(
                output_dir, this_name + '.csv'), this_optdict, this_llhdict)
            sfl = open(summary_file_name, 'a')
            sfl.write(this_name + ',')
            for p in cv_parameters:
                if do_drop:
                    if p in this_ip:
                        sfl.write('T,')
                    else:
                        sfl.write('F,')
                else:
                    if p in this_ap:
                        sfl.write('T,')
                    else:
                        sfl.write('F,')
            sfl.write(str(this_fitllh) + ',' + str(this_cvllh) + '\n')
            sfl.close()
            tfl = open(trace_file_name, 'a')
            tfl.write('Fit of ' + this_name + '\nFit LLH: ' +
                      str(this_fitllh) + '\nCV LLH: ' + str(this_cvllh) + '\n\n')
            tfl.close()

        # Test which of the parameters hurts the most
        t_word = "drop" if do_drop else "add"
        print "Performing", t_word, "test"
        worst_p = ''
        worst_llh = -99999999999
        for p, llh in new_llhs.items():
            if llh > worst_llh:
                worst_p = p
                worst_llh = llh
        print "Easiest parameter to", t_word, ":", worst_p
        print "New CV LLH:", worst_llh
        print "Baseline CV LLH:", old_llh
        tfl = open(trace_file_name, 'a')
        tfl.write('Performing ' + t_word + ' test\nEasiest to ' +
                  t_word + ': ' + worst_p + '\nNewLLH: ' + str(worst_llh))
        tfl.write('\nBase LLH: ' + str(old_llh) + '\n')

        # Is dropping/adding the worst parameter no big deal? Then dos it
        if worst_llh > (old_llh - llh_delta):
            print 'So do it!'
            tfl.write('Changing!\n\n\nNext Iteration:\n\n')
            tfl.close()
            old_llh = worst_llh
            old_name = old_name + "NO_" + worst_p + "___"
            old_pardict = new_pardicts[worst_p]
            if do_drop:
                ind_param.remove(worst_p)
                agg_param.append(worst_p)
                iter_param = ind_param
            else:
                ind_param.append(worst_p)
                agg_param.remove(worst_p)
                iter_param = agg_param
        else:
            print 'We are done!'
            tfl.write("No " + t_word + " -- finishing!\n\n")
            running = False

    tfl.write('Final aggregate parameters:\n\t' + str(agg_param))
    tfl.write('\n\nFinal individual parameters:\n\t' + str(ind_param))
    tfl.write('\n\n')
    tfl.close()

    # Now write out the final set of parameters and fits
    final_llhdict = _get_ind_cv_llhs(model_function, model_trials, empirical_data, old_pardict, parameter_names,
                                     msprt_samples, n_reps, use_full_balance, True, modtype, use_mass_ratio)
    write_ind_params(os.path.join(output_dir, 'cv_params.csv'),
                     old_pardict, final_llhdict)
    write_ind_outcomes(os.path.join(output_dir, 'cv_model_fits.csv'), model_function, empirical_data, old_pardict,
                       model_trials, fittrs, msprt_samples, n_reps, use_full_balance, modtype, use_mass_ratio)


def optimize_for_individual_rules(shapes_trials, materials_trials, balance_trials, shapes_emp, materials_emp,
                                  balance_emp, fit_parameters, parameter_names_dict,
                                  initial_params=None, fixed_params=None, max_iters=500, print_iters=50,
                                  n_reps=500, beam_type='exponential', all_rules=['1', '2', '3', '3a', '4']):
    models = {'sh': rules_shapes_full,
              'mat': rules_materials_full, 'bal': rules_balance_full}
    model_trials_dict = {'sh': shapes_trials,
                         'mat': materials_trials, 'bal': balance_trials}
    empirical_data_dict = {'sh': shapes_emp,
                           'mat': materials_emp, 'bal': balance_emp}

    all_parameter_names = list(
        set([pname for plist in parameter_names_dict.values() for pname in plist]))

    if initial_params is None:
        initial_params = get_init_params(fit_parameters)
    else:
        def_inits = get_init_params(fit_parameters)
        for p in fit_parameters:
            if p not in initial_params:
                initial_params[p] = def_inits[p]
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

    all_bounds = get_param_bounds(fit_parameters)
    firstps = [initial_params[p] for p in fit_parameters]
    bnds = [all_bounds[p] for p in fit_parameters]

    ruletypes = copy.copy(all_rules)
    ruletypes.remove('4')

    def hfnc_f(params):
        tot_llh = 0

        rule_assignments = dict()

        for mod in models:
            parameter_names = parameter_names_dict[mod]
            model_function = models[mod]
            model_trials = model_trials_dict[mod]
            empirical_data = empirical_data_dict[mod]

            use_param = dict([(fp, prm) for fp, prm in zip(
                fit_parameters, params) if fp in parameter_names])
            p = fill_params(use_param, parameter_names, def_params=fill_dict)

            trnms = model_trials.keys()
            wids = empirical_data.keys()

            def getoutcomes(trnm):
                tr = model_trials[trnm]
                ocms = model_function(tr, n_times=n_reps,
                                      beam_type=beam_type, **p)[0]
                return ocms

            def getllh(outcomes, w_empdat):
                return -sum([np.log(outcomes[tr][w_empdat[tr]]) for tr in trnms if tr in w_empdat.keys()])

            #all_outcomes = async_map(getoutcomes, trnms)
            all_outcomes = []
            for tr in trnms:
                all_outcomes.append(getoutcomes(tr))
            rule_outcomes = dict([(r, dict()) for r in all_rules])
            for i, tr in enumerate(trnms):
                for rule in all_rules:
                    rule_outcomes[rule][tr] = all_outcomes[i]['R' + rule]
            # Find the best fitting rule category by individual
            for w in wids:
                edat = empirical_data[w]

                llhs = dict()
                for rule in all_rules:
                    llhs[rule] = getllh(rule_outcomes[rule], edat)

                mlh = min(llhs.values())
                this_r = [r for r in all_rules if llhs[r] == mlh][0]
                tot_llh += mlh
                rule_assignments[w] = this_r

        return tot_llh, rule_assignments

    def hfnc(params):
        return hfnc_f(params)[0]

    optparams, llh, _ = SPSA(
        hfnc, firstps, 1e-4, .1, bounds=bnds, maxiter=max_iters, print_iters=print_iters)
    _, assignments = hfnc_f(optparams)
    optdict = dict([k, p] for k, p in zip(fit_parameters, optparams))
    if fixed_params is not None:
        for p in fixed_params:
            optdict[p] = initial_params[p]
    return optdict, assignments, -llh
    """
