from __future__ import division, print_function
from models import *
import h5py
import filelock
import argparse
import os
import random
import json
import warnings
import numpy as np
import copy
from run_combined_strats_all import strategy_types
from OptimTools import async_map

fitpars = ['com_u', 'com_range', 'beam_mass_mean', 'com_range_s',
           'com_range_vs', 'com_range_m', 'com_range_l', 'mass_jnd',
           'dist_jnd', 'distance_u', 'tower_mass_u']
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

"""Some clean-up to get commonly used variables"""
mod_names = models.keys()
fit_choices = ["all", "joint_strats", "joint_percept", "individual"]
def_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                       "output", "strat_crossval")
def_hdf5 = os.path.join(def_dir, "strat_crossval.hdf5")
def_of_strat = os.path.join(def_dir, "crossval_all_strat.csv")

"""Create the parser for reading in arguments"""
parser = argparse.ArgumentParser(description="Run cross-validation on balance beam models")
parser.add_argument('-s', '--strategy', help="Name of strategy collection",
                    default="comb", choices=strategy_types.keys())
parser.add_argument('-f', '--fit_type', default="all",
                    help="How to fit (by all, individual, joint)",
                    choices=fit_choices)
parser.add_argument('--hdf', help="Name of hdf5 database",
                    default=def_hdf5)
parser.add_argument('-d', '--nodebug', help="Turns off debug mode",
                    action="store_true")
parser.add_argument('-c', '--command', required=True,
                    help="What to do with this cross-validation",
                    choices=["create", "run", "read_llh", "read_params",
                             "read_llh_ind_vs_rules"])
parser.add_argument('-n', '--number', type=int, default=50,
                    help=("Either the size of the database (create) or the" +
                          "index to run cross-validation on"))
parser.add_argument('-i', '--initialization', default=None,
                    type=argparse.FileType('rU'),
                    help="File with initialization parameters")
parser.add_argument('-o', '--output', default=def_of_strat,
                    help="Output csv file (only for read_llh)")
parser.add_argument('--single_strat', action="store_true",
                    help="Whether to fix individuals to use a single strategy")

"""Functions that are used for getting into / out of the hdf5 database"""


def _split_trials(triallist):
    nsamp = int(np.floor(len(triallist)/2.))
    idx_in = random.sample(range(len(triallist)), nsamp)
    items_in = [triallist[i] for i in range(len(triallist)) if i in idx_in]
    items_out = [triallist[i] for i in range(len(triallist)) if i not in idx_in]
    return items_in, items_out

# Shortcut for turning dict keys into something writeable by hdf55
def _keywrite(d):
    return np.array(list(d.keys()), dtype='string_')

def _make_database(filename, size=50):
    assert not os.path.exists(filename), "Database exists -- cannot overwrite"
    with h5py.File(filename, 'w') as f:
        spl = f.create_group("splits")
        fit = f.create_group("fitting")
        f["size"] = size
        # Make the list of cross-validated trials
        for i in range(size):
            si = str(i)
            subspl = spl.create_group(si)
            tlist_fit = subspl.create_group('fit')
            tlist_cv = subspl.create_group('crossval')
            for m in mod_names:
                tin, tout = _split_trials(list(trials[m].keys()))
                tlist_fit[m] = np.array(tin, dtype='string_')
                tlist_cv[m] = np.array(tout, dtype='string_')
        # Make the shell for the fitting
        for fch in fit_choices:
            ch_group = fit.create_group(fch)
            for stype in strategy_types.keys():
                st_group = ch_group.create_group(stype)
                for i in range(size):
                    si = str(i)
                    st_group.create_group(si)
            fg = f['fitting'][fch]
            for stype in strategy_types:
                if stype not in fg.keys():
                    st_group = fg.create_group(stype)
                    for i in range(size):
                        si = str(i)
                        st_group.create_group(si)
    # Close the database and return
    return


def _check_file_consistency(filename, n, fit_choice, strat_type):
    f = h5py.File(filename, 'r')
    assert 0 <= n < f['size'].value
    si = str(n)
    ks = f['fitting'][fit_choice][strat_type][si].keys()
    assert len(ks) == 0, "Data already written to this branch"
    return f


def _retrieve_trial_splits(f, n):
    si = str(n)
    grp = f['splits'][si]
    fit_dict = dict()
    cv_dict = dict()
    for m in mod_names:
        fit_dict[m] = grp['fit'][m].value
        cv_dict[m] = grp['crossval'][m].value
    return fit_dict, cv_dict


def _check_one_database(hdf_file, fit_choice, strat_type):
    g = hdf_file['fitting'][fit_choice][strat_type]
    empties = [i for i, v in g.items() if len(v.keys())==0]
    assert len(empties) != hdf_file['size'].value, "All indices empty!"
    if len(empties) > 0:
        warnings.warn("Empty indices: " + str(empties), RuntimeWarning)
        return empties
    else:
        return []


def _check_completed_database(filename, fit_choice, strat_type):
    f = h5py.File(filename, 'r')
    _check_one_database(f, fit_choice, strat_type)
    return f

"""Reads likelihoods (both fit and cross-validated) for a model"""
def cv_read_llh(filename, fit_choice, strat_type):
    with _check_completed_database(filename, fit_choice, strat_type) as f:
        g = f['fitting'][fit_choice][strat_type]
        fit_llh = []
        cv_llh = []
        for dat in [dat for dat in g.values() if len(dat.keys()) > 0]:
            fit_llh.append(dat['llh_fit'].value)
            cv_llh.append(dat['llh_cv'].value)
    return {'fit_llh': fit_llh, 'cv_llh': cv_llh}

"""Reads likelihood for individual workers (just cv)"""
def cv_read_ind_llh(filename, fit_choice, strat_type):
    assert fit_choice != "all", "Cannot read individual llh from all fitting"
    rets = {}
    with _check_completed_database(filename, fit_choice, strat_type) as f:
        g = f['fitting'][fit_choice][strat_type]
        for i, dat in [[i, dat] for i, dat in g.items() if len(dat.keys()) > 0]:
            for w, d in dat['participants'].items():
                if w not in rets:
                    rets[w] = {}
                try:
                    rets[w][i] = d['llh_cv'].value
                except:
                    print("error not found: ", strat_type, i, w)
    return rets


"""Checks for missing indices"""
def find_all_missing_indices(database, strats=strategy_types.keys(),
                             fit_types=fit_choices):
    f = h5py.File(database, 'r')
    for s in strats:
        ret = _check_one_database(f, "all", s)
        if len(ret) > 0:
            print ("Missing items in strategy: " + s)
            print (ret)
            print ("")
    for ft in fit_types:
        ret = _check_one_database(f, ft, "base_strats")
        if len(ret) > 0:
            print ("Missing items in fit type: " + ft)
            print (ret)
            print ("")

"""Writes out llhs from cross-validated strategies to a csv file"""
def write_cv_llh_all(database, output_file, n=50):
    strats = strategy_types.keys()
    fch = copy.copy(fit_choices)
    fch.remove('all')
    with open(output_file, 'w') as ofl:
        ofl.write("Strategy,FitType,CVOrder,FitLLH,CVLLH\n")
        for s in strats:
            llhs = cv_read_llh(database, "all", s)
            cv_llh = llhs['cv_llh']
            fit_llh = llhs['fit_llh']
            if len(cv_llh) != n:
                print ("Malformed length: " + s)
                print ("Skipping\n")
            else:
                for i in range(n):
                    ofl.write(s + ',all,' + str(i) + ',' + str(fit_llh[i]) +
                              ',' + str(cv_llh[i]) + '\n')
        for ft in fch:
            if ft == 'individual':
                strats2 = ['base_strats']
            else:
                strats2 = ['base_strats', 'rules']
            for s in strats2:
                llhs = cv_read_llh(database, ft, s)
                cv_llh = llhs['cv_llh']
                fit_llh = llhs['fit_llh']
                if len(cv_llh) != n:
                    print ("Malformed length: " + ft)
                    print ("Skipping\n")
                else:
                    for i in range(n):
                        ofl.write(s + ',' + ft + ',' + str(i) + ',' +
                                  str(fit_llh[i]) + ',' + str(cv_llh[i]) + '\n')

def write_cv_llh_ind_vs_rules(database, output_file, fit_type, n=50):
    strats = ['base_strats', 'rules'] # This is only comparing base vs rules
    with open(output_file, 'w') as ofl:
        ofl.write("WID,Strategy,CVOrder,CVLLH\n")
        for s in strats:
            llhs = cv_read_ind_llh(database, fit_type, s)
            for wid, dat in llhs.items():
                for idx, cvllh in dat.items():
                    ofl.write(wid + ',' + s + ',' + idx + ',' +
                              str(cvllh) + '\n')


"""Reads the perceptual parameters from the cv database"""
def cv_read_percept_params(filename, fit_choice, strat_type):
    with _check_completed_database(filename, fit_choice, strat_type) as f:
        g = f['fitting'][fit_choice][strat_type]
        if fit_choice in ['individual', 'joint_percept']:
            participants = dict()
            for dat in [dat for dat in g.values() if len(dat.keys()) > 0]:
                parts = dat['participants']
                for ids, pdat in parts.items():
                    if ids not in participants:
                        participants[ids] = dict()
                        participants[ids]['model'] = pdat['model'].value
                    sdict = participants[ids]
                    for k, v in zip(pdat['percept_keys'], pdat['percept_values']):
                        if k in sdict:
                            sdict[k].append(v)
                        else:
                            sdict[k] = [v]
            return participants
        else:
            params = dict()
            for dat in [dat for dat in g.values() if len(dat.keys()) > 0]:
                for k, v in zip(dat['percept_keys'], dat['percept_values']):
                    if k in params:
                        params[k].append(v)
                    else:
                        params[k] = [v]
            return params

def cv_read_strat_params(filename, fit_choice, strat_type):
    with _check_completed_database(filename, fit_choice, strat_type) as f:
        g = f['fitting'][fit_choice][strat_type]
        if fit_choice in ['individual', 'joint_strats']:
            participants = dict()
            for dat in [dat for dat in g.values() if len(dat.keys()) > 0]:
                parts = dat['participants']
                for ids, pdat in parts.items():
                    if ids not in participants:
                        participants[ids] = dict()
                        participants[ids]['model'] = pdat['model'].value
                    sdict = participants[ids]
                    for k, v in zip(pdat['strat_keys'], pdat['strat_values']):
                        if k in sdict:
                            sdict[k].append(v)
                        else:
                            sdict[k] = [v]
            return participants
        else:
            params = dict()
            for dat in [dat for dat in g.values() if len(dat.keys()) > 0]:
                for k, v in zip(dat['strat_keys'], dat['strat_values']):
                    if k in params:
                        params[k].append(v)
                    else:
                        params[k] = [v]
            return params


"""Functions to lock the file for safe writing"""
def acquire_lock(locknm, callback, timeout=180.):
    lock = filelock.FileLock(locknm)
    try:
        with lock.acquire(timeout=timeout):
            callback()
            lock.release()
            return
    except filelock.Timeout:
        raise IOError("File not available for writing")


"""Functions for applying trial models"""
def _make_trial_models(model_function_dict, trial_dict, perc_params,
                       n_reps=500, n_cores=1):
    ind_stims = []
    for m in model_function_dict.keys():
        parameter_names = params[m]
        model_function = model_function_dict[m]
        trials = trial_dict[m]
        use_param = dict([(fp, prm) for fp, prm in perc_params.items()
                          if fp in parameter_names])
        ind_stims += [(trnm, tr, model_function, use_param) for
                      trnm, tr in trials.items()]

    def _get_itrial(itp):
        nm, tr, fn, param = itp
        return (nm, fn(tr, strategies, n_times=n_reps, **param))

    if n_cores == 1:
        return dict(map(_get_itrial, ind_stims))
    else:
        return dict(async_map(_get_itrial, ind_stims, ncpu=n_cores))

def _get_ind_llh(model_function, trials, perc_params, mix_params, param_names,
                 empdat, n_reps=500, n_cores=1):
    use_param = dict([(fp, prm) for fp, prm in perc_params.items()
                      if fp in param_names])
    use_trials = list(set(trials.keys()) & set(empdat.keys()))
    strategies = list(mix_params.keys())
    if 'guess' in strategies:
        strategies.remove('guess')
    def _get_itrial(trnm):
        mdat = model_function(trials[trnm], strategies, n_times=n_reps,
                              **use_param)
        return get_ind_trial_llh_from_mix(mdat, mix_params, empdat[trnm])

    if n_cores == 1:
        llhs = map(_get_itrial, use_trials)
    else:
        llhs = async_map(_get_itrial, use_trials, ncpu=n_cores)

    return sum(llhs)

"""Run the program"""
if __name__ == "__main__":
    args = parser.parse_args()
    print("Arguments:")
    print(args)

    lock_name = '.'.join(args.hdf.split('.')[:-1]) + "_LOCK"
    strat_nm = args.strategy
    comd = args.command

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
        maxiter = 5
        printiter = 1
        ncore = 1

    if strat_nm == 'rules':
        strategies = ['R1', 'R2', 'R3', 'R3a', 'R4']
        params = {
            'sh': params_mod_2_rules(params['sh'], False),
            'mat': params_mod_2_rules(params['mat'], False),
            'bal': params_mod_2_rules(params['bal'], False)
        }
        models = rule_models # Overwrite the models to use
    else:
        strategies = strategy_types[strat_nm]

    if args.initialization is not None:
        init = json.load(args.initialization)
        iperc = init['perceptual']
        istrat = init['strategies']
    else:
        iperc = None
        istrat = None

    sstate = os.path.join(def_dir, "sstate_"+args.fit_type+"_"+
                          strat_nm+"_"+str(args.number)+".json")

    if comd == "create":
        _make_database(args.hdf, args.number)

    elif comd == "read_llh":
        write_cv_llh_all(args.hdf, args.output)

    elif comd == "read_llh_ind_vs_rules":
        write_cv_llh_ind_vs_rules(args.hdf, args.output, args.fit_type)

    elif comd == "run":
        h5_file = _check_file_consistency(args.hdf, args.number,
                                          args.fit_type, args.strategy)
        fit_trials, cv_trials = _retrieve_trial_splits(h5_file, args.number)
        # Needs to encode as strings
        for m in mod_names:
            fit_trials[m] = np.array(fit_trials[m], dtype='str')
            cv_trials[m] = np.array(cv_trials[m], dtype='str')

        h5_file.close()

        if not args.nodebug:
            for m in mod_names:
                fit_trials[m] = fit_trials[m][:4]
                cv_trials[m] = cv_trials[m][:4]

        fit_tr_dict = dict([(m, {}) for m in mod_names])
        cv_tr_dict = dict([(m, {}) for m in mod_names])
        for m in mod_names:
            for trnm, tr in trials[m].items():
                if trnm in fit_trials[m]:
                    fit_tr_dict[m][trnm] = tr
                elif trnm in cv_trials[m]:
                    cv_tr_dict[m][trnm] = tr

        if args.fit_type == "all":
            perc_params, strat_mix, fit_llh = optimize_complete(models, fit_tr_dict,
                                                            empdat, fitpars,
                                                            params, strategies,
                                                            iperc, istrat,
                                                            max_iters=maxiter,
                                                            print_iters=printiter,
                                                            n_reps=nsims,
                                                            n_cores=ncore,
                                                            savestate=sstate)
            fit_llh *= -1
            cv_tr_mod = _make_trial_models(models, cv_tr_dict, perc_params,
                                           nsims, ncore)
            flat_emp_dat = dict()
            for mod in mod_names:
                flat_emp_dat.update(pack_empirical_data(empdat[mod],
                                                        cv_trials[mod]))
            cv_llh = get_all_trial_llh_from_mix(cv_tr_mod, strat_mix, flat_emp_dat)

            def _write():
                with h5py.File(args.hdf, 'a') as f:
                    g = f['fitting'][args.fit_type][strat_nm][str(args.number)]
                    g['llh_fit'] = fit_llh
                    g['llh_cv'] = cv_llh
                    g['percept_keys'] = _keywrite(perc_params)
                    g['percept_values'] = list(perc_params.values())
                    g['strat_keys'] = _keywrite(strat_mix)
                    g['strat_values'] = list(strat_mix.values())

            acquire_lock(lock_name, _write)

        elif args.fit_type == "joint_percept":
            perc_params, ind_strats, fit_llh = optimize_joint_perc_complete(models,
                                                            fit_tr_dict,
                                                            empdat, fitpars,
                                                            params, strategies,
                                                            iperc, istrat,
                                                            max_iters=maxiter,
                                                            print_iters=printiter,
                                                            n_reps=nsims,
                                                            n_cores=ncore,
                                                            savestate=sstate,
                                                            enforce_single_strat=args.single_strat)
            fit_llh *= -1
            cv_tr_mod = _make_trial_models(models, cv_tr_dict, perc_params,
                                           nsims, ncore)
            cv_by_ind = {}
            flat_cv_trs = [trnm for cvmod in cv_trials.values() for trnm in cvmod]
            for mod in mod_names:
                for wid, emp in empdat[mod].items():
                    utrs = list(set(emp.keys()) & set(flat_cv_trs))
                    cv_by_ind[wid] = 0
                    for tr in utrs:
                        cv_by_ind[wid] += get_ind_trial_llh_from_mix(cv_tr_mod[tr],
                                                             ind_strats[wid]['strat_params'],
                                                             emp[tr])
            cv_llh = sum(cv_by_ind.values())
            def _write():
                with h5py.File(args.hdf, 'a') as f:
                    g = f['fitting'][args.fit_type][strat_nm][str(args.number)]
                    g['llh_fit'] = fit_llh
                    g['llh_cv'] = cv_llh
                    g['percept_keys'] = _keywrite(perc_params)
                    g['percept_values'] = list(perc_params.values())
                    sm_g = g.create_group("participants")
                    for mod in mod_names:
                        for wid in empdat[mod].keys():
                            w_g = sm_g.create_group(wid)
                            w_g['model'] = mod
                            w_g['strat_keys'] = _keywrite(ind_strats[wid]['strat_params'])
                            w_g['strat_values'] = list(ind_strats[wid]['strat_params'].values())
                            w_g['llh_cv'] = cv_by_ind[wid]

            acquire_lock(lock_name, _write)

        elif args.fit_type == "joint_strats":
            ind_perc, strat_mix, fit_llh = optimize_joint_strats_complete(models,
                                                            fit_tr_dict,
                                                            empdat, fitpars,
                                                            params, strategies,
                                                            iperc, istrat,
                                                            max_iters=maxiter,
                                                            print_iters=printiter,
                                                            n_reps=nsims,
                                                            n_cores=ncore,
                                                            savestate=sstate)
            fit_llh *= -1
            cv_tasks = []
            attch_wids = []
            for mod in mod_names:
                iperc_mod = ind_perc[mod]
                for wid, perc_params in iperc_mod.items():
                    cv_tasks.append([models[mod], cv_tr_dict[mod], perc_params,
                                     params[mod], empdat[mod][wid]])
                    attch_wids.append(wid)
            def _do_cv(pars):
                mod_fn, cv_trs, ppars, pnames, edat = pars
                return _get_ind_llh(mod_fn, cv_trs, ppars, strat_mix, pnames,
                                    edat, nsims, 1)
            if ncore == 1:
                icv_llh = list(map(_do_cv, cv_tasks))
            else:
                icv_llh = list(async_map(_do_cv, cv_tasks, ncore))
            cv_by_ind = dict([w, l] for w, l in zip(attch_wids, icv_llh))
            cv_llh = sum(icv_llh)
            def _write():
                with h5py.File(args.hdf, 'a') as f:
                    g = f['fitting'][args.fit_type][strat_nm][str(args.number)]
                    g['llh_fit'] = fit_llh
                    g['llh_cv'] = cv_llh
                    g['strat_keys'] = _keywrite(strat_mix)
                    g['strat_values'] = list(strat_mix.values())
                    sm_g = g.create_group("participants")
                    for mod in mod_names:
                        for wid in empdat[mod].keys():
                            w_g = sm_g.create_group(wid)
                            w_g['model'] = mod
                            w_g['percept_keys'] = _keywrite(ind_perc[mod][wid])
                            w_g['percept_values'] = list(ind_perc[mod][wid].values())
                            w_g['llh_cv'] = cv_by_ind[wid]
            acquire_lock(lock_name, _write)

        elif args.fit_type == "individual":
            ind_perc, ind_strat, fit_llh = optimize_ind(models,
                                                        fit_tr_dict,
                                                        empdat, fitpars,
                                                        params, strategies,
                                                        iperc, istrat,
                                                        max_iters=maxiter,
                                                        print_iters=printiter,
                                                        n_reps=nsims,
                                                        n_cores=ncore,
                                                        savestate=sstate,
                                                        enforce_single_strat=args.single_strat)
            fit_llh *= -1
            cv_tasks = []
            attch_wids = []
            for mod in mod_names:
                iperc_mod = ind_perc[mod]
                for wid, perc_params in iperc_mod.items():
                    cv_tasks.append([models[mod], cv_tr_dict[mod], perc_params,
                                     ind_strat[wid], params[mod],
                                     empdat[mod][wid]])
                    attch_wids.append(wid)
            def _do_cv(pars):
                mod_fn, cv_trs, ppars, pstrat, pnames, edat = pars
                return _get_ind_llh(mod_fn, cv_trs, ppars, pstrat, pnames,
                                    edat, nsims, 1)
            if ncore == 1:
                icv_llh = map(_do_cv, cv_tasks)
            else:
                icv_llh = async_map(_do_cv, cv_tasks, ncore)
            icv_llh = list(icv_llh)
            cv_by_ind = dict([w, l] for w, l in zip(attch_wids, icv_llh))
            cv_llh = sum(icv_llh)
            def _write():
                with h5py.File(args.hdf, 'a') as f:
                    g = f['fitting'][args.fit_type][strat_nm][str(args.number)]
                    g['llh_fit'] = fit_llh
                    g['llh_cv'] = cv_llh
                    sm_g = g.create_group("participants")
                    for mod in mod_names:
                        for wid in empdat[mod].keys():
                            #if wid == 'A10DEO061A6L3O:33CUSNVVNOWUUPE407HDJUVPXQV883':
                            #    import pdb; pdb.set_trace
                            w_g = sm_g.create_group(wid)
                            w_g['model'] = mod
                            w_g['percept_keys'] = _keywrite(ind_perc[mod][wid])
                            w_g['percept_values'] = list(ind_perc[mod][wid].values())
                            w_g['strat_keys'] = _keywrite(ind_strat[wid])
                            w_g['strat_values'] = list(ind_strat[wid].values())
                            w_g['llh_cv'] = cv_by_ind[wid]
            acquire_lock(lock_name, _write)
