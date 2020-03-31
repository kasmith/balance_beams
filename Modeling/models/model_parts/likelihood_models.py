from __future__ import division
import numpy as np
from .common_parts import ALL_MOD_TYPES

##########
#
# Calculate individual likelihoods
#
##########

def pack_empirical_data(emp_dat, trnames):
    choices = dict([(tr, {'L':0, 'B':0, 'R':0}) for tr in trnames])
    for w in emp_dat.keys():
        for tr in trnames:
            if tr in emp_dat[w].keys():
                choices[tr][emp_dat[w][tr]] += 1
    return choices


def apply_strategy_selection(trial_model, strategy_mix):
    l = 0.
    b = 0.
    r = 0.
    for sk, spct in strategy_mix.items():
        if sk == 'guess':
            l += 1./3. * spct
            b += 1./3. * spct
            r += 1./3. * spct
        else:
            l += trial_model[sk][0]['L'] * spct
            b += trial_model[sk][0]['B'] * spct
            r += trial_model[sk][0]['R'] * spct
    return {'L': l, 'R': r, 'B': b}

def get_one_trial_llh_from_mix(trial_model, strategy_mix, emp_dat, prior=.01):
    tot_choice = apply_strategy_selection(trial_model, strategy_mix)
    return sum([emp_dat[ch]*np.log((tot_choice[ch] + prior) / (1.+3.*prior))
                for ch in ['L', 'B', 'R']])

def get_ind_trial_llh_from_mix(trial_model, strategy_mix, emp_dat, prior=.01):
    tot_choice = apply_strategy_selection(trial_model, strategy_mix)
    return np.log((tot_choice[emp_dat] + prior) / (1.+3.*prior))

def get_all_trial_llh_from_mix(trial_models, strategy_mix, all_emp_dat,
                               prior=.01):
    return sum([get_one_trial_llh_from_mix(trial_models[tr],
                                           strategy_mix,
                                           all_emp_dat[tr],
                                           prior)
                for tr in trial_models.keys()])

def get_one_trial_llh_by_strat(trial_model, strategy_list, emp_dat, prior=.01):
    llh_per = {}
    for strat in strategy_list:
        choices = trial_model[strat][0]
        llh_per[strat] = sum([emp_dat[ch]*np.log((choices[ch] + prior) /
                                                 (1.+ 3.*prior))
                              for ch in ['L', 'B', 'R']])
    return llh_per

def get_ind_trial_llh_by_strat(trial_model, strategy_list, emp_dat, prior=.01):
    llh_per = {}
    for strat in strategy_list:
        llh_per[strat] = np.log((choices[empdat] + prior) / 1. + 3.*prior)
    return llh_per

def get_all_trial_llh_by_strat(trial_models, strategy_list, all_emp_dat,
                               prior=.01):
    llh_per = dict([(s, 0.) for s in strategy_list])
    for tr in trial_models.keys():
        next_llh = get_one_trial_llh_by_strat(trial_models[tr],
                                              strategy_list,
                                              all_emp_dat[tr],
                                              prior)
        for s in strategy_list:
            llh_per[s] += next_llh[s]
    return llh_per

def get_best_llh_by_single_strat(trial_models, strategy_list, all_emp_dat,
                                 prior=.01):
    llhs = get_all_trial_llh_by_strat(trial_models, strategy_list, all_emp_dat,
                                      prior)
    best_strat = strategy_list[0]
    best_llh = llhs[best_strat]
    for s in strategy_list[1:]:
        this_llh = llhs[s]
        if this_llh < best_llh:
            best_strat = s
            best_llh = this_llh
    return best_llh, best_strat


def make_trial_models(model_function, model_trials, percept_params,
                      model_types=ALL_MOD_TYPES, n_times=500):
    return dict([(tr, model_function(model_trials[tr], strats=model_types,
                                     n_times=n_times, **percept_params))
                 for tr in model_trials.keys()])


def get_all_llh(model_function, model_trials, empirical_data, percept_params={},
                strategy_params={'guess': 1.}):
    assert np.isclose(sum(strategy_params.values()), 1.), \
        "Stragegy mix must sum to 1"
    trnms = model_trials.keys()

    # Load the empirical data
    tr_chs = pack_empirical_data(empirical_data, trnms)

    trial_models = make_trial_models(model_function, model_trials,
                                     percept_params)
    return get_all_trial_llh_from_mix(trial_models, strategy_params,
                                      empirical_data)


def get_ind_llh(model_function, model_trials, ind_data, fit_all_trials=True, **kwdargs):
    trnms = model_trials.keys()
    if not fit_all_trials:
        trnms = [t for t in trnms if t in ind_data.keys()]
    usetrs = [model_trials[tr] for tr in trnms]
    modout = model_function(usetrs, **kwdargs)
    modps = [m[0] for m in modout]
    llh = 0
    for mps, tr in zip(modps, trnms):
        if tr in ind_data.keys():
            edat = ind_data[tr]
            llh += np.log(mps[edat])
    return llh
