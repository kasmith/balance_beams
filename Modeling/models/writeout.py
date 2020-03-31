from __future__ import division
import numpy as np
from .model_parts import pack_empirical_data, make_trial_models, apply_strategy_selection
from .config import fill_params
import copy
from OptimTools import async_map

def write_all_params(flnm, paramdict, llh = 0):
    fl = open(flnm,'w')
    pnames = paramdict.keys()
    for pn in pnames:
        fl.write(pn + ',')
    fl.write('LLH\n')
    for pn in pnames:
        fl.write(str(paramdict[pn]) + ',')
    fl.write(str(llh))
    fl.write('\n')
    fl.close()

def write_ind_params(flnm, worker_opts, worker_llhs):
    fl = open(flnm,'w')
    firstwrite = True
    for wid in worker_opts.keys():
        thisdict = worker_opts[wid]
        thisllh = worker_llhs[wid]
        if firstwrite:
            pnames = thisdict.keys()
            fl.write('WID,')
            for pn in pnames:
                fl.write(pn + ',')
            fl.write('LLH\n')
            firstwrite = False
        fl.write(wid + ',')
        for pn in pnames:
            fl.write(str(thisdict[pn]) + ',')
        fl.write(str(thisllh))
        fl.write('\n')
    fl.close()


def _write_single_out(ofl, model_function, empdat, perc_param_dict,
                      strat_mix_dict, model_trials, fit_trials, n_sims=500,
                      wid=None, model_trials_as_fits=False):
    #trs = [model_trials[tr] for tr in model_trials.keys()]
    strats = list(strat_mix_dict.keys())
    strats.remove('guess')
    if model_trials_as_fits:
        tr_mods = model_trials
    else:
        tr_mods = make_trial_models(model_function, model_trials,
                                    perc_param_dict, strats, n_times=n_sims)
    tr_names = list(set(model_trials.keys()) & set(empdat.keys()))
    for tr in tr_names:
        mo = tr_mods[tr]
        mpct = apply_strategy_selection(mo, strat_mix_dict)
        edat = empdat[tr]
        if wid is not None:
            ofl.write(wid + ',')
            udat = {'L': 0, 'R': 0, 'B': 0}
            udat[edat] += 1
            edat = udat
        llh = edat['R']*np.log(mpct['R']) + edat['L']*np.log(mpct['L']) + edat['B']*np.log(mpct['B'])
        ofl.write(tr + ',')
        if tr in fit_trials:
            ofl.write('Y,')
        else:
            ofl.write('N,')
        ofl.write(str(edat['L']) + ',' + str(edat['B']) + ',' +
                  str(edat['R']) + ',' + str(mpct['L']) + ',' + str(mpct['B']) +
                  ',' + str(mpct['R']) + ',' + str(llh) + '\n')

def _write_single_halfstrat(ofl, model_function, empdat, perc_param_dict,
                      strat_mix_dict, model_trials, fit_trials, n_sims=500,
                      wid=None, model_trials_as_fits=False, is_first_half=True):
    #trs = [model_trials[tr] for tr in model_trials.keys()]
    strats = strat_mix_dict.keys()
    strats.remove('guess')
    if model_trials_as_fits:
        tr_mods = model_trials
    else:
        tr_mods = make_trial_models(model_function, model_trials,
                                    perc_param_dict, strats, n_times=n_sims)
    tr_names = list(set(model_trials.keys()) & set(empdat.keys()))
    for tr in tr_names:
        mo = tr_mods[tr]
        mpct = apply_strategy_selection(mo, strat_mix_dict)
        edat = empdat[tr]
        if wid is not None:
            ofl.write(wid + ',')
            udat = {'L': 0, 'R': 0, 'B': 0}
            udat[edat] += 1
            edat = udat
        llh = edat['R']*np.log(mpct['R']) + edat['L']*np.log(mpct['L']) + edat['B']*np.log(mpct['B'])
        ofl.write(tr + ',')
        if tr in fit_trials:
            ofl.write('Y,')
        else:
            ofl.write('N,')
        ofl.write(str(edat['L']) + ',' + str(edat['B']) + ',' +
                  str(edat['R']) + ',' + str(mpct['L']) + ',' + str(mpct['B']) +
                  ',' + str(mpct['R']) + ',' + str(llh) + ',')
        if is_first_half:
            ofl.write('first\n')
        else:
            ofl.write('second\n')


def write_all_outcomes(filename, model_function, empirical_data,
                       perc_param_dict, strat_mix_dict, model_trials,
                       fit_trials=None, n_sims=500):
    if fit_trials is None:
        fit_trials = model_trials.keys()
    empdat = pack_empirical_data(empirical_data, model_trials.keys())
    with open(filename, 'w') as ofl:
        ofl.write("Trial,WasFit,EmpLeft,EmpBal,EmpRight,ModLeft,ModBal," +
                  "ModRight,LLH\n")
        _write_single_out(ofl, model_function, empdat, perc_param_dict,
                          strat_mix_dict, model_trials, fit_trials, n_sims)


def write_joint_strat_outcomes(filename, model_function, empirical_data,
                               perc_param_dict, ind_strat_mix_dict,
                               model_trials, fit_trials=None, n_sims=500):
    if fit_trials is None:
        fit_trials = model_trials.keys()
    strats = ind_strat_mix_dict.values()[0]['strat_params'].keys()
    strats.remove('guess')
    tr_mods = make_trial_models(model_function, model_trials, perc_param_dict,
                                strats, n_times=n_sims)
    with open(filename, 'w') as ofl:
        ofl.write("WID,Trial,WasFit,EmpLeft,EmpBal,EmpRight,ModLeft,ModBal," +
                  "ModRight,LLH\n")
        for wid, indempdat in empirical_data.items():
            mix = ind_strat_mix_dict[wid]['strat_params']
            _write_single_out(ofl, model_function, indempdat, perc_param_dict,
                              mix, tr_mods, fit_trials, n_sims, wid, True)

def write_half_strat_outcomes(filename, model_function, empirical_data,
                               perc_param_dict, ind_strat_mix_dict, order_dict,
                               model_trials, fit_trials=None, n_sims=500):
    if fit_trials is None:
        fit_trials = model_trials.keys()
    strats = ind_strat_mix_dict.values()[0]['strat_params_first'].keys()
    strats.remove('guess')
    tr_mods = make_trial_models(model_function, model_trials, perc_param_dict,
                                strats, n_times=n_sims)
    with open(filename, 'w') as ofl:
        ofl.write("WID,Trial,WasFit,EmpLeft,EmpBal,EmpRight,ModLeft,ModBal," +
                  "ModRight,LLH,HalfOrder\n")
        for wid, indempdat in empirical_data.items():
            mix_0 = ind_strat_mix_dict[wid]['strat_params_first']
            mix_1 = ind_strat_mix_dict[wid]['strat_params_second']
            ord = order_dict[wid]
            idx_cutoff = int(len(model_trials.keys()) / 2)
            emp_0 = dict([(tr, r) for tr, r in indempdat.items()
                          if ord[tr] < idx_cutoff])
            emp_1 = dict([(tr, r) for tr, r in indempdat.items()
                          if ord[tr] >= idx_cutoff])
            _write_single_halfstrat(ofl, model_function, emp_0, perc_param_dict,
                              mix_0, tr_mods, fit_trials, n_sims, wid, True,
                              is_first_half=True)
            _write_single_halfstrat(ofl, model_function, emp_1, perc_param_dict,
                              mix_1, tr_mods, fit_trials, n_sims, wid, True,
                              is_first_half=False)


def write_ind_outcomes(filename, model_function, empirical_data,
                               ind_perc_param_dict, ind_strat_mix_dict,
                               model_trials, fit_trials=None, n_sims=500):
    if fit_trials is None:
        fit_trials = model_trials.keys()
    strats = ind_strat_mix_dict.values()[0].keys()
    strats.remove('guess')
    with open(filename, 'w') as ofl:
        ofl.write("WID,Trial,WasFit,EmpLeft,EmpBal,EmpRight,ModLeft,ModBal," +
                  "ModRight,LLH\n")
        for wid, indempdat in empirical_data.items():
            mix = ind_strat_mix_dict[wid]
            _write_single_out(ofl, model_function, indempdat,
                              ind_perc_param_dict[wid], mix, model_trials,
                              fit_trials, n_sims, wid)


'''
def write_ind_outcomes(filename, model_function, empirical_data, parameter_dict, model_trials, fit_trials = None,
                       msprt_samples=1, n_sims = 500, use_full_balance = True, modtype = 'full',
                       beam_type = 'trunc_norm', use_dist_max = True):
    if fit_trials is None:
        fit_trials = model_trials.keys()
    ofl = open(filename, 'w')
    ofl.write("WID,Trial,WasFit,EmpLeft,EmpBal,EmpRight,ModLeft,ModBal,ModRight,LLH,SymmetryHeuristic,MassHeuristic,SimLeft,SimBal,SimRight\n")
    for w in empirical_data.keys():
        edat = dict()
        for tr in model_trials.keys():
            if tr in empirical_data[w].keys():
                edat[tr] = {'R':0,'L':0,'B':0}
                edat[tr][empirical_data[w][tr]] += 1
        params = parameter_dict[w]
        _write_single_out(ofl, model_function, edat, params, model_trials, fit_trials, msprt_samples, w,
                          n_sims, use_full_balance, modtype, beam_type, use_dist_max)
    ofl.close()
'''

def write_ind_rule_outcomes(filename, model_function, empirical_data, parameter_dict, model_trials, n_sims=500,
                        beam_type = 'exponential'):

    trnms = model_trials.keys()
    wids = empirical_data.keys()
    ruletypes = ['1', '2', '3', '3a']
    #ruletypes = ['1', '2', '3']
    rule_dict = dict()

    def getoutcomes(trnm, rule):
        tr = model_trials[trnm]
        thisp = copy.copy(parameter_dict)
        for rt in ruletypes:
            if rule == rt:
                thisp['use_rule' + rt] = 1.
            else:
                thisp['use_rule' + rt] = 0.
        return model_function(tr, n_times=n_sims, beam_type=beam_type, **thisp)[0]

    def getllh(outcomes, w_empdat):
        return -sum([np.log(outcomes[tr][w_empdat[tr]]) for tr in trnms if tr in w_empdat.keys()])

    r1_out = async_map(lambda t: getoutcomes(t, '1'), trnms)
    r1_out = dict([(tr, o) for tr, o in zip(trnms, r1_out)])
    r2_out = async_map(lambda t: getoutcomes(t, '2'), trnms)
    r2_out = dict([(tr, o) for tr, o in zip(trnms, r2_out)])
    r3_out = async_map(lambda t: getoutcomes(t, '3'), trnms)
    r3_out = dict([(tr, o) for tr, o in zip(trnms, r3_out)])
    r3a_out = async_map(lambda t: getoutcomes(t, '3a'), trnms)
    r3a_out = dict([(tr, o) for tr, o in zip(trnms, r3a_out)])
    r4_out = async_map(lambda t: getoutcomes(t, '4'), trnms)
    r4_out = dict([(tr,o) for tr, o in zip(trnms, r4_out)])

    with open(filename,'w') as ofl:

        ofl.write("WID,Trial,EmpLeft,EmpBal,EmpRight,ModLeft,ModBal,ModRight,LLH,Rule\n")
        for w in wids:
            edat = empirical_data[w]
            lh_1 = getllh(r1_out, edat)
            lh_2 = getllh(r2_out, edat)
            lh_3 = getllh(r3_out, edat)
            lh_3a = getllh(r3a_out, edat)
            lh_4 = getllh(r4_out, edat)

            lhs = [lh_1, lh_2, lh_3, lh_3a, lh_4]
            #lhs = [lh_1, lh_2, lh_3, lh_4]
            mlh = min(lhs)

            if mlh == lh_1:
                rule = 'Rule1'
                ocms = r1_out
            elif mlh == lh_2:
                rule = 'Rule2'
                ocms = r2_out
            elif mlh == lh_3:
                rule = 'Rule3'
                ocms = r3_out
            elif mlh == lh_3a:
                rule = 'Rule3a'
                ocms = r3a_out
            else:
                rule = 'Rule4'
                ocms = r4_out

            rule_dict[w] = rule

            for tr in trnms:
                if tr in edat.keys():
                    ofl.write(w + ',' + tr + ',')
                    if edat[tr] == 'L':
                        ofl.write('1,')
                    else:
                        ofl.write('0,')
                    if edat[tr] == 'B':
                        ofl.write('1,')
                    else:
                        ofl.write('0,')
                    if edat[tr] == 'R':
                        ofl.write('1,')
                    else:
                        ofl.write('0,')
                    ofl.write(str(ocms[tr]['L']) + ',' + str(ocms[tr]['B']) + ',' + str(ocms[tr]['R']) + ',')
                    ofl.write(str(np.log(ocms[tr][edat[tr]])) + ',' + rule + '\n')

    return rule_dict


__all__ = ['write_all_params', 'write_ind_params', 'write_all_outcomes', 'write_ind_outcomes',
           'write_ind_rule_outcomes', 'write_half_strat_outcomes']
