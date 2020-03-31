from models import *
import argparse, copy, random
from OptimTools import async_map
import numpy as np

parser = argparse.ArgumentParser(description="Run model to determine whether rule assignments are consistent.")
parser.add_argument('-o', '--output', help="Destination file for model output", required=True)
parser.add_argument('-i', '--init_parameters', help="Filename for parameter initialization (optional)",
                    required=True, type=argparse.FileType('rU'))
parser.add_argument('-n', '--n_cv', help='Number of CV repeats', default=1000, type=int)
parser.add_argument('--n_times', help='Number of simulations', default=500, type=int)

ruletypes = ['1','2','3']

def get_trial_results(trial, rule, model, params, ntimes):
    thisp = copy.copy(params)
    for rt in ruletypes:
        if rule == rt:
            thisp['use_rule'+rt] = 1.
        else:
            thisp['use_rule'+rt] = 0.
    return model(trial, n_times=ntimes, beam_type='exponential', **thisp)[0]

def get_all_trial_results(trial_dict, rule, model, params, ntimes):

    def hfnc(trnm):
        tr = trial_dict[trnm]
        return get_trial_results(tr, rule, model, params, ntimes)

    trnms = trial_dict.keys()
    res = async_map(hfnc, trnms)
    return dict([(trnm, r) for trnm, r in zip(trnms, res)])

def getllh(outcomes, worker_emp, trials):
    return -sum([np.log(outcomes[tr][worker_emp[tr]]) for tr in trials if tr in worker_emp.keys()])

def match_rule(rule_predictions, worker_emp, trials):
    rs = ['Rule1','Rule2','Rule3','Rule4']
    llhs = [getllh(rule_predictions[r], worker_emp, trials) for r in rs]
    return rs[llhs.index(min(llhs))]

def split_trials(all_trials):
    n_trials = len(all_trials)
    trs = copy.copy(all_trials.keys())
    random.shuffle(trs)
    fittrs = copy.copy(trs[:int(n_trials / 2)])
    cvtrs = copy.copy(trs[int(n_trials / 2):])

    #fit_trials = dict([(tr, all_trials[tr]) for tr in fittrs])
    #cv_trials = dict([(tr, all_trials[tr]) for tr in cvtrs])

    return fittrs, cvtrs

if __name__ == '__main__':
    args = parser.parse_args()

    fixed_params = ['block_mass_u', 'brick_density', 'brick_u', 'iron_density', 'iron_u']

    # Fit all of the trials using the model
    sh_params = params_mod_2_rules(shapes_param_names, False)
    mat_params = params_mod_2_rules(materials_param_names, False)
    bal_params = params_mod_2_rules(balance_param_names, False)

    sh_wids = shapes_empirical.keys()
    mat_wids = materials_empirical.keys()
    bal_wids = balance_empirical.keys()

    initpars = read_params(args.init_parameters)

    print 'Initialized'

    # Shapes
    sh_outcomes = dict()
    sh_inputs = dict([(pnm, val) for pnm, val in initpars.items() if pnm in sh_params])
    for r in (ruletypes + ['4']):
        sh_outcomes['Rule'+r] = get_all_trial_results(shapes_trials, r, rules_shapes, sh_inputs, args.n_times)
        print 'Done -- Shapes', r

    # Materials
    mat_outcomes = dict()
    mat_inputs = dict([(pnm, val) for pnm, val in initpars.items() if pnm in mat_params])
    for r in (ruletypes + ['4']):
        mat_outcomes['Rule' + r] = get_all_trial_results(materials_trials, r, rules_materials, mat_inputs, args.n_times)
        print 'Done -- Materials', r

    # Balance
    bal_outcomes = dict()
    bal_inputs = dict([(pnm, val) for pnm, val in initpars.items() if pnm in bal_params])
    for r in (ruletypes + ['4']):
        bal_outcomes['Rule' + r] = get_all_trial_results(balance_trials, r, rules_balance, bal_inputs, args.n_times)
        print 'Done -- Balance', r



    # Randomly split the trials and get best rules for each split
    with open(args.output, 'w') as ofl:
        ofl.write('WID,Experiment,Iteration,Set1,Set2\n')

        for itr in xrange(args.n_cv):
            # Shapes
            set1, set2 = split_trials(shapes_trials)
            for w in sh_wids:
                r1 = match_rule(sh_outcomes, shapes_empirical[w], set1)
                r2 = match_rule(sh_outcomes, shapes_empirical[w], set2)
                ofl.write(w + ',Shapes,' + str(itr) + ',' + r1 + ',' + r2 + '\n')

            # Materials
            set1, set2 = split_trials(materials_trials)
            for w in mat_wids:
                r1 = match_rule(mat_outcomes, materials_empirical[w], set1)
                r2 = match_rule(mat_outcomes, materials_empirical[w], set2)
                ofl.write(w + ',Materials,' + str(itr) + ',' + r1 + ',' + r2 + '\n')

            # Balance
            set1, set2 = split_trials(balance_trials)
            for w in bal_wids:
                r1 = match_rule(bal_outcomes, balance_empirical[w], set1)
                r2 = match_rule(bal_outcomes, balance_empirical[w], set2)
                ofl.write(w + ',Balance,' + str(itr) + ',' + r1 + ',' + r2 + '\n')

    print 'All done'

