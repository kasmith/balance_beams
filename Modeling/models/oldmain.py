
from shapes import shapes_trials, shapes_empirical
from materials import materials_trials, materials_empirical
from balance import balance_trials, balance_empirical
from geomat import geomat_trials, geomat_empirical
from combined import combined_trials, combined_empirical
from ferretti import ferretti_trials, ferretti_empirical
from model_parts import *
from config import *
from optimization import *
from writeout import *
import argparse, copy

shapes_package = [model_all_shapes, shapes_trials, shapes_empirical, shapes_param_names, rules_all_shapes]
materials_package = [model_all_materials, materials_trials, materials_empirical, materials_param_names, rules_all_materials]
balance_package = [model_all_balance, balance_trials, balance_empirical, balance_param_names, rules_all_balance]
geomat_package = [model_all_materials, geomat_trials, geomat_empirical, materials_param_names, rules_all_materials]
combined_package = [model_all_combined, combined_trials, combined_empirical, combined_param_names, rules_all_combined]
ferretti_package = [model_all_clean, ferretti_trials, ferretti_empirical, clean_param_names, rules_all_clean]

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

modtype_args = ALL_MOD_TYPES

# Set up the parser to take in command line arguments
parser = argparse.ArgumentParser(description="Run model fitting for balance beams models.")
parser.add_argument('-o', '--output', help="Destination file for model output (leave out .csv)", required=True)
parser.add_argument('-t', '--type', help="The type of beams model to fit",
                    default='shapes', choices=['shapes','materials','balance','geomat','combined', 'ferretti'])
parser.add_argument('-f', '--fit', help="Whether to fit all data simultaneously, everyone individual, or joint fitting",
                    default='all', choices=['all','individual','joint','drop_cv','add_cv', 'complete', 'complete_joint',
                                            'rules_complete', 'rules_byind', 'rules'])
parser.add_argument('-i', '--init_parameters', help="Filename for parameter initialization (optional)",
                    default=None, type=argparse.FileType('rU'))
#parser.add_argument('-m', '--model_order', help="Model reductions to excise various contributions",
#                    default='smp', choices=modtype_args)
parser.add_argument('-s', '--strategy', help="Name of strategy collection",
                    default="comb", choices=strategy_types.keys())
parser.add_argument('--phys_samples', type=int, default = 1,
                    help="Number of samples to take for physics MSPRT (defaults to 1)")
parser.add_argument('--beam_type', help = "Whether to calculate beam mass using truncated normal, exponential, or probabilistic",
                    default='trunc_norm', choices=['trunc_norm', 'exponential', 'probabilistic'])
parser.add_argument('--fit_parameters', help="A list of parameter names to fit (by default will be all allowable",
                    nargs='*', default=None)
parser.add_argument('--remove_parameters', help="A list of parameter names to NOT fit (does not work with --fit parameters)",
                    nargs='*', default=None)
parser.add_argument('--individual_parameters', help="A list of parameter names to fit individually (only with joint)",
                    nargs='*', default=['mass_heur_p'])
parser.add_argument('--fixed_parameters', help="A list of parameters that can be initialized from -i but won't be fit",
                    nargs='*', default=None)
parser.add_argument('--no_fit', help="Does not fit the parameters but instead simply runs the model and outputs",
                    action='store_true')
parser.add_argument('--simulations', help="The number of simulations to run", type=int, default=cfg_dict['n_samples'])
parser.add_argument('--max_iter', help="Maximum number of optimization iterations", type=int, default=cfg_dict['max_iter'])
parser.add_argument('--print_iter', help="Number of iterations between printouts in optimization", type=int, default=cfg_dict['print_iter'])
parser.add_argument('--loop_repeats', help="Number of repeats of aggregate/individual fitting (for joint fitting only",
                    type=int, default=cfg_dict['loop_repeats'])
parser.add_argument('--loop_iter', help="Maximum number of optimization iterations for joint fitting", type=int, default=None)
parser.add_argument('--llh_delta', type=float, default=5,
                    help="The slack in allowing a parameter to go from individual to joint (only for drop_cv)")

args = parser.parse_args()

print 'Arguments:'
print args

modtype = args.model_order
umaxdist = True # Here as a placeholder -- may have a switch later


if args.type == 'shapes':
    package = shapes_package
elif args.type == 'materials':
    package = materials_package
elif args.type == 'balance':
    package = balance_package
elif args.type == 'geomat':
    package = geomat_package
elif args.type == 'combined':
    package = combined_package
elif args.type == 'ferretti':
    package = ferretti_package
else:
    raise Exception('Improper package type')

mod_fn, trs, emp, pnames, rule_fn = package
pnames = copy.copy(pnames)

if args.beam_type != 'trunc_norm':
    safe_remove(pnames, 'beam_mass_sd')

if args.beam_type == 'probabilistic' and args.type in ['balance', 'combined']:
    pnames.append('beam_use_prob')

if 'p' not in modtype:
    safe_remove(pnames, ['com_u'])

if 's' not in modtype:
    safe_remove(pnames, 'dist_jnd')

if 'm' not in modtype:
    safe_remove(pnames, 'mass_heur_p')

if 'd' not in modtype:
    safe_remove(pnames, 'dist_heur_p')

if modtype == 'p':
    safe_remove(pnames, 'mass_jnd')

if modtype == 's':
    safe_remove(pnames, ['com_range', 'com_range_vs', 'com_range_s', 'com_range_m', 'com_range_l'])


if args.phys_samples > 1:
    pnames.append('msprt_thresh')

# Set up the simulation
ofl_trials = args.output + '.csv'
ofl_params = args.output + '_params.csv'

nsims = args.simulations
maxiter = args.max_iter
printiter = args.print_iter
use_fb = cfg_dict['use_full_balance']
if args.loop_iter is None:
    loop_iter = maxiter
else:
    loop_iter = args.loop_iter

if args.fit_parameters is None:
    fitpars = copy.copy(pnames)
    if args.remove_parameters is not None:
        for p in args.remove_parameters:
            safe_remove(fitpars, p)
else:
    fitpars = args.fit_parameters

if args.init_parameters is None:
    initpars = None
else:
    initpars = read_params(args.init_parameters)

if args.fixed_parameters is None:
    fixpars = None
else:
    fixpars = args.fixed_parameters
    for p in fixpars:
        if p in fitpars:
            fitpars.remove(p)

if args.fit == 'all':
    if args.no_fit:
        if initpars is None:
            initpars = dict()
        optdict = fill_params(initpars, pnames)
        llh = get_all_llh(mod_fn, trs, emp, use_dist_max=umaxdist, type=modtype, beam_type=args.beam_type, **optdict)
    else:
        optdict, llh = optimize_all(mod_fn, trs, emp, fitpars, pnames, initpars, fixpars, args.phys_samples,
                                      maxiter, printiter, nsims, use_fb, modtype, args.beam_type, umaxdist)
    write_all_params(ofl_params, optdict, llh)
    ps = fill_params(optdict, pnames)
    write_all_outcomes(ofl_trials, mod_fn, emp, ps, trs, None, args.phys_samples, nsims, use_fb, modtype,
                       args.beam_type, umaxdist)
    print 'Fitting complete'

elif args.fit == 'individual':
    optdict, llhdict = optimize_ind(mod_fn, trs, emp, fitpars, pnames, initpars, fixpars, args.phys_samples,
                                    maxiter, printiter, nsims, use_fb, modtype, args.beam_type, umaxdist)
    write_ind_params(ofl_params, optdict, llhdict)
    used_params = dict([(w, fill_params(optdict[w], pnames)) for w in optdict.keys()])
    write_ind_outcomes(ofl_trials, mod_fn, emp, used_params, trs, None, args.phys_samples, nsims, use_fb, modtype,
                       args.beam_type, umaxdist)
    print 'Fitting complete'

elif args.fit == 'joint':
    agg_pars = [p for p in fitpars if p not in args.individual_parameters]
    worker_opt = optimize_joint(mod_fn, trs, emp, agg_pars, args.individual_parameters, pnames,
                                initpars, fixpars, args.phys_samples, loop_iter, printiter, args.loop_repeats, nsims,
                                use_fb, modtype, args.beam_type, umaxdist)
    optdict = dict([(w, worker_opt[w][0]) for w in worker_opt.keys()])
    llhdict = dict([(w, worker_opt[w][1]) for w in worker_opt.keys()])
    write_ind_params(ofl_params, optdict, llhdict)
    used_params = dict([(w, fill_params(optdict[w], pnames)) for w in worker_opt.keys()])
    write_ind_outcomes(ofl_trials, mod_fn, emp, used_params, trs, None, args.phys_samples, nsims, use_fb, modtype,
                       args.beam_type, umaxdist)
    print 'Fitting complete'

elif args.fit == 'rules':
    pnames = params_mod_2_rules(pnames, True)
    if args.no_fit:
        if initpars is None:
            initpars = dict()
        optdict = fill_params(initpars, pnames)
        llh = get_all_llh(rule_fn, trs, emp, beam_type=args.beam_type, **optdict)
    else:
        optdict, llh = optimize_all(rule_fn, trs, emp, fitpars, pnames, initpars, fixpars, args.phys_samples,
                                      maxiter, printiter, nsims, use_fb, modtype, args.beam_type, umaxdist)
    write_all_params(ofl_params, optdict, llh)
    ps = fill_params(optdict, pnames)
    write_all_outcomes(ofl_trials, rule_fn, emp, ps, trs, None, args.phys_samples, nsims, use_fb, modtype,
                       args.beam_type, umaxdist)
    print 'Fitting complete'

elif args.fit in ['drop_cv', 'add_cv']:
    raise NotImplementedError('Changes have been made -- will need to redo these parts')
    do_drop = args.fit == 'drop_cv'
    cross_validate_individual_params(mod_fn, trs, emp, fitpars, pnames, args.output, args.llh_delta,
                                     initpars, fixpars, loop_iter, args.loop_repeats, nsims, use_fb, modtype, umr, do_drop)
    print 'Fitting complete'

elif args.fit == 'rules_byind':
    params = {'sh': params_mod_2_rules(shapes_param_names,False),
              'mat': params_mod_2_rules(materials_param_names,False),
              'bal': params_mod_2_rules(balance_param_names, False)}

    if args.beam_type != 'trunc_norm':
        for pnames in params.values():
            safe_remove(pnames, 'beam_mass_sd')

    if args.beam_type == 'probabilistic':
        for pnames in params.values():
            pnames.append('beam_use_prob')

    if args.fit_parameters is None:
        fitpars = list(set([pname for plist in params.values() for pname in plist]))
        if args.remove_parameters is not None:
            for p in args.remove_parameters:
                safe_remove(fitpars,p)
    else:
        fitpars = args.fit_parameters

    optdict, _, llh = optimize_for_individual_rules(shapes_trials, materials_trials, balance_trials,
                                                 shapes_empirical, materials_empirical, balance_empirical,
                                                 fitpars, params, initpars, fixpars, maxiter, printiter,
                                                 nsims, args.beam_type)

    write_all_params(ofl_params, optdict, llh)

    pnames = params['sh']
    ps = fill_params(dict([(p, optdict[p]) for p in optdict.keys() if p in pnames]), pnames)
    sh_rules = write_ind_rule_outcomes(args.output+'_sh.csv', rules_shapes, shapes_empirical, ps, shapes_trials,
                                       nsims, args.beam_type)

    pnames = params['mat']
    ps = fill_params(dict([(p, optdict[p]) for p in optdict.keys() if p in pnames]), pnames)
    mat_rules = write_ind_rule_outcomes(args.output + '_mat.csv', rules_materials, materials_empirical, ps,
                                        materials_trials, nsims, args.beam_type)

    pnames = params['bal']
    ps = fill_params(dict([(p, optdict[p]) for p in optdict.keys() if p in pnames]), pnames)
    bal_rules = write_ind_rule_outcomes(args.output + '_bal.csv', rules_balance, balance_empirical, ps,
                                        balance_trials, nsims, args.beam_type)

    with open(args.output+'_rules.csv', 'w') as rfl:
        rfl.write('WID,Experiment,Rule\n')
        for w, r in sh_rules.items():
            rfl.write(w + ',Shapes,' + r + '\n')
        for w, r in mat_rules.items():
            rfl.write(w + ',Materials,' + r + '\n')
        for w, r in bal_rules.items():
            rfl.write(w + ',Balance,' + r + '\n')
    print "Fitting done"

elif args.fit == 'complete' or args.fit == 'complete_joint' or args.fit == 'rules_complete':

    if args.fit == 'rules_complete':
        models = {'sh': rules_all_shapes, 'mat': rules_all_materials, 'bal': rules_all_balance}
        params = {'sh': params_mod_2_rules(shapes_param_names),
                  'mat': params_mod_2_rules(materials_param_names),
                  'bal': params_mod_2_rules(balance_param_names)}
    else:
        models = {'sh': model_all_shapes, 'mat': model_all_materials, 'bal': model_all_balance}
        params = {'sh': shapes_param_names, 'mat': materials_param_names, 'bal': balance_param_names}
    trials = {'sh': shapes_trials, 'mat': materials_trials, 'bal': balance_trials}
    empdat = {'sh': shapes_empirical, 'mat': materials_empirical, 'bal': balance_empirical}


    if args.beam_type != 'trunc_norm':
        for pnames in params.values():
            safe_remove(pnames, 'beam_mass_sd')

    if args.beam_type == 'probabilistic':
        for pnames in params.values():
            pnames.append('beam_use_prob')

    if args.fit_parameters is None:
        fitpars = list(set([pname for plist in params.values() for pname in plist]))
        if args.remove_parameters is not None:
            for p in args.remove_parameters:
                safe_remove(fitpars,p)
    else:
        fitpars = args.fit_parameters


    if args.fit == 'complete' or args.fit == 'rules_complete':

        optdict, llh = optimize_complete(models, trials, empdat, fitpars, params, initpars, fixpars, args.phys_samples,
                                         maxiter, printiter, nsims, use_fb, modtype, args.beam_type, umaxdist)

        write_all_params(ofl_params, optdict, llh)
        for m in ['sh', 'mat', 'bal']:
            pnames = params[m]
            ps = fill_params(dict([(p, optdict[p]) for p in optdict.keys() if p in pnames]), pnames)
            write_all_outcomes(args.output + '_' + m + '.csv', models[m], empdat[m], ps, trials[m], None, args.phys_samples,
                               nsims, use_fb, modtype, args.beam_type, umaxdist)
        print 'Fitting complete'

    elif args.fit == 'complete_joint':

        if len(args.individual_parameters) != 1:
            raise Exception('Complete joint fitting requires just a single individual parameter!')
        ipar = args.individual_parameters[0]
        safe_remove(fitpars, ipar)

        worker_opt = optimize_one_joint_complete(models, trials, empdat, ipar, fitpars, params, initpars, fixpars,
                                                      args.phys_samples, maxiter, printiter, args.loop_repeats, nsims,
                                                      use_fb, modtype, args.beam_type, umaxdist, 10)
        optdict = dict([(w, worker_opt[w][0]) for w in worker_opt.keys()])
        llhdict = dict([(w, worker_opt[w][1]) for w in worker_opt.keys()])

        wrp_dict = copy.deepcopy(optdict)
        for m in models:
            emp = empdat[m]
            for w in emp.keys():
                if w in wrp_dict.keys():
                    wrp_dict[w]['Experiment'] = m
                else:
                    print "WID not found:", w, "; from model:", m

        write_ind_params(ofl_params, wrp_dict, llhdict)

        for m in models:
            emp = empdat[m]
            pnames = params[m]
            used_params = dict()
            for w in emp.keys():
                odw = optdict[w]
                u = dict([(pnm, odw[pnm]) for pnm in pnames if pnm in odw.keys()])
                used_params[w] = fill_params(u, pnames)
            write_ind_outcomes(args.output + '_' + m + '.csv', models[m], emp, used_params, trials[m], None,
                               args.phys_samples, nsims, use_fb, modtype, args.beam_type, umaxdist)

        print 'Fitting complete'
