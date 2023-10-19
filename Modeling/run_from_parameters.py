from models import *
import argparse, copy
import json

"""
    This model __main__ doesn't do fitting, but can take parameters and run the model for a given setting
"""


shapes_package = [model_shapes, shapes_trials, shapes_empirical, shapes_param_names, rules_shapes]
materials_package = [model_materials, materials_trials, materials_empirical, materials_param_names, rules_materials]
balance_package = [model_balance, balance_trials, balance_empirical, balance_param_names, rules_balance]
geomat_package = [model_materials, geomat_trials, geomat_empirical, materials_param_names, rules_materials]
combined_package = [model_combined, combined_trials, combined_empirical, combined_param_names, rules_combined]
ferretti_package = [model_clean, ferretti_trials, ferretti_empirical, clean_param_names, rules_clean]
learnbene_package = [model_clean, learnbene_trials, learnbene_empirical, clean_param_names, rules_clean]

strategy_types = {
    'comb': MSP_MOD_TYPES,
    'inc_dist': ALL_MOD_TYPES,
    'base_strats': ['smp', 'sp'],
    'no_phys': make_type_defs('sm'),
    'no_sym': make_type_defs('mp'),
    'no_weight': make_type_defs('sp'),
    'just_sp': ['sp'],
    'just_smp': ['smp'],
    'rules': [],
    'just_phys': ['p']
}

for stp in MSP_MOD_TYPES:
    if stp not in ['sp', 'smp']:
        strategy_types['add1_'+stp] = ['sp', 'smp', stp]

modtype_args = ALL_MOD_TYPES

# Set up the parser to take in command line arguments
parser = argparse.ArgumentParser(description="Run model fitting for balance beams models.")
parser.add_argument('-o', '--output', help="Destination file for model output",
                    default=None, type=str)
parser.add_argument('-i', '--init_parameters', help="Filename for parameter initialization (optional)",
                    required=True, type=argparse.FileType('r'))
parser.add_argument('-t', '--type', help="The type of beams model to fit",
                    default='shapes', choices=['shapes','materials','balance',
                                               'geomat','combined', 'ferretti',
                                               'learnbene'])
parser.add_argument('-s', '--strategy', help="Name of strategy collection",
                    default="comb", choices=strategy_types.keys())
parser.add_argument('--simulations', help="The number of simulations to run", type=int, default=cfg_dict['n_samples'])
args = parser.parse_args()

if __name__ == '__main__':
    print ('Arguments:')
    print (args)

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
    elif args.type == 'learnbene':
        package = learnbene_package
    else:
        raise Exception('Improper package type')

    mod_fn, trs, emp, pnames, rule_fn = package
    pnames = copy.copy(pnames)

    if args.strategy == 'rules':
        fn = rule_fn
        pnames = params_mod_2_rules(pnames, False)
    else:
        fn = mod_fn

    safe_remove(pnames, 'beam_mass_mean')

    # Set up the simulation
    ofl_trials = args.output
    nsims = args.simulations

    initpars = json.load(args.init_parameters)
    strat_pars = initpars['strategies']
    perc_pars = initpars['perceptual']

    use_perc = {}
    for p in pnames:
        use_perc[p] = perc_pars.get(p, default_params[p][0])

    write_all_outcomes(ofl_trials, fn, emp, use_perc, strat_pars, trs, n_sims=nsims)
