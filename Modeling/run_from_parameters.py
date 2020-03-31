from models import *
import os
import json
import argparse

exp_packages = {
    'shapes': [model_shapes, shapes_trials, shapes_empirical,
               shapes_param_names, rules_shapes],
    'materials': [model_materials, materials_trials, materials_empirical,
                  materials_param_names, rules_materials],
    'pivot': [model_balance, balance_trials, balance_empirical,
              balance_param_names, rules_balance],
    'combined': [model_combined, combined_trials, combined_empirical,
                 combined_param_names, rules_combined],
    'ferretti': [model_clean, ferretti_trials, ferretti_empirical,
                 clean_param_names, rules_clean],
    'geomat': [model_materials, geomat_trials, geomat_empirical,
               materials_param_names, rules_materials]
}

parser = argparse.ArgumentParser(
    description="Fits an experiment from parameters")
parser.add_argument('-e', '--experiment', choices=list(exp_packages.keys()),
                    help="The experiment data to fit", required=True)
parser.add_argument('-p', '--parameters', type=argparse.FileType('rU'),
                    help="The parameter file path", required=True)
parser.add_argument('-o', '--output',
                    help="The output file path", required=True)
parser.add_argument('-n', '--num_sims', type=int, default=500,
                    help="Number of simulations to run")

if __name__ == '__main__':
    args = parser.parse_args()
    model, trials, empdat, param_names, rules_model = \
        exp_packages[args.experiment]

    params = json.load(args.parameters)
    percept_params = params["perceptual"]
    strat_mix = params["strategies"]

    use_params = dict()
    for p in param_names:
        if p in percept_params.keys():
            use_params[p] = percept_params[p]
        else:
            use_params[p] = default_params[p][0]

    write_all_outcomes(args.output, model, empdat, use_params,
                       strat_mix, trials, n_sims=args.num_sims)
