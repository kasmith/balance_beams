# Fit parameters

### Start by fitting the default parameters for just the base strategies for E1-3:

`sbatch OM_fit_params_comb.sbatch -s base_strats -f all`

### Fit other parameters to E1-3 using the base strats as the baseline:

`sbatch OM_fit_params_comb.sbatch -s comb -f all -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s inc_dist -f all -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s no_phys -f all -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s no_sym -f all -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s no_weight -f all -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s just_phys -f all -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s just_smp -f all -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s just_sp -f all -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s rules -f all -i output/comb_strats/base_strats_all_params.json`

### Fit at more individual levels:

`sbatch OM_fit_params_comb.sbatch -s base_strats -f joint_strats -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s base_strats -f joint_percept -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s base_strats -f individual -i output/comb_strats/base_strats_all_params.json`

`sbatch OM_fit_params_comb.sbatch -s rules -f joint_strats -i output/comb_strats/rules_all_params.json --single_strat`


### Set up and run crossvalidation

(NOTE: very computationally expensive... set up to run on SLURM cluster)

`python run_crossval_strat.py -c create -n 50`

For each strategy {strat} in [base_strats, comb, inc_dist, no_phys, no_sym, no_weight, just_sp, just_smp, just_phys, rules]:

`python cluster_crossval_run.py -s {strat}`

For each fit type {ft} in [joint_percept, joint_strats, individual]:

`python cluster_crossval_run.py -f {ft}`
`python cluster_crossval_run.py -f {ft} -s rules`

(Possibly skip individual, rules... might take too long...)

Now run the 'add1' parameterizations -- for {strat} in {sm, spm, s, msp, ms, mps, mp, m, psm, ps, pms, pm, p}:

`python cluster_crossval_run.py -b -s add1_{strat}`


### Run additional experiment analyses (much faster -- don't need cluster)

`python run_from_parameters.py -i output/comb_strats/base_strats_all_params.json -o output/extensions_from_comb/base_strats_ferretti.csv -t ferretti -s base_strats`

`python run_from_parameters.py -i output/comb_strats/base_strats_all_params.json -o output/extensions_from_comb/base_strats_geomat.csv -t geomat -s base_strats`

`python run_from_parameters.py -i output/comb_strats/base_strats_all_params.json -o output/extensions_from_comb/base_strats_combined.csv -t combined -s base_strats`

`python run_from_parameters.py -i output/comb_strats/rules_all_params.json -o output/extensions_from_comb/rules_ferretti.csv -t ferretti -s rules`

`python run_from_parameters.py -i output/comb_strats/rules_all_params.json -o output/extensions_from_comb/rules_geomat.csv -t geomat -s rules`

`python run_from_parameters.py -i output/comb_strats/rules_all_params.json -o output/extensions_from_comb/rules_combined.csv -t combined -s rules`

`python run_rationality.py`
