# balance_beams
Code for running and modeling balance beam experiments

## Analysis

All analyses in the main paper and appendix can be recreated by running the `Analysis/BalBeamOutline.R` script

## Experiments

Each experiment can be found as a separate subfolder inside of `Experiments`. These are each set up as separate [psiTurk](https://psiturk.org/) projects and thus require the psiTurk framework to run.

## Modeling

The output of all models and aggregation with human data is included in this repo and can be found in the `Modeling/anonymized_output` folder.

If you wish to rerun the models, all code to do so is in the `Modeling` folder. To install dependencies run the following commands:

> pip install -r requirements.txt

Further notes can be found in the README in the `Modeling` folder.