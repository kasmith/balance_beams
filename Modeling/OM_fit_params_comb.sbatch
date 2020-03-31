#!/usr/bin/env bash
#SBATCH --job-name=BalanceModel
#SBATCH --output=SlurmOut/ModelStrat-%j.out
#SBATCH --error=SlurmOut/ModelStrat-%j.err
#SBATCH --mem=32000
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH -t 2-0

python run_combined_strats_all.py -d "$@"
