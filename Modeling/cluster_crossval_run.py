from __future__ import print_function, division
import os
import sys
import h5py
import tempfile
import time
import numpy as np
import argparse
import random
from math import ceil

from run_combined_strats_all import strategy_types

fittypes = ["all", "joint_strats", "joint_percept", "individual"]
def_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                       "output", "strat_crossval")
def_hdf5 = os.path.join(def_dir, "strat_crossval.hdf5")

input_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                         "output", "comb_strats")

shorten = {"all": "A", "joint_strats": "JS", "joint_percept": "JP",
           "individual": "I"}

""" Parse command line arguments """
parser = argparse.ArgumenParser(description="Run cross-validation on models on SLURM")
parser.add_argument('-s', '--strategy', help="Name of strategy collection",
                    default="base_strats", choices=strategy_types.keys())
parser.add_argument('-f', '--fittype', help="Name of the fit type",
                    default="all", choices=fittypes)
parser.add_argument('--hdf', help="Name of hdf5 database",
                    default=def_hdf5)
parser.add_argument('--dbsize', type=int, default=50,
                    help="Number of cross-validation runs to do")


def do_run(fit, strat, db, dbsize, cores=30):
    # Make the input file
    ifnm = os.path.join(input_dir, strat + "_all_params.json")

    # Make and run the batch script
    with tempfile.NamedTemporaryFile(mode="w+") as fl:
        fl.write("#!/usr/bin/env bash\n")
        fl.write("#SBATCH --job-name=BB_CV_" + shorten[fit] + "_" + "strat" + "\n")
        fl.write("#SBATCH --output=SlurmOut/BBCV_" + shorten[fit] + "_" +
                 "strat" + "%A_%a.out\n")
        fl.write("#SBATCH --output=SlurmOut/BBCV_" + shorten[fit] + "_" +
                 "strat" + "%A_%a.err\n")
        fl.write("#SBATCH --mem=32000\n")
        fl.write("#SBATCH -c " + str(cores) + "\n")
        if fit == "individual" or fit == "joint_percept" or strat == "inc_dist":
            fl.write("#SBATCH -t 7-0\n")
        else:
            fl.write("#SBATCH -t 2-0\n")
        fl.write("#SBATCH --array=[1-" + str(dbsize) + "]\n\n")
        fl.write("python run_crossval_strat.py -d -c run -f " + fit +
                 " -s " + strat + " -n $SLURM_ARRAY_TASK_ID -i " + ifnm)
        if strat == "rules" and (fit == "joint_percept" or fit == "individual"):
            fl.write(" --single_strat")
        fl.write(" --hdf " + db)
        fl.write('\n')
        fl.flush()
        fnm = fl.name
        os.system('cat ' + fnm)
        os.system('sbatch + ' fnm)

if __name__ == '__main__':
    args = parser.parse_args()
    do_run(args.fittype, args.strategy, args.hdf, args.dbsize)
