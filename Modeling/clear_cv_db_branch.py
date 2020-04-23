import h5py
import argparse
import os
from run_crossval_strat import def_hdf5, strategy_types, fit_choices


def clear_branch(dbnm, fit_choice, strat_type):
    with h5py.File(dbnm, 'a') as fl:
        g = fl['fitting'][fit_choice][strat_type]
        for k in g.keys():
            del g[k]
            g.create_group(k)


parser = argparse.ArgumentParser(description="Clear branches from db when errors are made")
parser.add_argument('-s', '--strategy', help="Name of strategy collection",
                    choices=strategy_types.keys(), required=True)
parser.add_argument('-f', '--fit_type', required=True,
                    help="How to fit (by all, individual, joint)",
                    choices=fit_choices)
parser.add_argument('--hdf', help="Name of hdf5 database",
                    default=def_hdf5)

if __name__ == '__main__':
    args = parser.parse_args()
    clear_branch(args.hdf, args.fit_type, args.strategy)
