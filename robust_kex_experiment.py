# function for running experiments with different robust kidney exchange formulations

import argparse
import time
import numpy as np

from kidney_graph_io import get_cmu_graphs
from utils import generate_filepath, realize_edge_existence, initialize_edges, score_after_edge_realization, \
    set_edge_realizations
from kidney_ip import OptConfig, optimize_picef, optimize_picef_heterogeneous_edge_prob, optimize_dro_saa_edge_existence


def robust_kex_experiment(args):
    rs = np.random.RandomState(seed=args.seed)

    # generate an output file
    output_file = generate_filepath(args.output_dir, 'edge_existence_robust_kex_experiment', 'csv')

    # list of output column names
    output_columns = ['graph_name',
                      'trial_num',
                      'realization_num',
                      'cycle_cap',
                      'chain_cap',
                      'method',
                      'realized_score',
                      'runtime'
                      ]

    # output file header, and write experiment parameters to the first line
    with open(output_file, 'w') as f:
        f.write(str(args) + '\n')
        f.write((','.join(len(output_columns) * ['%s']) + '\n') % tuple(output_columns))

        graph_generator = get_cmu_graphs(args.input_dir)

    # run the experiment for each graph
    for digraph, ndd_list, graph_name in graph_generator:

        print("running tests for graph: %s" % graph_name)

        for i_trial in range(args.num_trials):

            # store all solutions (matchings) in a dict
            sol_dict = {}
            time_dict = {}

            # initialize edge success probabilities
            initialize_edges(digraph, ndd_list, rs)

            # method 1: solve the non-robust approach with PICEF
            if args.use_nonrobust:
                t0 = time.time()
                sol_dict['nonrobust'] = optimize_picef(OptConfig(digraph, ndd_list, args.cycle_cap,
                                                                 args.chain_cap,
                                                                 verbose=args.verbose))
                time_dict['nonrobust'] = time.time() - t0

            # method 1a: solve the non-robust approach with a uniform success probability (still PICEF)
            mean_edge_prob = 0.5
            if args.use_nonrobust_fixed_p:
                t0 = time.time()
                sol_dict['nonrobust_fixed_p'] = optimize_picef(OptConfig(digraph, ndd_list, args.cycle_cap,
                                                                         args.chain_cap,
                                                                         verbose=args.verbose,
                                                                         edge_success_prob=mean_edge_prob))
                time_dict['nonrobust_fixed_p'] = time.time() - t0

            # method 2: solve inhomogeneous edge prob. PICEF formulation
            if args.use_nonrobust_variable_p:
                t0 = time.time()
                sol_dict['nonrobust_variable_p'] = optimize_picef_heterogeneous_edge_prob(
                    OptConfig(digraph, ndd_list, args.cycle_cap,
                              args.chain_cap,
                              verbose=args.verbose)
                )
                time_dict['nonrobust_variable_p'] = time.time() - t0

            # method 3: solve using DRO-SAA PICEF formulation
            if args.use_saa:
                t0 = time.time()
                set_edge_realizations(digraph, ndd_list, args.num_measurements, rs)
                sol_dict['saa'] = optimize_dro_saa_edge_existence(
                    OptConfig(digraph, ndd_list, args.cycle_cap,
                              args.chain_cap,
                              verbose=args.verbose),
                    args.num_measurements, args.gamma, args.alpha
                )
                time_dict['saa'] = time.time() - t0


            with open(output_file, 'a') as f:
                for i_realization in range(args.num_realizations):

                    # apply a realization to each edge
                    realize_edge_existence(digraph, ndd_list, rs)

                    if args.use_omniscient:
                        # solve for the (omniscient) optimal edge weight
                        sol_dict['omniscient'] = optimize_picef(OptConfig(digraph,
                                                                          ndd_list,
                                                                          args.cycle_cap,
                                                                          args.chain_cap,
                                                                          verbose=args.verbose),
                                                                check_edge_success=True)
                        time_dict['omniscient'] = 0.0

                    for sol_name, (cycles, chains) in sol_dict.items():
                        score = score_after_edge_realization(cycles, chains)

                        f.write((','.join(len(output_columns) * ['%s']) + '\n') %
                                (graph_name,
                                 '%d' % i_trial,
                                 '%d' % i_realization,
                                 '%d' % args.cycle_cap,
                                 '%d' % args.chain_cap,
                                 sol_name,
                                 '%.4e' % score,
                                 time_dict[sol_name]))


def main():
    # run the experiment. sample usage:
    # >>> python robust_kex_experiment.py  --num-measurements=5  --output-dir /output/dir/ --input-dir /path/to/grpah/dir/
    parser = argparse.ArgumentParser()

    parser.add_argument('--verbose',
                        action='store_true',
                        help='verbose output')
    parser.add_argument('--seed',
                        type=int,
                        help='random seed for experiments',
                        default=0)

    # experiment params
    parser.add_argument('--num-trials',
                        type=int,
                        help='number of times to randomly assign edge distributions & measurements',
                        default=1)
    parser.add_argument('--num-realizations',
                        type=int,
                        help='number of times to realize edge failures/successes',
                        default=1)
    parser.add_argument('--input-dir',
                        type=str,
                        default=None,
                        help='input directory, containing exchange graph files (in cmu format)')
    parser.add_argument('--output-dir',
                        type=str,
                        default=None,
                        help='output directory, where an output csv will be written')

    # matching params
    parser.add_argument('--chain-cap',
                        type=int,
                        default=4,
                        help='chain cap')
    parser.add_argument('--cycle-cap',
                        type=int,
                        default=3,
                        help='cycle cap')

    # which matching methods to use
    parser.add_argument('--use-nonrobust',
                        action='store_true',
                        help='if set, calculate non-robust matching (PICEF)',
                        default=False)
    parser.add_argument('--use-nonrobust-fixed-p',
                        action='store_true',
                        help='if set, calculate non-robust matching, using fixed edge success prob. (PICEF)',
                        default=False)
    parser.add_argument('--use-nonrobust-variable-p',
                        action='store_true',
                        help='if set, calculate non-robust matching, using variable edge success prob. (PICEF)',
                        default=False)
    parser.add_argument('--use-omniscient',
                        action='store_true',
                        help='if set, calculate the omniscient optimal matching, *for each edge realization*',
                        default=False)

    # DRO-SAA specific parameters
    parser.add_argument('--use-saa',
                        action='store_true',
                        help='if set, calculate DRO-SAA matching. requires num_measurements, gamma, and alpha',
                        default=False)
    parser.add_argument('--num-measurements',
                        type=int,
                        help='number of measurements observed by DRO-SAA',
                        default=3)
    parser.add_argument('--alpha',
                        type=float,
                        help='alpha parameter, for DRO-SAA',
                        default=0.5)
    parser.add_argument('--gamma',
                        type=float,
                        help='gamma parameter, for DRO-SAA',
                        default=10.0)

    parser.add_argument('--DEBUG',
                        action='store_true',
                        help='if set, use a fixed arg string for debugging. otherwise, parse args.',
                        default=False)

    args = parser.parse_args()

    if args.DEBUG:
        # fixed set of parameters, for debugging:
        arg_str = '--num-trials 1'
        arg_str += ' --num-realizations 10'
        arg_str += ' --use-nonrobust'
        arg_str += ' --use-omniscient'
        arg_str += ' --use-saa'
        arg_str += ' --num-measurements 3'
        arg_str += ' --alpha 0.5'
        arg_str += ' --gamma 10.0'
        arg_str += ' --use-nonrobust-fixed-p'
        arg_str += ' --use-nonrobust-variable-p'
        arg_str += ' --output-dir /output/dir/'
        arg_str += ' --input-dir /input/dir/'

        args_fixed = parser.parse_args(arg_str.split())
        robust_kex_experiment(args_fixed)

    else:
        args = parser.parse_args()
        robust_kex_experiment(args)


if __name__ == "__main__":
    main()
