import networkx as nx
from EoN import Gillespie_simple_contagion as gsp_alg
from collections import defaultdict
import pandas as pd
from tqdm import trange
import numpy as np
import argparse
import sys


def parse_args(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()

    # Simulation parameters
    parser.add_argument(
        "-n",
        "--netname",
        type=str,
        default='ER',
        help="name of the graph to be used"
    )

    parser.add_argument(
        '-b',
        '--beta',
        type=float,
        default=1.1,
        help='beta value'
    )

    parser.add_argument(
        "-g",
        "--gammap",
        type=float,
        help='prompted transition rate'
    )

    parser.add_argument(
        "-e",
        "--epsilon",
        type=float,
        help='ethics rate'
    )

    parser.add_argument(
        "-i",
        "--iterations",
        type=int,
        default=100,
        help='number of iterations'
    )

    return parser.parse_args(args)


def run_stoc_sim(params, G, G_name, IC, n_iter):
    #initiliaze spontaneous transition network

    nodes = ['V', 'C', 'W', 'B', 'C_W', 'W_B', 'P']
    edges = [
        ('C_B', 'P', {'rate': params['gammaP']}), ('C', 'P', {'rate':params['gammaP']}), #user update prompted
        ('W_B', 'P', {'rate': params['mu']}), ('W', 'P', {'rate': params['mu']}), #user update prompted or the worm updates 
        ('C', 'W', {'rate':params['epsilon']}), ('C_B', 'W_B', {'rate':params['epsilon']}) # white worm gets active
            ]

    sp_proc = nx.DiGraph()
    sp_proc.add_nodes_from(nodes)
    sp_proc.add_edges_from(edges)

    # initialize induced transition network

    nodes = [
        ('B', 'V'), ('B', 'B'), ('C_B', 'V'), ('C_B', 'B'), ('W_B', 'V'), ('W_B', 'B'),  # Transition from vulnerable to black infected
        ('W', 'V'), ('W', 'C'), ('W_B', 'V'), ('W_B', 'C'), # Transition from vulnerable to "cautious"
        ('W', 'B'), ('W', 'C_B'), ('W_B', 'B'), ('W_B', 'C_B'), # Transition from black infected to "cautious"+black infected
        ('B', 'C'), ('B', 'C_B'), ('C_B', 'C'), ('C_B', 'C_B'), ('W_B', 'C'), ('W_B', 'C_B'), # Transition from cautious to cautious+black infected
        ('B', 'W'), ('B', 'W_B'), ('C_B', 'W'), ('C_B', 'W_B'), ('W_B', 'W'), ('W_B', 'W_B') # Transition from White infected to double infected
            ]

    edges = [
        (nodes[0], nodes[1], {'rate':params['beta']}), 
        (nodes[2], nodes[3], {'rate':params['beta']}),
        (nodes[4], nodes[5], {'rate':params['beta']}),
        (nodes[6], nodes[7], {'rate':params['beta']}),
        (nodes[8], nodes[9], {'rate':params['beta']}),
        (nodes[10], nodes[11], {'rate':params['beta']}),
        (nodes[12], nodes[13], {'rate':params['beta']}),
        (nodes[14], nodes[15], {'rate':params['beta']}), 
        (nodes[16], nodes[17], {'rate':params['beta']}),
        (nodes[18], nodes[19], {'rate':params['beta']}),
        (nodes[20], nodes[21], {'rate':params['beta']}), 
        (nodes[22], nodes[23], {'rate':params['beta']}),
        (nodes[24], nodes[25], {'rate':params['beta']})
             ]

    ind_proc = nx.DiGraph()
    ind_proc.add_nodes_from(nodes)
    ind_proc.add_edges_from(edges)


    # initial conditions
    return_statuses = ('V', 'B', 'W_B', 'C_B', 'C', 'W', 'P')
    
    max_peack_black = np.zeros(n_iter)
    final_black = np.zeros(n_iter)
    t_window = np.zeros(n_iter)
    
    N = len(G.nodes)
    beta=params['beta']
    gamma = params['gammaP']
    epsilon = params['epsilon']
    
    for k in trange(n_iter):
        
        t, V, B, W_B, C_B, C, W, P  = gsp_alg(G, sp_proc, ind_proc, IC, return_statuses,
                                            tmax = float('Inf'))
        
        max_peack_black[k] = np.max(B+W_B+C_B)
        final_black[k] = B[-1]
        
        ts_indx = np.where((B+W_B+C_B)/N>0.5)[0]
        if len(ts_indx)>0:
            t_window[k] = t[ts_indx[-1]]-t[ts_indx[0]]
        else:
            t_window[k] = 0
        
        if k==0:
                returns = (V, B, W_B, C_B, C, W, P)
                dict_returns = {key : returns[i] for i, key in enumerate(return_statuses)}
                dict_returns['t'] = t
                pd.DataFrame.from_dict(dict_returns, orient='index').T.to_csv(f'results/sim_G{G_name}_beta{beta}_gammaP{gamma}_epsilon{epsilon}.csv', index=False)

    dict_results = {}
    dict_results['max_peack_black'] = max_peack_black
    dict_results['final_black'] = final_black
    dict_results['t_window'] = t_window
    pd.DataFrame.from_dict(dict_results, orient='index').to_csv(f'results/avg_sim_G{G_name}_beta{beta}_gammaP{gamma}_epsilon{epsilon}.csv')
    
    return max_peack_black, final_black, t_window


if __name__ == '__main__':
    args = parse_args()

    G_name = args.netname
    G = nx.read_edgelist(
        f'networks/network_{G_name}_10000.txt',
        create_using=nx.Graph(),
        delimiter=', ',
        nodetype=int)

    # Initial conditions
    IC = defaultdict(lambda: 'V')
    for node in range(1):
        IC[node] = 'B'
    for node in range(2, 3):
        IC[node] = 'C'

    params = {'beta': 1.1,
              'gammaP': args.gammap,
              'epsilon': args.epsilon,
              'mu': 1}

    run_stoc_sim(params, G, G_name, IC, args.iterations)
