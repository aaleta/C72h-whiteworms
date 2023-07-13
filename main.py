import sys
import pickle
import networkx as nx
from parser import parse_args
import stochastic.gillespie as gillespie


def get_protected(args):
    network = nx.read_edgelist(f'networks/{args.network}.edgelist', nodetype=int)
    parameters = {'beta_B': args.beta_b, 'beta_W': args.beta_w, 'epsilon': args.epsilon,
                  'gamma': args.gamma, 'mu': args.mu}

    initial_conditions = gillespie.random_seeds(len(network.nodes), args.n_black, args.n_white)

    protection = gillespie.estimate_protection(network, parameters, initial_conditions, args.iterations)

    with open(f'results/protected_{args.network}_bB{args.beta_b}_bW{args.beta_w}_e{args.epsilon}_g{args.gamma}_m{args.mu}.pickle', 'wb') as file:
        pickle.dump(protection, file)



if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    if args.analysis_type == 'protected':
        get_protected(args)
