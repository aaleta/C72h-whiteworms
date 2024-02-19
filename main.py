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

    if args.network.startswith('CG_'):
        parameters['beta_B'] /= len(network)
        parameters['beta_W'] /= len(network)

    protection = gillespie.estimate_protection(network, parameters, initial_conditions, args.iterations)

    with open(f'results/protected_{args.network}_bB{args.beta_b}_bW{args.beta_w}_e{args.epsilon}_g{args.gamma}_m{args.mu}.pickle',
              'wb') as file:
        pickle.dump(protection, file)


def get_botnet(args):
    network = nx.read_edgelist(f'networks/{args.network}.edgelist', nodetype=int)
    parameters = {'beta_B': args.beta_b, 'beta_W': args.beta_w, 'epsilon': args.epsilon,
                  'gamma': args.gamma, 'mu': args.mu}

    initial_conditions = gillespie.random_seeds(len(network.nodes), args.n_black, args.n_white)

    if args.network.startswith('CG_'):
        parameters['beta_B'] /= len(network)
        parameters['beta_W'] /= len(network)

    botnet = gillespie.estimate_botnet(network, parameters, initial_conditions, args.iterations)

    with open(f'results/botnet_{args.network}_bB{args.beta_b}_bW{args.beta_w}_e{args.epsilon}_g{args.gamma}_m{args.mu}.pickle',
              'wb') as file:
        pickle.dump(botnet, file)


def get_botnet_threshold(args):
    network = nx.read_edgelist(f'networks/{args.network}.edgelist', nodetype=int)
    parameters = {'beta_B': args.beta_b, 'beta_W': args.beta_w, 'epsilon': args.epsilon,
                  'gamma': args.gamma, 'mu': args.mu}

    initial_conditions = gillespie.random_seeds(len(network.nodes), args.n_black, args.n_white)

    if args.network.startswith('CG_'):
        parameters['beta_B'] /= len(network)
        parameters['beta_W'] /= len(network)

    botnet_threshold = gillespie.estimate_botnet_threshold(network, parameters, args.threshold, initial_conditions, args.iterations)

    with open(f'results/botnet_threshold_{args.network}_th{args.threshold}_bB{args.beta_b}_bW{args.beta_w}_e{args.epsilon}_g{args.gamma}_m{args.mu}.pickle',
              'wb') as file:
        pickle.dump(botnet_threshold, file)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    if args.analysis_type == 'protected':
        get_protected(args)
    elif args.analysis_type == 'botnet':
        get_botnet(args)
    elif args.analysis_type == 'threshold':
        get_botnet_threshold(args)
