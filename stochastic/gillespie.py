import numpy as np
import networkx as nx
from collections import defaultdict
from EoN import Gillespie_simple_contagion


def create_model(parameters):
    """Initializes the transition network of the model.

    :param parameters: dictionary containing beta_B, beta_W, epsilon, gamma, mu

    :returns a tuple containing the spontaneous and induced transitions networks
    """
    # Spontaneous transitions
    spontaneous_transitions = nx.DiGraph()
    spontaneous_transitions.add_nodes_from(['V', 'B', 'D', 'D_B', 'W', 'W_B', 'P_g', 'P_mu'])
    # White worm activation
    spontaneous_transitions.add_edge('D', 'W', rate=parameters['epsilon'])
    spontaneous_transitions.add_edge('D_B', 'W_B', rate=parameters['epsilon'])
    # User fix when prompted
    spontaneous_transitions.add_edge('D', 'P_g', rate=parameters['gamma'])
    spontaneous_transitions.add_edge('D_B', 'P_g', rate=parameters['gamma'])
    # White worm fix
    spontaneous_transitions.add_edge('W', 'P_mu', rate=parameters['mu'])
    spontaneous_transitions.add_edge('W_B', 'P_mu', rate=parameters['mu'])

    # Induced transitions
    induced_transitions = nx.DiGraph()
    induced_transitions.add_nodes_from(['V', 'B', 'D', 'D_B', 'W', 'W_B', 'P_g', 'P_mu'])
    # Black infections to V
    induced_transitions.add_edge(('B', 'V'), ('B', 'B'), rate=parameters['beta_B'])
    induced_transitions.add_edge(('D_B', 'V'), ('D_B', 'B'), rate=parameters['beta_B'])
    induced_transitions.add_edge(('W_B', 'V'), ('W_B', 'B'), rate=parameters['beta_B'])
    # White infections to V
    induced_transitions.add_edge(('W', 'V'), ('W', 'D'), rate=parameters['beta_W'])
    induced_transitions.add_edge(('W_B', 'V'), ('W_B', 'D'), rate=parameters['beta_W'])
    # White infections to B
    induced_transitions.add_edge(('W', 'B'), ('W', 'W_B'), rate=parameters['beta_W'])
    induced_transitions.add_edge(('W_B', 'B'), ('W_B', 'D_B'), rate=parameters['beta_W'])
    # Black infections to D
    induced_transitions.add_edge(('B', 'D'), ('B', 'D_B'), rate=parameters['beta_B'])
    induced_transitions.add_edge(('D_B', 'D'), ('D_B', 'D_B'), rate=parameters['beta_B'])
    induced_transitions.add_edge(('W_B', 'D'), ('W_B', 'D_B'), rate=parameters['beta_B'])
    # Black infections to W
    induced_transitions.add_edge(('B', 'W'), ('B', 'W_B'), rate=parameters['beta_B'])
    induced_transitions.add_edge(('D_B', 'W'), ('D_B', 'W_B'), rate=parameters['beta_B'])
    induced_transitions.add_edge(('W_B', 'W'), ('W_B', 'W_B'), rate=parameters['beta_B'])

    return spontaneous_transitions, induced_transitions


def stochastic_simulation(graph, parameters, initial_conditions, normalized=True):
    """Single run of the model using Gillespie's algorithm.

    :param graph: networkx graph
    :param parameters: dictionary containing beta_B, beta_W, epsilon, gamma, mu
    :param initial_conditions: dictionary containing the initial state of each node
    :param normalized: if true, the results are normalized over the number of nodes

    :return a tuple with the fraction of nodes in each state as a function of time
    """

    spontaneous_transitions, induced_transitions = create_model(parameters)
    return_statuses = ('V', 'B', 'D', 'D_B', 'W', 'W_B', 'P_g', 'P_mu')

    t, V, B, D, D_B, W, W_B, P_g, P_mu = Gillespie_simple_contagion(graph,
                                                                    spontaneous_transitions,
                                                                    induced_transitions,
                                                                    initial_conditions,
                                                                    return_statuses,
                                                                    tmax=float('Inf'))

    if normalized:
        nodes = len(graph.nodes)
        V, B, D, D_B, W, W_B, P_g, P_mu = V / nodes, B / nodes, D / nodes, D_B / nodes, W / nodes, W_B / nodes, P_g / nodes, P_mu / nodes

    return t, V, B, D, D_B, W, W_B, P_g, P_mu


def random_seeds(n_nodes, n_black, n_white):
    """Initializes the state of the nodes randomly.

    :param n_nodes: number of nodes in the network
    :param n_black: number of black worms
    :param n_white: number of white worms

    :return: a dictionary with the state of each node
    """
    initial_conditions = defaultdict(lambda: 'V')

    nodes = np.random.choice(n_nodes, n_black + n_white, replace=False)
    for node in nodes[0:n_black]:
        initial_conditions[node] = 'B'

    for node in nodes[n_black:n_black + n_white]:
        initial_conditions[node] = 'W'

    return initial_conditions


def estimate_protection(graph, parameters, initial_conditions, runs):
    """Obtain the fraction of protected nodes by the end of the simulation.

    :param graph: network on which to run the simulation
    :param parameters: dictionary containing beta_B, beta_W, epsilon, gamma and mu
    :param initial_conditions: dictionary containing the initial state of each node
    :param runs: number of iterations

    :return: the fraction of protected nodes for each run
    """

    protected = []
    for _ in range(runs):
        _, _, _, _, _, _, _, P_g, P_mu = stochastic_simulation(graph, parameters,
                                                               initial_conditions)
        protected.append(P_g[-1] + P_mu[-1])

    return protected


def estimate_botnet(graph, parameters, initial_conditions, runs):
    """Obtain the maximum fraction of nodes in the botnet.

    :param graph: network on which to run the simulation
    :param parameters: dictionary containing beta_B, beta_W, epsilon, gamma and mu
    :param initial_conditions: dictionary containing the initial state of each node
    :param runs: number of iterations

    :return: the maximum fraction of nodes in the botnet, and the final fraction of protected nodes
    """

    botnet = []
    for _ in range(runs):
        _, _, B, _, D_B, _, W_B, P_g, P_mu = stochastic_simulation(graph, parameters,
                                                                   initial_conditions)
        total_bots = B + D_B + W_B
        total_protected = P_g[-1] + P_mu[-1]
        botnet.append([np.max(total_bots), total_protected])

    return botnet


def estimate_botnet_threshold(graph, parameters, threshold, initial_conditions, runs):
    """Obtain the time evolution of the black worm.

    :param graph: network on which to run the simulation
    :param parameters: dictionary containing beta_B, beta_W, epsilon, gamma and mu
    :param threshold: fraction of nodes in the botnet for it to be dangerous
    :param initial_conditions: dictionary containing the initial state of each node
    :param runs: number of iterations

    :return: fraction of protected nodes, of the simulation above the critical size, and total simulation time
    """

    results = []
    for _ in range(runs):
        t, _, B, _, D_B, _, W_B, P_g, P_mu = stochastic_simulation(graph, parameters,
                                                                   initial_conditions)
        total_protected = P_g[-1] + P_mu[-1]

        evolution = np.column_stack((t, B + D_B + W_B))
        total_time = t[-1] - t[0]

        above_size = evolution[evolution[:, 1] > threshold]
        time = above_size[-1, 0] - above_size[0, 0] if len(above_size) else 0
        time = time / max(t)

        results.append([total_protected, time, total_time])

    return results
