import numpy as np
from scipy.integrate import odeint

t, V, B, D, DB, W, WB, P = list(range(8))


def system(state, time, parameters):
    """Computes the derivative of the current state.

    :param state: contains the current fraction of V, B, D, DB, W, WB, P, in this order
    :param time: current time step
    :param parameters: should contain the parameters beta, epsilon, gamma, mu

    :returns the derivative of the states in the same order
    """
    rho_V, rho_B, rho_D, rho_DB, rho_W, rho_WB, rho_P = state

    # Auxiliary parameters
    beta, epsilon, gamma, mu = parameters['beta'], parameters['epsilon'], parameters['gamma'], parameters['mu']
    bphi = {'B': beta['B'] * (rho_B + rho_DB + rho_WB),
            'W': beta['W'] * (rho_W + rho_WB)}

    # System
    drho_V = -(bphi['B'] + bphi['W']) * rho_V
    drho_B = bphi['B'] * rho_V - bphi['W'] * rho_B
    drho_D = bphi['W'] * rho_V - (bphi['B'] + epsilon + gamma) * rho_D
    drho_DB = bphi['B'] * rho_D + bphi['W'] * rho_B - (epsilon + gamma) * rho_DB
    drho_W = epsilon * rho_D - (bphi['B'] + mu) * rho_W
    drho_WB = epsilon * rho_DB + bphi['B'] * rho_W - mu * rho_WB
    drho_P = mu * rho_WB + mu * rho_W + gamma * rho_DB + gamma * rho_D

    return [drho_V, drho_B, drho_D, drho_DB, drho_W, drho_WB, drho_P]


def solve_homogeneous_system(parameters, tmax, dt=0.01, initial_B=0.01, initial_W=0.01):
    """Solves the homogeneous system of equations.

    :param parameters: should contain the parameters beta (dict), epsilon, gamma, mu
    :param tmax: maximum time to integrate
    :param dt: timestep, defaults to 0.01
    :param initial_B: initial fraction of B agents, defaults to 0.01
    :param initial_W: initial fraction of W agents, defaults to 0.01

    :returns the solved system of equations
    """
    t = np.arange(0, tmax, dt)
    initial_state = [1.0 - initial_B - initial_W, initial_B, 0, 0, initial_W, 0, 0]

    result = odeint(system, initial_state, t, args=(parameters,))
    result = np.column_stack((t, result))

    return result


def compute_diagram(parameters, tmax, beta_B=None, beta_W=None, dt=0.01, initial_B=0.01, initial_W=0.01):
    """Computes the invasion diagram as a function of beta.

    :param parameters: should contain the parameters beta (dict), epsilon, gamma, mu
    :param tmax: maximum time to integrate
    :param beta_B: if not None, the diagram is computed over beta_B
    :param beta_W: if not None, the diagram is computed over beta_W
    :param dt: timestep, defaults to 0.01
    :param initial_B: initial fraction of B agents, defaults to 0.01
    :param initial_W: initial fraction of W agents, defaults to 0.01

    :returns the solved system of equations
    """

    if sum(param is not None for param in [beta_B, beta_W]) != 1:
        raise ValueError('The diagram can only be computed as a function of one value')

    results = np.empty((0, 2))
    if beta_W is not None:
        for beta in beta_W:
            parameters['beta']['W'] = beta
            last_value = solve_homogeneous_system(parameters, tmax, dt, initial_B, initial_W)[-1]

            results = np.vstack((results,
                                 np.array([[beta, last_value[P]]])))

    return results


def compute_2D_diagram(parameters, tmax, beta_W, epsilon_p, dt=0.01, initial_B=0.01, initial_W=0.01):
    """Computes the invasion diagram as a function of beta.

    :param parameters: should contain the parameters beta (dict), epsilon, gamma, mu
    :param tmax: maximum time to integrate
    :param beta_W: the diagram is computed over beta_W values
    :param epsilon_p: the diagram is computed over epsilon_p values
    :param dt: timestep, defaults to 0.01
    :param initial_B: initial fraction of B agents, defaults to 0.01
    :param initial_W: initial fraction of W agents, defaults to 0.01

    :returns the solved system of equations
    """

    results = np.empty((0, 3))
    for beta in beta_W:
        parameters['beta']['W'] = beta
        for epsilon in epsilon_p:
            parameters['epsilon'] = epsilon
            last_value = solve_homogeneous_system(parameters, tmax, dt, initial_B, initial_W)[-1]

            results = np.vstack((results,
                                 np.array([[beta, epsilon, last_value[P]]])))

    return results


def compute_2D_botnet(parameters, tmax, beta_B_list, beta_W_list, dt=0.01, initial_B=0.01, initial_W=0.01):
    """Computes the botnet size as a function of the betas.

    :param parameters: should contain the parameters beta (dict), epsilon, gamma, mu
    :param tmax: maximum time to integrate
    :param beta_B_list: the botnet is computed over beta_B values
    :param beta_W_list: the botnet is computed over beta_W values
    :param dt: timestep, defaults to 0.01
    :param initial_B: initial fraction of B agents, defaults to 0.01
    :param initial_W: initial fraction of W agents, defaults to 0.01

    :returns the solved system of equations
    """

    results = np.empty((0, 3))
    for beta_B in beta_B_list:
        parameters['beta']['B'] = beta_B
        for beta_W in beta_W_list:
            parameters['beta']['W'] = beta_W
            evolution = solve_homogeneous_system(parameters, tmax, dt, initial_B, initial_W)
            botnet_size = np.sum(evolution[:, [B, DB, WB]], axis=1)

            results = np.vstack((results,
                                 np.array([[beta_B, beta_W, np.max(botnet_size)]])))

    return results


def compute_botnet_threshold(parameters, tmax, force_rate_list, threshold_list, dt=0.01, initial_B=0.01, initial_W=0.01):
    """Computes the time the botnet is above the threshold as a function of epsilon/gamma.

    :param parameters: should contain the parameters beta (dict), epsilon, gamma, mu
    :param tmax: maximum time to integrate
    :param force_rate_list: the botnet is computed over force_rate values
    :param threshold_list: list of threshold size (fraction) for the botnet
    :param dt: timestep, defaults to 0.01
    :param initial_B: initial fraction of B agents, defaults to 0.01
    :param initial_W: initial fraction of W agents, defaults to 0.01

    :returns the solved system of equations
    """

    results = np.empty((0, 4))
    for threshold in threshold_list:
        for force_rate in force_rate_list:
            parameters['epsilon'] = force_rate * parameters['gamma']
            evolution = solve_homogeneous_system(parameters, tmax, dt, initial_B, initial_W)
            black_evolution = np.column_stack((evolution[:, [t]], np.sum(evolution[:, [B, DB, WB]], axis=1)))

            total_time = black_evolution[-1, 0] - black_evolution[0, 0]

            above_size = black_evolution[black_evolution[:, 1] > threshold]
            time = above_size[-1, 0] - above_size[0, 0] if len(above_size) else 0
            time = time / tmax

            results = np.vstack((results,
                                 np.array([[force_rate, threshold, time, total_time]])))

    return results
