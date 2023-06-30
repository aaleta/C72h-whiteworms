import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import random
import itertools


# Auxiliar functions 

# Define theta_b function 
def calculate_theta_b(rho_b, rho_cb, rho_wb, degrees, prob_degree, avg_degree):
    theta_b = 0
    for i, k in enumerate(degrees):
        theta_b += (k*prob_degree[i]/avg_degree)*(rho_b[i] + rho_cb[i] + rho_wb[i])

    return theta_b

# Define theta_w function
def calculate_theta_w(rho_w, rho_wb, degrees, prob_degree, avg_degree): 
    theta_w = 0
    for i, k in enumerate(degrees):
        theta_w += (k*prob_degree[i]/avg_degree)*(rho_w[i] + rho_wb[i])

    return theta_w


# Define the system of ODEs
def system(y, t, p, degrees, prob_degree, num_degree, avg_degree):
    # Define compartments
    V, B, C, CB, W, WB, PG, PM = range(8)

    # Calculate the densities dividing by the number of nodes with a given degree
    matrix = np.reshape(y, (8, len(degrees)))/num_degree

    # Calculate theta_b
    theta_b = calculate_theta_b(matrix[B], matrix[CB], matrix[WB], degrees, prob_degree, avg_degree)

    # Calculate theta_w
    theta_w = calculate_theta_w(matrix[W], matrix[WB], degrees, prob_degree, avg_degree)


    matrix_derivatives = np.zeros((8, len(degrees)))
    for i, k in enumerate(degrees):
        # First equations
        matrix_derivatives[V][i] = -p['beta_b'] * matrix[V][i] * k * theta_b - p['beta_w'] * matrix[V][i] * k * theta_w
        # Second equations
        matrix_derivatives[B][i] = p['beta_b'] * matrix[V][i] * k * theta_b - p['beta_w'] * matrix[B][i] * k * theta_w
        # Third equations
        matrix_derivatives[C][i] = p['beta_w'] * matrix[V][i] * k * theta_w - p['beta_b'] * matrix[C][i] * k * theta_b - p['epsilon'] * matrix[C][i] - p['gamma_p'] * matrix[C][i]
        # Fourth equations 
        matrix_derivatives[CB][i] = p['beta_b'] * matrix[C][i] * k * theta_b + p['beta_w'] * matrix[B][i] *k * theta_w - p['gamma_p'] * matrix[CB][i] - p['epsilon'] * matrix[CB][i]
        # Fifth equations
        matrix_derivatives[W][i] = p['epsilon'] * matrix[C][i] - p['beta_b'] * matrix[W][i] * k * theta_b - p['mu'] * matrix[W][i]
        # Sixth equations
        matrix_derivatives[WB][i] = p['epsilon'] * matrix[CB][i] + p['beta_b'] * matrix[W][i] * k * theta_b - p['mu'] * matrix[WB][i]
        # Seventh equations
        matrix_derivatives[PG][i] = p['gamma_p'] * matrix[C][i] + p['gamma_p'] * matrix[CB][i] 
        # Eighth equations
        matrix_derivatives[PM][i] = p['mu'] * matrix[W][i] + p['mu'] * matrix[WB][i]
    # Return array 
    return matrix_derivatives.flatten()

# Main function
def mean_field_approximation_function(G, p):
    
    # Count number of nodes
    n = len(G.nodes()) 

    # Calculate the average degree of the network
    avg_degree = np.mean([G.degree(n) for n in G.nodes()])

    # Calculate the list of unique values of degrees
    degrees = sorted(list(set([G.degree(n) for n in G.nodes()])))

    # Number of nodes with a given degree
    num_degree = [list(d for n,d in G.degree()).count(deg) for deg in degrees]

    # Probability of a node having a given degree
    prob_degree = [d/n for d in num_degree]

    # Time vector
    duration = 3000
    t = np.linspace(0, duration, duration*10)

    # Initalize nodes states 
    classes = ['V', 'B', 'C', 'W', 'C_b', 'W_b', 'PG', 'PM']
    num_nodes_per_class = [n-2, 1, 1, 0, 0, 0, 0]  # Number of nodes for each class

    # Create a list of classes, each class repeated the number of times specified
    nodes_classes = [cls for cls, num in zip(classes, num_nodes_per_class) for _ in range(num)]

    # Shuffle the list to get a random order of classes
    random.shuffle(nodes_classes)

    # Create a dictionary with nodes as keys and classes as values
    node_class_dict = {node: cls for node, cls in zip(G.nodes(), nodes_classes)}

    # Assign the classes to the nodes
    nx.set_node_attributes(G, node_class_dict, 'class')

    # Initialize a dictionary with classes as keys and empty lists as values
    class_degree_dict = {cls: [] for cls in classes}

    # Iterate over the nodes of the graph
    for node in G.nodes():
        # Get the degree and class of the node
        degree = G.degree[node]
        cls = G.nodes[node]['class']

        # Append the degree to the list corresponding to the class
        class_degree_dict[cls].append(degree)

    # Count degree frequencies for each class intially
    V0 = [class_degree_dict['V'].count(d) for d in degrees]
    B0 = [class_degree_dict['B'].count(d) for d in degrees]
    C0 = [class_degree_dict['C'].count(d) for d in degrees]
    W0 = [class_degree_dict['W'].count(d) for d in degrees]
    C_b0 = [class_degree_dict['C_b'].count(d) for d in degrees]
    W_b0 = [class_degree_dict['W_b'].count(d) for d in degrees]
    PG0 = [class_degree_dict['PG'].count(d) for d in degrees]
    PM0 = [class_degree_dict['PM'].count(d) for d in degrees]

    y0 = [V0, B0, C0, W0, C_b0, W_b0, PG0, PM0]

    y0 = list(itertools.chain(*y0))

    solution = odeint(system, y0, t, args=(p,degrees,prob_degree,num_degree,avg_degree))

    # Determine empty matrices for each class
    V = np.zeros((len(t), len(degrees)))
    B = np.zeros((len(t), len(degrees)))
    C = np.zeros((len(t), len(degrees)))
    W = np.zeros((len(t), len(degrees)))
    C_b = np.zeros((len(t), len(degrees)))
    W_b = np.zeros((len(t), len(degrees)))
    PM = np.zeros((len(t), len(degrees)))
    PG = np.zeros((len(t), len(degrees)))

    # Fill the matrices with the solution
    for i in range(len(t)): 
        V[i] = solution[i][0:len(degrees)]
        B[i] = solution[i][len(degrees):2*len(degrees)]
        C[i] = solution[i][2*len(degrees):3*len(degrees)]
        W[i] = solution[i][3*len(degrees):4*len(degrees)]
        C_b[i] = solution[i][4*len(degrees):5*len(degrees)]
        W_b[i] = solution[i][5*len(degrees):6*len(degrees)]
        PM[i] = solution[i][6*len(degrees):7*len(degrees)]
        PG[i] = solution[i][7*len(degrees):8*len(degrees)]


    # Sum over all degrees to get total numbers
    V_total = np.sum(V, axis=1)
    B_total = np.sum(B, axis=1)
    C_total = np.sum(C, axis=1)
    W_total = np.sum(W, axis=1)
    C_b_total = np.sum(C_b, axis=1)
    W_b_total = np.sum(W_b, axis=1)
    PM_total = np.sum(PM, axis=1)
    PG_total = np.sum(PG, axis=1)

    # Percentage of nodes in max_peak_black
    #max_peak_black = np.max(B_total + C_b_total + W_b_total)/n*100

    # Percentage of nodes in final_black
    final_pm = PM_total[-1]/n*100
    final_pg = PG_total[-1]/n*100

    return final_pm, final_pg

#G = nx.complete_graph(1000)

#beta_b = 1.1
#beta_w = 1.1
#gamma_p = 0.5
#epsilon = 0.5
#mu = 1

# Parameters
#p = {'beta_b': beta_b, 'beta_w': beta_w, 'gamma_p': gamma_p, 'epsilon': epsilon, 'mu': mu}

#max_peak_black, final_black = mean_field_approximation_function(G, p)

#print(max_peak_black, final_black)