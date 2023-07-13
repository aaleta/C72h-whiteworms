import math
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def plot_evolution(data, dt=0.01):
    """Plot the time evolution of all states.

    :param data: a numpy array of size (timesteps, 7)
    :param dt: size of the integration step
    """

    timesteps = np.arange(data.shape[0]) * dt

    # Define the colors and labels for each state
    colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black']
    labels = ['V', 'B', 'D', 'D_B', 'W', 'W_B', 'P']

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Iterate over each state and plot its time evolution
    for i in range(data.shape[1]):
        ax.plot(timesteps, data[:, i], color=colors[i], label=labels[i])

    # Set the x and y labels
    ax.set_xlabel('Time')
    ax.set_ylabel('State Value')

    # Set the title
    ax.set_title('Time Evolution of States')

    # Add a legend
    ax.legend()

    # Display the plot
    plt.show()


def plot_diagram(data_list, label_x='beta_W', label_y='P'):
    """Plot the invasion diagram."""

    # Create a figure and axis
    fig, ax = plt.subplots()

    for data_dict in data_list:
        legend_text = data_dict['legend']
        data = data_dict['data']
        threshold = data_dict['threshold']

        x, y = data[:, 0], data[:, 1]
        ax.plot(x, y, label=legend_text)

        ax.axvline(x=threshold, ymin=0, ymax=0.2, color='black', linestyle='--', linewidth=2)

    # Set the y limits
    ax.set_ylim(bottom=0)

    # Set the x and y labels
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)

    # Add a legend
    ax.legend()

    # Display the plot
    plt.show()


def plot_2D_diagram(data, mu=1, gamma=1):
    """Plots the diagram as a heatmap."""

    # Extract beta, epsilon, and P values from the data array
    beta_values = data[:, 0]
    epsilon_values = data[:, 1]
    P_values = data[:, 2]

    # Determine the grid size
    grid_size = int(np.sqrt(len(P_values)))

    # Reshape the P values to create a 2D grid
    P_grid = P_values.reshape((grid_size, grid_size))

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Create a heatmap using imshow
    heatmap = ax.imshow(P_grid, extent=[np.min(epsilon_values), np.max(epsilon_values),
                                        np.min(beta_values), np.max(beta_values)],
                        origin='lower', aspect='auto', cmap='viridis')

    # Calculate the theoretical line values
    epsilon_unique = np.unique(epsilon_values)
    beta_theo = mu * ((epsilon_unique/gamma) + 1) / (epsilon_unique/gamma)

    # Plot the theoretical line
    ax.plot(epsilon_unique, beta_theo, 'r--', linewidth=2, label='Theoretical Line')

    # Set the x and y labels
    ax.set_xlabel('epsilon')
    ax.set_ylabel('beta')

    # Add a colorbar
    cbar = fig.colorbar(heatmap)

    # Add limits
    ax.set_xlim(np.min(epsilon_values), np.max(epsilon_values))
    ax.set_ylim(np.min(beta_values), np.max(beta_values))

    # Add a legend
    ax.legend()

    # Display the plot
    plt.show()


def plot_3D_protected(data):
    """Plot the final fraction of protected nodes in 3D, compared to the theoretical one."""

    def R0(beta, epsilon, gamma=1, mu=1):
        ep = epsilon / gamma
        return (beta / mu) * (ep / (ep + 1))

    def func(x, R0):
        return x - 1 + math.exp(-R0 * x)

    # Extract the columns from the data array
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    # Create a 3D figure and axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the data points as a scatter plot
    ax.scatter(x, y, z, label='numerical')

    unique_values = np.unique(x)

    # Iterate over unique values and plot lines
    for value in unique_values:
        indices = np.where(x == value)
        z_theo = [optimize.root(func, [0.5], args=(R0(value, epsilon),)).x[0]
                  for epsilon in y[indices]]
        ax.plot(x[indices], y[indices], z_theo, color='black')

    # Set the x, y, and z labels
    ax.set_xlabel('beta')
    ax.set_ylabel('epsilon')
    ax.set_zlabel('protected')

    # Chaging the view
    ax.view_init(elev=30, azim=-120)

    # Create a proxy artist for the legend
    legend_lines = [Line2D([0], [0], color='black', linewidth=2)]
    legend_labels = ['analytical']

    # Add a legend
    ax.legend(legend_lines, legend_labels)

    # Display the plot
    plt.show()


def plot_2D_botnet(data, epsilon=1, gamma=1, mu=1):
    """Plots the botnet size as a heatmap."""

    # Extract beta, epsilon, and P values from the data array
    beta_B_values = data[:, 0]
    beta_W_values = data[:, 1]
    botnet_size = data[:, 2]

    # Determine the grid size
    grid_size = int(np.sqrt(len(botnet_size)))

    # Reshape the P values to create a 2D grid
    P_grid = botnet_size.reshape((grid_size, grid_size))

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Create a heatmap using imshow
    heatmap = ax.imshow(P_grid, extent=[np.min(beta_W_values), np.max(beta_W_values),
                                        np.min(beta_B_values), np.max(beta_B_values)],
                        origin='lower', aspect='auto', cmap='viridis',
                        vmin=0, vmax=1)

    # Calculate the theoretical line values
    beta_W_unique = np.unique(beta_W_values)
    beta_theo = 0.5*(-(epsilon + mu + gamma) + np.sqrt((epsilon + gamma - mu)**2 + 4 * beta_W_unique * epsilon))

    # Plot the theoretical line
    ax.plot(beta_W_unique, beta_theo, 'r--', linewidth=2, label='Theoretical Line')

    # Set the x and y labels
    ax.set_xlabel('beta_W')
    ax.set_ylabel('beta_B')

    # Add a colorbar
    cbar = fig.colorbar(heatmap)

    # Set the title
    ax.set_title(f'epsilon={epsilon} gamma={gamma}')

    # Add limits
    ax.set_xlim(np.min(beta_W_values), np.max(beta_W_values))
    ax.set_ylim(np.min(beta_B_values), np.max(beta_B_values))

    # Add a legend
    ax.legend()

    # Display the plot
    plt.show()


def plot_botnet_threshold(beta_W_max, epsilon, gamma, mu):
    """Plot the botnet threshold."""

    values_list = epsilon if len(epsilon) > 1 else gamma

    # Create a color gradient based on the values
    colors = plt.cm.viridis(np.linspace(np.min(values_list), np.max(values_list), len(values_list)))

    # Calculate beta_B values for each value
    beta_W_values = np.linspace(0, beta_W_max, 100)
    for value, color in zip(values_list, colors):
        ep = value if len(epsilon) > 1 else epsilon
        ga = value if len(gamma) > 1 else gamma

        beta_B_values = 0.5 * (-(ep + mu + ga) +
                               np.sqrt((ep + ga - mu) ** 2 + 4 * beta_W_values * ep))
        mask = beta_B_values >= 0
        plt.plot(beta_W_values[mask], beta_B_values[mask], color=color)

    # Set x-axis and y-axis labels
    plt.xlabel('beta_W')
    plt.ylabel('beta_B')

    # Create a color bar for epsilon values
    label = 'epsilon' if len(epsilon) > 1 else 'gamma'

    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    sm.set_array(values_list)
    cbar = plt.colorbar(sm)
    cbar.set_label(label)

    # Add a legend
    plt.legend()

    # Show the plot
    plt.show()
