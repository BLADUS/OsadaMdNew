import matplotlib.pyplot as plt
import numpy as np

# Define particle colors
particle1_color = 'blue'
particle2_color = 'red'

# Initialize empty lists to store particle positions
particle1_x = []
particle1_y = []
particle2_x = []
particle2_y = []

# Read data from pars_output.txt
with open('pars_output.txt', 'r') as f:
    for line in f:
        # Split line into particle identifier and coordinates
        particle_data = line.split(':')
        particle_id = particle_data[0]
        coordinates = particle_data[1].split(';')
        x = float(coordinates[0])
        y = float(coordinates[1])

        # Store coordinates based on particle identifier
        if particle_id == 'r1':
            particle1_x.append(x)
            particle1_y.append(y)
        elif particle_id == 'r2':
            particle2_x.append(x)
            particle2_y.append(y)

# Convert lists to NumPy arrays
particle1_x = np.array(particle1_x)
particle1_y = np.array(particle1_y)
particle2_x = np.array(particle2_x)
particle2_y = np.array(particle2_y)

# Create the plot
plt.figure(figsize=(10, 6))

# Plot particle 1 trajectory
plt.plot(particle1_x, particle1_y, color=particle1_color, label='Particle 1')

# Plot particle 2 trajectory
plt.plot(particle2_x, particle2_y, color=particle2_color, label='Particle 2')

# Set plot limits
plt.xlim(0, 2)
plt.ylim(0, 2)

# Set axis labels and titles
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Particle Movement')

# Add legend
plt.legend()

# Set grid
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Set x-axis ticks with 0.1 increments
plt.xticks(np.arange(0, 2.1, 0.1))

# Set y-axis ticks with 0.25 increments
plt.yticks(np.arange(0, 2.1, 0.25))

# Add markers for initial positions
plt.scatter(particle1_x[0], particle1_y[0], marker='o', color=particle1_color, label='Particle 1 Initial Position')
plt.scatter(particle2_x[0], particle2_y[0], marker='o', color=particle2_color, label='Particle 2 Initial Position')

# Add legend for initial positions
plt.legend()

# Show the plot
plt.show()
