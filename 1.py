import numpy as np
import matplotlib.pyplot as plt

# Data for numerical and analytical solutions
positions = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])
numerical_temperatures = np.array([100, 140, 180, 220, 260, 300, 340, 380, 420, 460, 500])
analytical_temperatures = np.array([100, 140, 180, 220, 260, 300, 340, 380, 420, 460, 500])

# Plotting the temperature distributions
plt.figure(figsize=(10, 6))
plt.plot(positions, numerical_temperatures, 'o-', label='Numerical Solution', color='blue')
plt.plot(positions, analytical_temperatures, 's--', label='Analytical Solution', color='red')

plt.xlabel('Position (m)')
plt.ylabel('Temperature (°C)')
plt.title('Temperature Distribution in the Wall')
plt.legend()
plt.grid(True)

# Save the plot as a PNG file
plt.savefig('1_temperature_distribution.png')

# Display the plot
plt.show()
