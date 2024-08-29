import numpy as np
import matplotlib.pyplot as plt

# Data for numerical and analytical solutions
positions = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
numerical_temperatures = np.array([298.077, 273.269, 248.462, 223.654, 198.846, 174.038, 149.231, 124.423, 99.6154, 74.8077, 50])
analytical_temperatures = np.array([298.077, 273.269, 248.462, 223.654, 198.846, 174.038, 149.231, 124.423, 99.6154, 74.8077, 50])

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
plt.savefig('2_temperature_distribution.png')

# Display the plot
plt.show()
