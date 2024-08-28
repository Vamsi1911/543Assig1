import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 0.1  # Thickness of the wall in meters
k = 0.1  # Thermal conductivity in W/m-K
q_flux = 200  # Heat flux in W/m²
h = 25  # Heat transfer coefficient in W/m²-K
T_f = 300  # Temperature of the hot fluid in °C
T_r = 50  # Temperature at the right face in °C
N = 10  # Number of divisions
TOLERANCE = 1e-6  # Convergence tolerance

# Step size
dx = L / N

# Analytical Solution
def analytical_temperature(x):
    T0 = (T_r + (q_flux * L / k) + h * T_f) / (1 + (h * L / k))
    return (-q_flux / k + h * (T0 - T_f) / k) * x + T0

# Numerical Solution using Finite Difference Method
def numerical_solution():
    T = np.zeros(N + 1)
    T[-1] = T_r  # Boundary condition at x = L

    # Initial guess
    T[0] = T_f

    iteration = 0
    max_change = TOLERANCE + 1

    while max_change >= TOLERANCE and iteration < 100000:
        T_new = np.copy(T)
        
        # Apply boundary condition at x = 0
        heat_flux_conduction = -k * (T[1] - T[0]) / dx
        heat_flux_total = q_flux + h * (T_f - T[0])
        T_new[0] = T[0] + (heat_flux_total - heat_flux_conduction) / h
        
        # Update interior points
        for i in range(1, N):
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1])
        
        # Update boundary condition
        T_new[-1] = T_r
        
        # Check for convergence
        max_change = np.max(np.abs(T_new - T))
        T = T_new
        iteration += 1
    
    if max_change >= TOLERANCE:
        print("Did not converge after", iteration, "iterations.")
    
    return T

# Define the positions
x_values = np.linspace(0, L, N + 1)

# Compute the analytical solution
T_analytical = analytical_temperature(x_values)

# Compute the numerical solution
T_numerical = numerical_solution()

# Print the arrays
print("Position (m) | Analytical Temperature (°C) | Numerical Temperature (°C)")
print("-------------------------------------------------------------")
for x, T_a, T_n in zip(x_values, T_analytical, T_numerical):
    print(f"{x: .3f}      | {T_a: .2f}                      | {T_n: .2f}")

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_values, T_analytical, label='Analytical Solution', color='blue')
plt.plot(x_values, T_numerical, label='Numerical Solution', linestyle='--', color='red')
plt.xlabel('Position (m)')
plt.ylabel('Temperature (°C)')
plt.title('Temperature Distribution in the Wall')
plt.legend()
plt.grid(True)
plt.show()
