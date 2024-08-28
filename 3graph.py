import numpy as np
import matplotlib.pyplot as plt

# Parameters (replace these with actual values or read from file)
L = 0.02  # Thickness of the plate in meters
k = 0.5   # Thermal conductivity in W/m-K
Q = 1000 * 1e3  # Heat generation in W/m³ (converted from kW/m³)
T_L = 100  # Temperature at the left face in °C
T_R = 200  # Temperature at the right face in °C
N = 10    # Number of divisions
TOLERANCE = 1e-6  # Convergence tolerance
MAX_ITERATIONS = 100000

# Step size
dx = L / N

# Analytical Solution
def analytical_temperature(x):
    C2 = T_L
    C1 = (T_R - T_L + (Q / (2 * k)) * L**2) / L
    return - (Q / (2 * k)) * x**2 + C1 * x + C2

# Numerical Solution using Finite Difference Method
def numerical_solution():
    T = np.zeros(N + 1)
    T[-1] = T_R  # Boundary condition at x = L

    # Initial guess
    T[0] = T_L

    iteration = 0
    max_change = TOLERANCE + 1

    while max_change >= TOLERANCE and iteration < MAX_ITERATIONS:
        T_new = np.copy(T)
        
        # Update interior points
        for i in range(1, N):
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1]) + (dx**2 * Q) / (2 * k)
        
        # Apply boundary condition at x = 0
        T_new[0] = T_L
        
        # Update boundary condition
        T_new[-1] = T_R
        
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
plt.title('Temperature Distribution in the Plate')
plt.legend()
plt.grid(True)
plt.show()
