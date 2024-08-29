#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// Problem parameters
const double L = 1.0;  // Wall width (m)
const double alpha = 1.0;  // Thermal diffusivity (m²/hr)
const double T_initial = 100.0;  // Initial temperature (°C)
const double T_boundary = 300.0;  // Boundary temperature (°C)
const double t_max = 0.5;  // Maximum time (hr)
const double dx = 0.05;  // Grid size (m)

void solve_and_output(double dt, const std::string& filename) {
    int nx = static_cast<int>(L / dx) + 1;  // Number of grid points
    int nt = static_cast<int>(t_max / dt) + 1;  // Number of time steps

    // Initialize temperature vector
    std::vector<std::vector<double> > T(nt, std::vector<double>(nx, T_initial));
    for (int n = 0; n < nt; ++n) {
        T[n][0] = T_boundary;  // Left boundary condition
        T[n][nx-1] = T_boundary;  // Right boundary condition
    }

    // Stability criterion
    double stability_factor = alpha * dt / (dx * dx);
    std::cout << "dt = " << dt << ", Stability factor: " << stability_factor << std::endl;
    if (stability_factor > 0.5) {
        std::cout << "Warning: The solution may be unstable!" << std::endl;
    }

    // Solve using explicit method
    for (int n = 1; n < nt; ++n) {
        for (int i = 1; i < nx - 1; ++i) {
            T[n][i] = T[n-1][i] + stability_factor * (T[n-1][i+1] - 2*T[n-1][i] + T[n-1][i-1]);
        }
    }

    // Output data to file
    std::ofstream outfile(filename);
    outfile << std::setprecision(6) << std::fixed;
    
    // Write header
    outfile << "x";
    for (int n = 0; n < nt; n += std::max(1, nt / 6)) {
        outfile << ",t=" << n * dt;
    }
    outfile << std::endl;

    // Write data
    for (int i = 0; i < nx; ++i) {
        outfile << i * dx;
        for (int n = 0; n < nt; n += std::max(1, nt / 6)) {
            outfile << "," << T[n][i];
        }
        outfile << std::endl;
    }

    std::cout << "Data written to " << filename << std::endl;
}

int main() {
    // Test different time steps
    solve_and_output(0.1, "heat_transfer_dt_0.1.csv");   // Original time step
    solve_and_output(0.05, "heat_transfer_dt_0.05.csv");  // Smaller time step
    solve_and_output(0.025, "heat_transfer_dt_0.025.csv"); // Even smaller time step
    solve_and_output(0.2, "heat_transfer_dt_0.2.csv");   // Larger time step (potentially unstable)
    solve_and_output(0.001, "heat_transfer_dt_0.001.csv"); // Very small time step (should be stable)

    return 0;
}