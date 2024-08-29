#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

// Problem parameters
const double L = 1.0;  // Wall width (m)
const double alpha = 1.0;  // Thermal diffusivity (m²/hr)
const double T_initial = 100.0;  // Initial temperature (°C)
const double T_boundary = 300.0;  // Boundary temperature (°C)
const double t_max = 0.5;  // Maximum time (hr)

void solve_and_output(double dt, double dx, const std::string& filename) {
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
    std::cout << "dt = " << dt << ", dx = " << dx << ", Stability factor: " << stability_factor << std::endl;
    std::cout << "Stability criterion: " << (stability_factor <= 0.5 ? "Satisfied" : "Not satisfied") << std::endl;

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

    std::cout << "Data written to " << filename << std::endl << std::endl;
}

int main() {
    double dt = 0.01;  // Fixed time step

    // Test different dx values
    solve_and_output(dt, 0.10, "heat_transfer_dx_0.10.csv");  // Unstable
    solve_and_output(dt, 0.15, "heat_transfer_dx_0.15.csv");  // Marginally stable
    solve_and_output(dt, 0.20, "heat_transfer_dx_0.20.csv");  // Stable
    solve_and_output(dt, 0.25, "heat_transfer_dx_0.25.csv");  // More stable

    return 0;
}