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

// Analytical solution
double analytical_solution(double x, double t, int terms = 100) {
    double sum = 0.0;
    for (int m = 1; m <= terms; ++m) {
        double term = std::exp(-std::pow(m * M_PI / L, 2) * alpha * t) *
                      (1 - std::pow(-1, m)) / m *
                      std::sin(m * M_PI * x / L);
        sum += term;
        
        // Check for convergence
        if (std::abs(term) < 1e-10) break;
    }
    return T_boundary + 2 * (T_initial - T_boundary) / M_PI * sum;
}

void compute_and_output(double dt, double dx, const std::string& filename) {
    int nx = static_cast<int>(L / dx) + 1;  // Number of grid points
    int nt = static_cast<int>(t_max / dt) + 1;  // Number of time steps

    std::cout << "Computing analytical solution: dt = " << dt << ", dx = " << dx << std::endl;

    // Output data to file
    std::ofstream outfile(filename);
    outfile << std::setprecision(6) << std::fixed;
    
    // Write header
    outfile << "x";
    for (int n = 0; n < nt; ++n) {
        outfile << ",t=" << n * dt;
    }
    outfile << std::endl;

    // Compute and write data
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        outfile << x;
        for (int n = 0; n < nt; ++n) {
            double t = n * dt;
            double T = analytical_solution(x, t);
            outfile << "," << T;
        }
        outfile << std::endl;
    }

    std::cout << "Data written to " << filename << std::endl << std::endl;
}

int main() {
    // Compute for different dt and dx values
    compute_and_output(0.1, 0.05, "heat_transfer_analytical_dt_0.1_dx_0.05.csv");
    compute_and_output(0.05, 0.05, "heat_transfer_analytical_dt_0.05_dx_0.05.csv");
    compute_and_output(0.01, 0.05, "heat_transfer_analytical_dt_0.01_dx_0.05.csv");
    
    compute_and_output(0.1, 0.1, "heat_transfer_analytical_dt_0.1_dx_0.1.csv");
    compute_and_output(0.1, 0.025, "heat_transfer_analytical_dt_0.1_dx_0.025.csv");

    return 0;
}