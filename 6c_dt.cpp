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

// Thomas algorithm for tridiagonal matrix
std::vector<double> thomas_algorithm(const std::vector<double>& a, 
                                     const std::vector<double>& b, 
                                     const std::vector<double>& c, 
                                     const std::vector<double>& d) {
    int n = d.size();
    std::vector<double> c_star(n, 0.0);
    std::vector<double> d_star(n, 0.0);
    std::vector<double> x(n, 0.0);

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double m = 1.0 / (b[i] - a[i] * c_star[i-1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i-1]) * m;
    }

    x[n-1] = d_star[n-1];

    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_star[i] - c_star[i] * x[i+1];
    }

    return x;
}

void solve_and_output(double dt, const std::string& filename) {
    int nx = static_cast<int>(L / dx) + 1;  // Number of grid points
    int nt = static_cast<int>(t_max / dt) + 1;  // Number of time steps

    // Initialize temperature vector
    std::vector<std::vector<double> > T(nt, std::vector<double>(nx, T_initial));
    for (int n = 0; n < nt; ++n) {
        T[n][0] = T_boundary;  // Left boundary condition
        T[n][nx-1] = T_boundary;  // Right boundary condition
    }

    double r = alpha * dt / (dx * dx);
    
    std::vector<double> a(nx-2, -r);
    std::vector<double> b(nx-2, 1 + 2*r);
    std::vector<double> c(nx-2, -r);

    std::cout << "dt = " << dt << ", r = " << r << std::endl;

    // Solve using implicit method
    for (int n = 0; n < nt - 1; ++n) {
        std::vector<double> d(nx-2);
        for (int i = 1; i < nx-1; ++i) {
            d[i-1] = T[n][i];
        }
        d[0] += r * T[n+1][0];
        d[nx-3] += r * T[n+1][nx-1];

        std::vector<double> solution = thomas_algorithm(a, b, c, d);
        for (int i = 1; i < nx-1; ++i) {
            T[n+1][i] = solution[i-1];
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
    solve_and_output(0.1, "heat_transfer_implicit_dt_0.1.csv");   
    solve_and_output(0.05, "heat_transfer_implicit_dt_0.05.csv");  
    solve_and_output(0.025, "heat_transfer_implicit_dt_0.025.csv"); 
    solve_and_output(0.2, "heat_transfer_implicit_dt_0.2.csv");   
    solve_and_output(0.001, "heat_transfer_implicit_dt_0.001.csv"); 

    return 0;
}