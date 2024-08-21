#include <iostream>
#include <vector>
#include <cmath>    // For fabs

using namespace std;

// Constants
const double L = 0.02;            // Thickness of the plate (in meters)
const double k = 0.5;             // Thermal conductivity (W/m-K)
const double Q = 1000000;         // Heat generation rate (W/m^3)
const double T_L = 100;           // Temperature at the left face (°C)
const double T_R = 200;           // Temperature at the right face (°C)
const int N = 10;                // Number of divisions
const double TOLERANCE = 1e-6;    // Convergence tolerance
const int MAX_ITERATIONS = 10000; // Maximum number of iterations

int main() {
    // Step size
    double dx = L / N;

    // Initialize temperature values
    vector<double> T(N + 1, 0);

    // Set boundary conditions
    T[0] = T_L;   // Left boundary
    T[N] = T_R;   // Right boundary

    // Initial guess for interior points
    for (int i = 1; i < N; ++i) {
        T[i] = (T_L + T_R) / 2.0;  // Average of boundary temperatures
    }

    int iteration = 0;
    double max_change = 0.0;

    // Perform iterations to solve the system
    for (iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
        vector<double> T_new = T;

        max_change = 0.0; // Reset max_change for this iteration

        // Update temperatures using finite difference
        for (int i = 1; i < N; ++i) {
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1]) + (dx * dx * Q) / (2 * k);
            
            // Track the maximum change for convergence
            double change = fabs(T_new[i] - T[i]);
            if (change > max_change) {
                max_change = change;
            }
        }

        // Update the temperature values with the new values
        T = T_new;

        // Check for convergence
        if (max_change < TOLERANCE) {
            cout << "Converged after " << iteration + 1 << " iterations." << endl;
            break;
        }
    }

    if (max_change >= TOLERANCE) {
        cout << "Did not converge after " << iteration << " iterations." << endl;
    }

    // Print the final temperature distribution
    cout << "Temperature distribution:" << endl;
    for (int i = 0; i <= N; ++i) {
        double x = i * dx;
        cout << "Position: " << x << " meters, Temperature: " << T[i] << " °C" << endl;
    }

    return 0;
}
