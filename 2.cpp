#include <iostream>
#include <vector>
#include <cmath>    // For fabs

using namespace std;

// Constants
const double L = 0.1;            // Thickness of the wall (in meters)
const double k = 0.1;            // Thermal conductivity (W/m-K)
const double q_flux = 200;       // Heat flux at left face (W/m^2)
const double h = 25;             // Heat transfer coefficient (W/m^2-K)
const double T_f = 300;          // Temperature of the hot fluid (°C)
const double T_r = 50;           // Temperature at the right face (°C)
const int N = 20;                // Number of divisions in the wall
const double TOLERANCE = 1e-6;   // Convergence tolerance

int main() {
    // Step size
    double dx = L / N;

    // Initialize temperature values
    vector<double> T(N + 1, 0);

    // Set boundary condition at the right face (x = L)
    T[N] = T_r;

    // Initial guess for the temperature at the left face (x = 0)
    double T0 = T_f; // Initial guess close to the hot fluid temperature

    int iteration = 0;
    double max_change = 0.0;

    // Perform iterations to solve the system
    for (iteration = 0; iteration < 100000; ++iteration) {
        // Copy the current temperature values
        vector<double> T_new = T;

        // Apply boundary condition at the left face (x = 0)
        double heat_flux_conduction = -k * (T[1] - T0) / dx;
        double heat_flux_total = q_flux + h * (T_f - T0);  
        T_new[0] = T0 + (heat_flux_total - heat_flux_conduction) / h;

        // Update the temperature at each interior point using finite difference
        max_change = 0.0; // Reset max_change for this iteration
        for (int i = 1; i < N; ++i) {
            double old_temp = T[i];
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1]);

            // Track the maximum change for convergence
            double change = fabs(T_new[i] - old_temp);
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
    cout << "Final temperature distribution:" << endl;
    for (int i = 0; i <= N; ++i) {
        double x = i * dx;
        cout << "Position: " << x << " meters, Temperature: " << T[i] << " °C" << endl;
    }

    return 0;
}
