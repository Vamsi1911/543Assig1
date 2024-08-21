#include <iostream>
#include <vector>
#include <cmath>    // For fabs
using namespace std;

// Constants
const int MAX_ITERATIONS = 10000;
const double TOLERANCE = 1e-6;
const double L = 0.5;  // Length of the rod (meters)
const int N = 10;      // Number of divisions (N+1 points)
const double T0 = 100; // Temperature at the left end (x=0)
const double T1 = 500; // Temperature at the right end (x=L)

int main() {
    // Step size
    double dx = L / N;

    // Initialize temperature values
    vector<double> T(N + 1, 0);

    // Set boundary conditions
    T[0] = T0;
    T[N] = T1;

    // Initial guess: Set all interior points to an average value between T0 and T1
    for (int i = 1; i < N; ++i) {
        T[i] = (T0 + T1) / 2.0;
    }

    // Perform iterations to solve the system
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        // Copy the current temperature values
        vector<double> T_new = T;

        // Update the temperature at each interior point
        double max_change = 0.0;
        for (int i = 1; i < N; ++i) {
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1]);

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
            cout << "Converged after " << iter + 1 << " iterations." << endl;
            break;
        }
    }

    // Print the final temperature distribution
    cout << "Final temperature distribution:" << endl;
    for (int i = 0; i <= N; ++i) {
        double x = i * dx;
        cout << "Position: " << x << " meters, Temperature: " << T[i] << " C" << endl;
    }

    return 0;
}
