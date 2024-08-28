#include <iostream>
#include <vector>
#include <cmath>    // For fabs
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Function to read parameters from input file
void read_parameters(const string& filename, double& L, double& k, double& Q, double& T_L, double& T_R, int& N, double& TOLERANCE, int& MAX_ITERATIONS) {
    ifstream inputFile(filename);
    if (!inputFile) {
        cerr << "Error opening " << filename << " file!" << endl;
        exit(1);
    }

    string line;
    while (getline(inputFile, line)) {
        stringstream ss(line);
        string key;
        ss >> key;
        
        if (key == "L") ss >> L;
        else if (key == "k") ss >> k;
        else if (key == "Q") ss >> Q;
        else if (key == "T_L") ss >> T_L;
        else if (key == "T_R") ss >> T_R;
        else if (key == "N") ss >> N;
        else if (key == "Tolerance") ss >> TOLERANCE;
        else if (key == "MaxIterations") ss >> MAX_ITERATIONS;
    }

    inputFile.close();
}

int main() {
    // Variables to hold input values
    double L, k, Q, T_L, T_R, TOLERANCE;
    int N, MAX_ITERATIONS;

    // Read parameters from input file
    read_parameters("input3.dat", L, k, Q, T_L, T_R, N, TOLERANCE, MAX_ITERATIONS);

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
