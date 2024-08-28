#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Function to read parameters from input file
void read_parameters(const string& filename, double& k1, double& k2, double& k3, double& T_L, double& T_R, int& N1, int& N2, int& N3, double& TOLERANCE, int& MAX_ITERATIONS) {
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
        
        if (key == "k1") ss >> k1;
        else if (key == "k2") ss >> k2;
        else if (key == "k3") ss >> k3;
        else if (key == "T_L") ss >> T_L;
        else if (key == "T_R") ss >> T_R;
        else if (key == "N1") ss >> N1;
        else if (key == "N2") ss >> N2;
        else if (key == "N3") ss >> N3;
        else if (key == "Tolerance") ss >> TOLERANCE;
        else if (key == "MaxIterations") ss >> MAX_ITERATIONS;
    }

    inputFile.close();
}

int main() {
    // Variables to hold input values
    double k1, k2, k3, T_L, T_R, TOLERANCE;
    int N1, N2, N3, MAX_ITERATIONS;

    // Read parameters from input file
    read_parameters("input4.dat", k1, k2, k3, T_L, T_R, N1, N2, N3, TOLERANCE, MAX_ITERATIONS);

    // Total number of divisions
    int N = N1 + N2 + N3;

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

        // Update temperatures in Plasterboard
        for (int i = 1; i < N1; ++i) {
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1]);
        }

        // Update temperatures in Fiberglass
        for (int i = N1; i < N1 + N2; ++i) {
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1]);
        }

        // Update temperatures in Plywood
        for (int i = N1 + N2; i < N; ++i) {
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1]);
        }

        // Track the maximum change for convergence
        for (int i = 1; i < N; ++i) {
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
        double x = i; // You can adjust this for actual position
        cout << "Position: " << x << " mm, Temperature: " << T[i] << " K" << endl;
    }

    return 0;
}
