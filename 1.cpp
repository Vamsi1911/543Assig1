#include <iostream>
#include <vector>
#include <cmath>    // For fabs
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Function to read parameters from the input file
void read_parameters(const string& filename, double& L, int& N, double& T0, double& T1, double& TOLERANCE, int& MAX_ITERATIONS) {
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
        else if (key == "N") ss >> N;
        else if (key == "T0") ss >> T0;
        else if (key == "T1") ss >> T1;
        else if (key == "TOLERANCE") ss >> TOLERANCE;
        else if (key == "MAX_ITERATIONS") ss >> MAX_ITERATIONS;
    }

    inputFile.close();
}

int main() {
    // Variables to hold input values
    double L, T0, T1, TOLERANCE;
    int N, MAX_ITERATIONS;

    // Read parameters from input file
    read_parameters("input1.dat", L, N, T0, T1, TOLERANCE, MAX_ITERATIONS);

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
