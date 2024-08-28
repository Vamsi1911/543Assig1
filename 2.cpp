#include <iostream>
#include <vector>
#include <cmath>    // For fabs
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Function to read parameters from input file
void read_parameters(const string& filename, double& L, double& k, double& q_flux, double& h, double& T_f, double& T_r, int& N, double& TOLERANCE) {
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
        else if (key == "q_flux") ss >> q_flux;
        else if (key == "h") ss >> h;
        else if (key == "T_f") ss >> T_f;
        else if (key == "T_r") ss >> T_r;
        else if (key == "N") ss >> N;
        else if (key == "Tolerance") ss >> TOLERANCE;
    }

    inputFile.close();
}

int main() {
    // Variables to hold input values
    double L, k, q_flux, h, T_f, T_r, TOLERANCE;
    int N;

    // Read parameters from input file
    read_parameters("input2.dat", L, k, q_flux, h, T_f, T_r, N, TOLERANCE);

    // Step size
    double dx = L / N;

    // Initialize temperature values
    vector<double> T(N + 1, 0);

    // Set boundary condition at the right face (x = L)
    T[N] = T_r;

    // Perform iterations to solve the system
    int iteration = 0;
    double max_change = TOLERANCE + 1;

    while (max_change >= TOLERANCE && iteration < 100000) {
        // Copy the current temperature values
        vector<double> T_new = T;

        // Apply boundary condition at the left face (x = 0)
        double heat_flux_conduction = -k * (T[1] - T[0]) / dx;
        double heat_flux_total = q_flux + h * (T_f - T[0]);
        T_new[0] = T[0] + (heat_flux_total - heat_flux_conduction) / h;

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
        iteration++;
    }

    if (max_change >= TOLERANCE) {
        cout << "Did not converge after " << iteration << " iterations." << endl;
    } else {
        cout << "Converged after " << iteration << " iterations." << endl;
    }

    // Print the final temperature distribution
    cout << "Final temperature distribution:" << endl;
    for (int i = 0; i <= N; ++i) {
        double x = i * dx;
        cout << "Position: " << x << " meters, Temperature: " << T[i] << " °C" << endl;
    }

    return 0;
}
