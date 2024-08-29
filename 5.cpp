#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Function to calculate thermal conductivity based on temperature
double calculate_k(double T, double beta, double k0) {
    return k0 * (1 + beta * T);
}

// Function to solve temperature distribution with temperature-dependent thermal conductivity
void solve_temperature_distribution(double beta, vector<double>& T_res, double k0, double tolerance, int max_iterations) {
    int iteration;
    double max_change;

    int N = T_res.size() - 1;

    for (iteration = 0; iteration < max_iterations; ++iteration) {
        vector<double> T_new = T_res;
        max_change = 0.0;

        for (int i = 1; i < N; ++i) {
            double k_im1_2 = calculate_k(0.5 * (T_res[i - 1] + T_res[i]), beta, k0);
            double k_ip1_2 = calculate_k(0.5 * (T_res[i + 1] + T_res[i]), beta, k0);
            T_new[i] = (k_im1_2 * T_res[i - 1] + k_ip1_2 * T_res[i + 1]) / (k_im1_2 + k_ip1_2);

            double change = fabs(T_new[i] - T_res[i]);
            if (change > max_change) max_change = change;
        }
        
        T_res = T_new;
        if (max_change < tolerance) break;
    }
}

int main() {
    // Variables to hold input values
    int N;
    double beta_positive, beta_negative, tolerance, left_temp, right_temp, rod_length, k0;
    int max_iterations;

    // Open input file
    ifstream inputFile("input5.dat");
    if (!inputFile) {
        cerr << "Error opening input.dat file!" << endl;
        return 1;
    }

    // Read the input values from the file
    string line;
    while (getline(inputFile, line)) {
        stringstream ss(line);
        string key;
        ss >> key;

        if (key == "GridSize") ss >> N;
        else if (key == "BetaPositive") ss >> beta_positive;
        else if (key == "BetaNegative") ss >> beta_negative;
        else if (key == "Tolerance") ss >> tolerance;
        else if (key == "MaxIterations") ss >> max_iterations;
        else if (key == "LeftTemp") ss >> left_temp;
        else if (key == "RightTemp") ss >> right_temp;
        else if (key == "RodLength") ss >> rod_length;
        else if (key == "k0") ss >> k0;
    }

    // Close the file
    inputFile.close();

    // Initialize temperature values
    vector<double> T(N + 1, 0);

    // Set boundary conditions
    T[0] = left_temp;   // Left boundary
    T[N] = right_temp;  // Right boundary

    // Initial guess for interior points (linear profile)
    for (int i = 1; i < N; ++i) {
        T[i] = left_temp + (right_temp - left_temp) * i / N;
    }

    // Solve for temperature distribution with temperature-dependent thermal conductivity
    vector<double> T_beta_pos = T, T_beta_neg = T;

    solve_temperature_distribution(beta_positive, T_beta_pos, k0, tolerance, max_iterations);
    solve_temperature_distribution(beta_negative, T_beta_neg, k0, tolerance, max_iterations);

    // Output results
    cout << "Temperature-Dependent Thermal Conductivity Distribution (beta = " << beta_positive << "):" << endl;
    for (int i = 0; i <= N; ++i) {
        cout << "Position: " << i * (rod_length / N) << " m, Temperature: " << T_beta_pos[i] << " C" << endl;
    }

    cout << "\nTemperature-Dependent Thermal Conductivity Distribution (beta = " << beta_negative << "):" << endl;
    for (int i = 0; i <= N; ++i) {
        cout << "Position: " << i * (rod_length / N) << " m, Temperature: " << T_beta_neg[i] << " C" << endl;
    }

    return 0;
}
