#include <iostream>
#include <cmath> // For std::max

using namespace std;

// Function to compute the analytical temperature distribution
void compute_analytical_temperature(double L, double k, double q_flux, double h, double T_f, double T_r, int N) {
    // Compute T0
    double T0 = (T_r + (q_flux * L / k) + h * T_f) / (1 + (h * L / k));

    // Compute the slope C1
    double C1 = - (q_flux + h * (T_f - T0)) / k;

    // Output temperature distribution
    cout << "Analytical temperature distribution:" << endl;
    for (int i = 0; i <= N; ++i) {
        double x = i * L / N;
        double T = C1 * x + T0;
        cout << "Position: " << x << " meters, Analytical Temperature: " << T << " °C" << endl;
    }
}

int main() {
    // Define parameters
    double L = 0.1;         // Thickness of the wall (meters)
    double k = 0.1;         // Thermal conductivity (W/m-K)
    double q_flux = 200;    // Heat flux at the left face (W/m^2)
    double h = 25;          // Heat transfer coefficient (W/m^2-K)
    double T_f = 300;       // Temperature of the hot fluid (°C)
    double T_r = 50;        // Temperature at the right face (°C)
    int N = 10;             // Number of divisions

    // Compute and print the analytical temperature distribution
    compute_analytical_temperature(L, k, q_flux, h, T_f, T_r, N);

    return 0;
}
