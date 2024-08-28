#include <iostream>
#include <iomanip>  // For setting output precision

using namespace std;

int main() {
    // Constants
    double k1 = 0.15;  // W/m-K
    double k2 = 0.038; // W/m-K
    double k3 = 0.1;   // W/m-K
    double T_L = 600;  // K
    double T_R = 400;  // K

    // Lengths of each material (in meters)
    double L1 = 10.0e-3; // 10 mm
    double L2 = 100.0e-3; // 100 mm
    double L3 = 20.0e-3; // 20 mm

    double L = L1 + L2 + L3; // Total length

    // Discretization
    int N1 = 10; // Number of points in region 1
    int N2 = 100; // Number of points in region 2
    int N3 = 20; // Number of points in region 3
    int N = N1 + N2 + N3; // Total number of points

    // Calculate constants
    double x1 = L1 / L; // End of region 1
    double x2 = (L1 + L2) / L; // End of region 2

    // Compute temperature gradients and offsets
    double A1 = (T_R - T_L) / (k1 * L1 + k2 * L2 + k3 * L3);
    double B1 = T_L - A1 * 0; // Temperature at x=0

    double A2 = (T_R - (A1 * x1 * k1 * L1 + B1)) / (k2 * L2);
    double B2 = (T_L - A1 * 0 * k1 * L1 - B1) / k2; // Temperature at x=x1

    double A3 = (T_R - (A1 * x1 * k1 * L1 + B1 + A2 * (x2 - x1) * k2 + B2)) / (k3 * L3);
    double B3 = T_R - A3 * (x2 - x1) * k3; // Temperature at x=x2

    // Print the temperature distribution
    cout << "Position (m) | Temperature (K)" << endl;
    cout << "-----------------------------" << endl;
    for (int i = 0; i <= N; ++i) {
        double x = i * (L / N);

        double T;
        if (x <= x1) {
            T = A1 * x + B1;
        } else if (x <= x2) {
            T = A2 * (x - x1) + B2;
        } else {
            T = A3 * (x - x2) + B3;
        }

        cout << fixed << setprecision(3) << x << "       | " << fixed << setprecision(2) << T << endl;
    }

    return 0;
}
