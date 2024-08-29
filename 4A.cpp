#include <iostream>
#include <iomanip>

int main() {
    const int n = 130; // Number of segments
    double dx = 1.0;   // Segment width in mm
    double T[n+1];     // Array to store temperatures at each position

    // Constants based on the calculated heat flux and integration constants
    double q = 0.069;                  // Heat flux in W/m^2
    double C1 = 600;                   // Constant for Region 1 (Plasterboard)
    double C2 = 613.56;                // Constant for Region 2 (Fiberglass)
    double k1 = 0.15;                  // Thermal conductivity for Plasterboard
    double k2 = 0.038;                 // Thermal conductivity for Fiberglass
    double k3 = 0.1;                   // Thermal conductivity for Plywood

    // Adjust C3 to ensure T(130 mm) = 400 K
    double C3 = 400 + (q / k3) * 130;  // Adjusted constant for Region 3 (Plywood)

    // Calculate temperatures for each segment
    for (int i = 0; i <= n; i++) {
        double x = i * dx; // Current position in mm
        if (x <= 10) {
            T[i] = -q/k1 * x + C1;
        } else if (x <= 110) {
            T[i] = -q/k2 * x + C2;
        } else {
            T[i] = -q/k3 * x + C3;
        }
    }

    // Print the results
    std::cout << std::fixed << std::setprecision(3);
    for (int i = 0; i <= n; i++) {
        std::cout << "Position: " << i << " mm, Temperature: " << T[i] << " K" << std::endl;
    }

    return 0;
}
