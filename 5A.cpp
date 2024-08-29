#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    // Coefficients for the first equation
    const double a1 = 1e-4;
    const double b1 = -0.2;

    // Coefficients for the second equation
    const double a2 = 1e-4;
    const double b2 = 0.2;

    // Range and step size for x
    const double x_start = 0.0;
    const double x_end = 0.5;
    const double x_step = 0.005;

    // Print header
    std::cout << std::setw(10) << "x (m)" << std::setw(15) << "T2_1 (°C)" << std::setw(15) << "T2_2 (°C)" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    // Iterate over x values
    for (double x = x_start; x <= x_end + 1e-8; x += x_step) {
        // Discriminant and T2 for the first equation
        double c1 = 68.296 * x + 60.706;
        double discriminant1 = b1 * b1 - 4 * a1 * c1;

        double T2_1 = NAN; // Default to NaN if discriminant is negative
        if (discriminant1 >= 0) {
            T2_1 = (-b1 - std::sqrt(discriminant1)) / (2 * a1);
        }

        // Discriminant and T2 for the second equation
        double c2 = - (251.7 * x + 88.54);
        double discriminant2 = b2 * b2 - 4 * a2 * c2;



        double T2_2 = NAN; // Default to NaN if discriminant is negative
        if (discriminant2 >= 0) {
            T2_2 = (-b2 + std::sqrt(discriminant2)) / (2 * a2);
        }

        // Output results to console with formatting
        std::cout << std::setw(10) << std::fixed << std::setprecision(3) << x
                  << std::setw(15) << std::fixed << std::setprecision(2) << T2_1
                  << std::setw(15) << std::fixed << std::setprecision(2) << T2_2
                  << std::endl;
    }

    return 0;
}
