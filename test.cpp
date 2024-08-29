#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    // Coefficients
    const double a = 1e-4;
    const double b = -0.2;

    // Range and step size for x
    const double x_start = 0.0;
    const double x_end = 0.5;
    const double x_step = 0.005;

    // Print header
    std::cout << std::setw(10) << "x (m)" << std::setw(15) << "T1 (°C)" << std::setw(15) << "T2 (°C)" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    // Iterate over x values
    for (double x = x_start; x <= x_end + 1e-8; x += x_step) {
        double c = 68.296 * x + 60.706;
        
        // Discriminant
        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            std::cout << "No real solution for x = " << x << std::endl;
            continue;
        }

        // Solutions for T
        double T1 = (-b + std::sqrt(discriminant)) / (2 * a);
        double T2 = (-b - std::sqrt(discriminant)) / (2 * a);

        // Output results to console with formatting
        std::cout << std::setw(10) << std::fixed << std::setprecision(3) << x
                  << std::setw(15) << std::fixed << std::setprecision(2) << T1
                  << std::setw(15) << std::fixed << std::setprecision(2) << T2
                  << std::endl;
    }

    return 0;
}
