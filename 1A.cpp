#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Function to read parameters from the input file
void read_parameters(const string& filename, double& L, int& N, double& T0, double& T1) {
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
    }

    inputFile.close();
}

// Function to calculate the analytical solution
double analytical_solution(double x, double T0, double T1, double L) {
    return T0 + (T1 - T0) * (x / L);
}

int main() {
    // Variables to hold parameters
    double L, T0, T1;
    int N;

    // Read parameters from input file
    read_parameters("input1.dat", L, N, T0, T1);

    // Step size
    double dx = L / N;

    // Print the analytical temperature distribution
    cout << "Analytical temperature distribution:" << endl;
    for (int i = 0; i <= N; ++i) {
        double x = i * dx;
        double T_analytical = analytical_solution(x, T0, T1, L);
        cout << "Position: " << x << " meters, Analytical Temperature: " << T_analytical << " °C" << endl;
    }

    return 0;
}
