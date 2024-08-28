#include <iostream>
#include <iomanip>  // For setting output precision
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Function to read parameters from input file
void read_parameters(const string& filename, double& L, double& k, double& Q, double& T_L, double& T_R, int& N) {
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
        else if (key == "Q") ss >> Q;
        else if (key == "T_L") ss >> T_L;
        else if (key == "T_R") ss >> T_R;
        else if (key == "N") ss >> N;
    }

    inputFile.close();
}

int main() {
    // Variables to hold input values
    double L, k, Q, T_L, T_R;
    int N;

    // Read parameters from input file
    read_parameters("input3.dat", L, k, Q, T_L, T_R, N);

    // Calculate constants for the analytical solution
    double C2 = T_L;
    double C1 = (T_R - T_L + (Q / (2 * k)) * L * L) / L;
    
    // Calculate the temperature distribution
    cout << "Position (m) | Temperature (°C)" << endl;
    cout << "------------------------------" << endl;
    for (int i = 0; i <= N; ++i) {
        double x = i * L / N;
        double T = - (Q / (2 * k)) * x * x + C1 * x + C2;
        cout << fixed << setprecision(3) << x << "       | " << fixed << setprecision(2) << T << endl;
    }

    return 0;
}
