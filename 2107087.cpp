#include <bits/stdc++.h>
using namespace std;
#define MAX_ITER 1000
#define EPSILON 0.001

double evaluatePoly(const vector<double>& coeffs, double x) {
    double result = 0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * pow(x, coeffs.size() - i - 1);
    }
    return result;
}

double evaluatePolyDerivative(const vector<double>& coeffs, double x) {
    double result = 0;
    for (size_t i = 0; i < coeffs.size() - 1; ++i) {
        result += coeffs[i] * (coeffs.size() - i - 1) * pow(x, coeffs.size() - i - 2);
    }
    return result;
}

bool isRootClose(double root, const vector<double>& roots) {
    for (double r : roots) {
        if (fabs(root - r) < EPSILON * 10) {
            return true;
        }
    }
    return false;
}

// Bisection Method
void bisection(const vector<double>& coeffs, double a, double b, double tol, vector<double>& roots) {
    double mid;
    int iter = 0;

    while (fabs(b - a) > tol && iter < MAX_ITER) {
        mid = (a + b) / 2;

        if (evaluatePoly(coeffs, mid) == 0.0 || fabs(evaluatePoly(coeffs, mid)) < tol) {
            if (!isRootClose(mid, roots)) {
                roots.push_back(mid);
            }
            return;
        } else if (evaluatePoly(coeffs, mid) * evaluatePoly(coeffs, a) < 0) {
            b = mid;
        } else {
            a = mid;
        }

        iter++;
    }
}

// False Position Method
void falsePos(const vector<double>& coeffs, double a, double b, double tol, vector<double>& roots) {
    double c;
    int iter = 0;

    while (fabs(b - a) > tol && iter < MAX_ITER) {
        c = (a * evaluatePoly(coeffs, b) - b * evaluatePoly(coeffs, a)) / (evaluatePoly(coeffs, b) - evaluatePoly(coeffs, a));

        if (evaluatePoly(coeffs, c) == 0.0 || fabs(evaluatePoly(coeffs, c)) < tol) {
            if (!isRootClose(c, roots)) {
                roots.push_back(c);
            }
            return;
        } else if (evaluatePoly(coeffs, c) * evaluatePoly(coeffs, a) < 0) {
            b = c;
        } else {
            a = c;
        }

        iter++;
    }
}

// Secant Method
void secant(const vector<double>& coeffs, double x0, double x1, double tol, vector<double>& roots) {
    int iter = 0;
    double x2;

    while (iter < MAX_ITER) {
        double f0 = evaluatePoly(coeffs, x0);
        double f1 = evaluatePoly(coeffs, x1);
        if (fabs(f1 - f0) < 1e-10) return; 

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        if (fabs(evaluatePoly(coeffs, x2)) < tol && !isRootClose(x2, roots)) {
            roots.push_back(x2);
            return;
        }

        x0 = x1;
        x1 = x2;
        iter++;
    }
}

// Newton-Raphson Method
void newtonRaphson(const vector<double>& coeffs, double x0, vector<double>& roots) {
    int iterCount = 0;

    while (iterCount < MAX_ITER) {
        double fx = evaluatePoly(coeffs, x0);
        double fPrime = evaluatePolyDerivative(coeffs, x0);

        if (abs(fPrime) < EPSILON) return;

        double nextX = x0 - fx / fPrime;
        if (abs(nextX - x0) < EPSILON && !isRootClose(nextX, roots)) {
            roots.push_back(nextX);
            return;
        }
        x0 = nextX;
        iterCount++;
    }
}

// Function to find roots using selected method
void findRoots(const vector<double>& coeffs, int choice) {
    vector<double> roots;

    switch (choice) {
        case 1: // Bisection Method
            cout << "\nFinding roots using Bisection Method:\n";
            for (double i = -10; i <= 10; i += 0.5) { 
                double a = i;
                double b = a + 1; 
                bisection(coeffs, a, b, EPSILON, roots);
            }
            break;

        case 2: // False Position Method
            cout << "\nFinding roots using False Position Method:\n";
            for (double i = -10; i <= 10; i += 0.5) {
                double a = i;
                double b = a + 1; 
                falsePos(coeffs, a, b, EPSILON, roots);
            }
            break;

        case 3: // Secant Method
            cout << "\nFinding roots using Secant Method:\n";
            for (double i = -10; i <= 10; i += 0.5) {
                double x0 = i;
                double x1 = x0 + 0.1;  // Small step for Secant
                secant(coeffs, x0, x1, EPSILON, roots);
            }
            break;

        case 4: // Newton-Raphson Method
            cout << "\nFinding roots using Newton-Raphson Method:\n";
            for (double i = -10; i <= 10; i += 0.5) {
                double x0 = i;
                newtonRaphson(coeffs, x0, roots);
            }
            break;

        default:
            cout << "Invalid choice!" << endl;
            return;
    }

    // Output found roots
    if (roots.empty()) {
        cout << "No roots found.\n";
    } else {
        for (double root : roots) {
            cout << "Root: " << root << "\n";
        }
    }
}

int main() {
    vector<double> coeffs;
    int degree;

    cout << "Enter the degree of the polynomial: ";
    cin >> degree;
    coeffs.resize(degree + 1);

    cout << "Enter " << degree + 1 << " coefficients in decreasing order of powers: ";
    for (double& coef : coeffs) {
        cin >> coef;
    }

    int choice;
    do {
        cout << "\nChoose the root-finding method (0 to exit):\n";
        cout << "    [1] Bisection Method\n";
        cout << "    [2] False Position Method\n";
        cout << "    [3] Secant Method\n";
        cout << "    [4] Newton-Raphson Method\n";
        cout << "Enter your choice: ";
        cin >> choice;

        if (choice != 0) {
            findRoots(coeffs, choice);
        }
    } while (choice != 0);

    cout << "Exiting program.\n";
    return 0;
}
