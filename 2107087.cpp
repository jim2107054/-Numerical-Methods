#include <bits/stdc++.h>
using namespace std;
#define MAX_ITER 1000
#define EPSILON 0.001

// Function to evaluate polynomial for a given x
double evaluatePoly(const vector<double>& coeffs, double x) {
    double result = 0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * pow(x, coeffs.size() - i - 1);
    }
    return result;
}

// Function to evaluate the derivative of the polynomial
double evaluatePolyDerivative(const vector<double>& coeffs, double x) {
    double result = 0;
    for (size_t i = 0; i < coeffs.size() - 1; ++i) {
        result += coeffs[i] * (coeffs.size() - i - 1) * pow(x, coeffs.size() - i - 2);
    }
    return result;
}

// Function to find an interval with a root by scanning a range
pair<double, double> findInterval(const vector<double>& coeffs, double range = 1000.0, double step = 0.01) {
    for (double i = -range; i < range; i += step) {
        if (evaluatePoly(coeffs, i) * evaluatePoly(coeffs, i + step) < 0) {
            return {i, i + step};
        }
    }
    throw runtime_error("No root found in the range [-" + to_string(range) + ", " + to_string(range) + "].");
}

// Secant Method
void secant(const vector<double>& coeffs, double x0, double x1, double tol) {
    cout << "\nSecant Method:\n";
    int iter = 0;
    double x2;
    do {
        double f0 = evaluatePoly(coeffs, x0);
        double f1 = evaluatePoly(coeffs, x1);
        if (fabs(f1 - f0) < 1e-10) {
            cout << "Division by nearly zero.\n";
            return;
        }
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        cout << "Iteration " << ++iter << ": x2 = " << x2 << ", f(x2) = " << evaluatePoly(coeffs, x2) << "\n";

        x0 = x1;
        x1 = x2;
    } while (fabs(evaluatePoly(coeffs, x2)) > tol && iter < MAX_ITER);

    cout << "Root (Secant): " << x2 << "\n";
}

// False Position Method
void falsePos(const vector<double>& coeffs, double a, double b, double tol) {
    cout << "\nFalse Position Method:\n";
    double c = a;
    int iter = 0;

    while (fabs(b - a) > tol && iter < MAX_ITER) {
        c = (a * evaluatePoly(coeffs, b) - b * evaluatePoly(coeffs, a)) / (evaluatePoly(coeffs, b) - evaluatePoly(coeffs, a));
        cout << "Iteration " << ++iter << ": c = " << c << ", f(c) = " << evaluatePoly(coeffs, c) << "\n";

        if (evaluatePoly(coeffs, c) == 0.0)
            break;
        else if (evaluatePoly(coeffs, c) * evaluatePoly(coeffs, a) < 0)
            b = c;
        else
            a = c;
    }
    cout << "Root (False Position): " << c << "\n";
}

// Bisection Method
void bisection(const vector<double>& coeffs, double a, double b, double tol) {
    cout << "\nBisection Method:\n";
    double mid;
    int iter = 0;

    while (fabs(b - a) > tol && iter < MAX_ITER) {
        mid = (a + b) / 2;
        cout << "Iteration " << ++iter << ": mid = " << mid << ", f(mid) = " << evaluatePoly(coeffs, mid) << "\n";

        if (evaluatePoly(coeffs, mid) == 0.0)
            break;
        else if (evaluatePoly(coeffs, mid) * evaluatePoly(coeffs, a) < 0)
            b = mid;
        else
            a = mid;
    }
    cout << "Root (Bisection): " << mid << "\n";
}

// Newton-Raphson Method
void newtonRaphson(const vector<double>& coeffs, double x0) {
    cout << "\nNewton-Raphson Method:\n";
    int iterCount = 0;

    while (iterCount < MAX_ITER) {
        double fx = evaluatePoly(coeffs, x0);
        double fPrime = evaluatePolyDerivative(coeffs, x0);

        if (abs(fPrime) < EPSILON) {
            cout << "Derivative too small.\n";
            return;
        }

        double nextX = x0 - fx / fPrime;
        cout << "Iteration " << ++iterCount << ": x = " << nextX << ", f(x) = " << evaluatePoly(coeffs, nextX) << "\n";

        if (abs(nextX - x0) < EPSILON) {
            cout << "Approximate root (Newton-Raphson): " << nextX << "\n";
            return;
        }
        x0 = nextX;
    }
    cout << "No suitable root found in " << MAX_ITER << "\n";
}

int main() {
    cout << "Enter the degree of the polynomial: ";
    int degree;
    cin >> degree;

    vector<double> coeffs(degree + 1);
    cout << "Enter " << degree + 1 << " coefficients in decreasing order of powers: ";
    for (double &coef : coeffs) {
        cin >> coef;
    }

    double tol;
    cout << "Enter tolerance level: ";
    cin >> tol;

    while (true) {
        try {
            auto [a, b] = findInterval(coeffs, 1000.0, 0.01);

            cout << "\nChoose the root-finding method:\n";
            cout << "1. Secant Method\n";
            cout << "2. False Position Method\n";
            cout << "3. Bisection Method\n";
            cout << "4. Newton-Raphson Method\n";
            cout << "5. Exit\n";
            int choice;
            cin >> choice;

            switch (choice) {
                case 1:
                    secant(coeffs, a, b, tol);
                    break;
                case 2:
                    falsePos(coeffs, a, b, tol);
                    break;
                case 3:
                    bisection(coeffs, a, b, tol);
                    break;
                case 4: {
                    double initialGuess;
                    cout << "Enter initial guess for Newton-Raphson: ";
                    cin >> initialGuess;
                    newtonRaphson(coeffs, initialGuess);
                    break;
                }
                case 5:
                    cout << "Exiting program.\n";
                    return 0;
                default:
                    cout << "Invalid choice!";
            }
        } catch (const exception& e) {
            cout << e.what() << "\n";
        }
    }
    return 0;
}
