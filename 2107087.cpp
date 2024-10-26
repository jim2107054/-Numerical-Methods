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
    }
    cout << "No suitable root found in " << MAX_ITER << "\n";
}
