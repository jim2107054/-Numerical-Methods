#include<bits/stdc++.h>
using namespace std;

// Jacobi Iterative Method
void jacobiMethod(  vector<vector<double>>& A,   vector<double>& B, int n) {
    vector<double> X(n, 0), tempX(n, 0);
    double epsilon = 1e-6;
    int maxIterations = 1000;

    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i][j] * X[j];
                }
            }
            tempX[i] = (B[i] - sum) / A[i][i];
        }

        bool converged = true;
        for (int i = 0; i < n; i++) {
            if (fabs(tempX[i] - X[i]) > epsilon) {
                converged = false;
            }
            X[i] = tempX[i];
        }

        if (converged) {
            break;
        }
    }

    cout << "Jacobi Method Solution: \n";
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << X[i] << endl;
    }
}

// Gauss-Seidel Iterative Method
void gaussSeidelMethod(  vector<vector<double>>& A,   vector<double>& B, int n) {
    vector<double> X(n, 0);
    double epsilon = 1e-6;
    int maxIterations = 1000;

    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i][j] * X[j];
                }
            }
            X[i] = (B[i] - sum) / A[i][i];
        }

        bool converged = true;
        for (int i = 0; i < n; i++) {
            if (fabs(X[i] - B[i]) > epsilon) {
                converged = false;
            }
        }

        if (converged) {
            break;
        }
    }

    cout << "Gauss-Seidel Method Solution: \n";
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << X[i] << endl;
    }
}

// Gauss Elimination
void gaussElimination(vector<vector<double>>& A, vector<double>& B, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= ratio * A[i][k];
            }
            B[j] -= ratio * B[i];
        }
    }

    vector<double> X(n);
    for (int i = n - 1; i >= 0; i--) {
        X[i] = B[i];
        for (int j = i + 1; j < n; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }

    cout << "Gauss Elimination Solution: \n";
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << X[i] << endl;
    }
}

// Gauss-Jordan Elimination
void gaussJordanElimination(vector<vector<double>>& A, vector<double>& B, int n) {
    for (int i = 0; i < n; i++) {
        double pivot = A[i][i];
        for (int j = 0; j < n; j++) {
            A[i][j] /= pivot;
        }
        B[i] /= pivot;

        for (int j = 0; j < n; j++) {
            if (i != j) {
                double factor = A[j][i];
                for (int k = 0; k < n; k++) {
                    A[j][k] -= factor * A[i][k];
                }
                B[j] -= factor * B[i];
            }
        }
    }

    cout << "Gauss-Jordan Elimination Solution: \n";
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << B[i] << endl;
    }
}

// LU Factorization
void luFactorization(vector<vector<double>>& A, vector<double>& B, int n) {
    vector<vector<double>> L(n, vector<double>(n, 0)), U(n, vector<double>(n, 0));
    vector<double> Y(n), X(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j < i)
                L[j][i] = 0;
            else {
                L[j][i] = A[j][i];
                for (int k = 0; k < i; k++) {
                    L[j][i] = L[j][i] - L[j][k] * U[k][i];
                }
            }
        }

        for (int j = 0; j < n; j++) {
            if (j < i)
                U[i][j] = 0;
            else if (j == i)
                U[i][j] = 1;
            else {
                U[i][j] = A[i][j] / L[i][i];
                for (int k = 0; k < i; k++) {
                    U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[i][j] * Y[j];
        }
        Y[i] /= L[i][i];
    }

    for (int i = n - 1; i >= 0; i--) {
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++) {
            X[i] -= U[i][j] * X[j];
        }
    }

    cout << "LU Factorization Solution: \n";
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << X[i] << endl;
    }
}

// Main function
int main() {
    int n;
    cout << "Enter the number of equations: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n);

    cout << "Enter the coefficients of the matrix A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the  ant terms of the matrix B:\n";
    for (int i = 0; i < n; i++) {
        cin >> B[i];
    }

    int choice;
    cout << "\nChoose the method to solve the system of equations:\n";
    cout << "1. Jacobi Iterative Method\n";
    cout << "2. Gauss-Seidel Iterative Method\n";
    cout << "3. Gauss Elimination Method\n";
    cout << "4. Gauss-Jordan Elimination Method\n";
    cout << "5. LU Factorization Method\n";
    cout << "Enter your choice: ";
    cin >> choice;

    switch (choice) {
        case 1:
            jacobiMethod(A, B, n);
            break;
        case 2:
            gaussSeidelMethod(A, B, n);
            break;
        case 3:
            gaussElimination(A, B, n);
            break;
        case 4:
            gaussJordanElimination(A, B, n);
            break;
        case 5:
            luFactorization(A, B, n);
            break;
        default:
            cout << "Invalid choice!";
    }

    return 0;
}
