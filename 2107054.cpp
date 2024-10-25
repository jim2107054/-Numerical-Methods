#include<bits/stdc++.h>
using namespace std;


#define TOLERANCE 1e-6  // Convergence tolerance
#define MAX_ITER 1000   // Maximum number of iterations

// Jacobi Iterative Method
vector<double> jacobiMethod(const vector<vector<double>>A, const vector<double>B) {
    int n = A.size();
    vector<double> x_old(n, 0);
    vector<double> x_new(n, 0);

    for (int iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += A[i][j] * x_old[j];
                }
            }
            x_new[i] = (B[i] - sum) / A[i][i];
        }

        bool converged = true;
        for (int i = 0; i < n; i++) {
            if (fabs(x_new[i] - x_old[i]) > TOLERANCE) {
                converged = false;
                break;
            }
        }

        if (converged) {
            cout << "Converged after " << iter + 1 << " iterations." << endl;
            break;
        }

        x_old = x_new;
    }

    return x_new; //Final solution
}

vector<double> gaussSeidelMethod(const vector<vector<double>>A, const vector<double>B) {
    int n = A.size();
    vector<double> x(n, 0); // Initial guess (all zeros)
    
    for (int iter = 0; iter < MAX_ITER; iter++) {
        vector<double> x_old = x; 
        
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += A[i][j] * x[j];  
                }
            }
            x[i] = (B[i] - sum) / A[i][i];
        }
        
        bool converged = true;
        for (int i = 0; i < n; i++) {
            if (fabs(x[i] - x_old[i]) > TOLERANCE) {
                converged = false;
                break;
            }
        }
        
        if (converged) {
            cout << "Converged after " << iter + 1 << " iterations." << endl;
            break;
        }
    }
    
    return x; //Final solution
}

vector<double> gaussElimination(vector<vector<double>>A, vector<double>B) {
    int n = A.size();

    // Forward Elimination Process
    for (int k = 0; k < n; k++) {
        // Partial Pivoting
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - factor * A[k][j];
            }
            B[i] = B[i] - factor * B[k];
        }
    }

    // Back Substitution
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = B[i];
        for (int j = i + 1; j < n; j++) {
            x[i] = x[i] - A[i][j] * x[j];
        }
        x[i] = x[i] / A[i][i];
    }

    return x; //Final solution
}

// Gauss-Jordan Elimination
vector<double> gaussJordanElimination(vector<vector<double>>A, vector<double>B) {
    int n = A.size();
    
    for (int i = 0; i < n; i++) {
        A[i].push_back(B[i]);
    }

    for (int i = 0; i < n; i++) {
        if (A[i][i] == 0) {
            for (int j = i + 1; j < n; j++) {
                if (A[j][i] != 0) {
                    swap(A[i], A[j]);
                    break;
                }
            }
        }

        // Normalize the current row 
        double diagElement = A[i][i];
        for (int j = 0; j <= n; j++) {
            A[i][j] /= diagElement;
        }

        // Eliminate all other entries in the current column
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = 0; k <= n; k++) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }
    }

    // Extract the solution from the augmented matrix
    vector<double> solution(n);
    for (int i = 0; i < n; i++) {
        solution[i] = A[i][n];  // The last column is the solution
    }

    return solution; //Final solution
}

// Print all the solution.
void print(vector<double>solution){
    cout << "Solution:" << endl;
    cout << fixed << setprecision(3);
    for (int i = 0; i < solution.size(); i++) {
      cout << "x" << i + 1 << " = " << solution[i] << endl;
    }
}

void performOperation(vector<vector<double>>&A, vector<double>&B){
    int choice;
    cout << "\nChoose the method to solve the system of equations:\n";
    cout << "1. Jacobi Iterative Method\n";
    cout << "2. Gauss-Seidel Iterative Method\n";
    cout << "3. Gauss Elimination Method\n";
    cout << "4. Gauss-Jordan Elimination Method\n";
    cout << "5. LU Factorization Method\n";
    cout << "6. Exit\n";
    cout << "Enter your choice: ";
    cin >> choice;

    switch (choice) {
        case 1:
            {
                vector<double> solution = jacobiMethod(A, B);

                // Print the solution
                print(solution);
                performOperation(A,B);
                break;
            }
        case 2:
            {
                vector<double> solution = gaussSeidelMethod(A, B);

                // Print the solution
                print(solution);
                performOperation(A,B);
                break;
            }
        case 3:
            {
                vector<double> solution = gaussElimination(A, B);

                // Print the solution
                print(solution);
                performOperation(A,B);
                break;
            }
        case 4:
            {
                vector<double> solution = gaussJordanElimination(A, B);

                // Print the solution
                print(solution);
                performOperation(A,B);
                break;
            }
        case 5:
            {
                //luFactorization(A, B);
                //performOperation(A,B);
                break;
            }
        case 6:
            return ;
        default:
            cout << "Invalid choice!";
            performOperation(A,B);
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
    
    performOperation(A,B);

    return 0;
}
