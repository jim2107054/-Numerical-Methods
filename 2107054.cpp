#include<bits/stdc++.h>
using namespace std;


#define TOLERANCE 1e-6  
#define MAX_ITER 1000   

vector<vector<double>> LU_A;
vector<double> LU_B;

void copyData(vector<vector<double>>&A,vector<double>&B) {
    LU_A = A;
    LU_B = B;
}


// Jacobi Iterative Method
vector<double> jacobiMethod(const vector<vector<double>>&A, const vector<double>&B) {
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

    return x_new; 
}

vector<double> gaussSeidelMethod(const vector<vector<double>>&A, const vector<double>&B) {
    int n = A.size();
    vector<double> x(n, 0); 
    
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
    
    return x; 
}

vector<double> gaussElimination(vector<vector<double>>&A, vector<double>&B) {
    int n = A.size();

    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - factor * A[k][j];
            }
            B[i] = B[i] - factor * B[k];
        }
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = B[i];
        for (int j = i + 1; j < n; j++) {
            x[i] = x[i] - A[i][j] * x[j];
        }
        x[i] = x[i] / A[i][i];
    }

    return x; 
}

// Gauss-Jordan Elimination
vector<double> gaussJordanElimination(vector<vector<double>>&A, vector<double>&B) {
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

        
        double diagElement = A[i][i];
        for (int j = 0; j <= n; j++) {
            A[i][j] /= diagElement;
        }

       
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = 0; k <= n; k++) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }
    }

    
    vector<double> solution(n);
    for (int i = 0; i < n; i++) {
        solution[i] = A[i][n]; 
    }

    return solution; 
}

bool luFactorization(vector<vector<double>> &A, vector<double> &b, vector<double>&x){
    int n = A.size();
    cout<<"n: "<<n<<endl;

    if(n == 2){
        double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        if(det == 0){
            cout << "Determinant is zero, cannot solve the equations" << endl;
            return false;
        }

        x[0] = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
        x[1] = (A[0][0] * b[1] - A[1][0] * b[0]) / det;
    }
    else if(n == 3){
        double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
        if(det == 0){
            cout << "Determinant is zero, cannot solve the equations" << endl;
            return false;
        }

        x[0] = (b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (b[1] * A[2][2] - A[1][2] * b[2]) + A[0][2] * (b[1] * A[2][1] - A[1][1] * b[2])) / det;
        x[1] = (A[0][0] * (b[1] * A[2][2] - A[1][2] * b[2]) - b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * b[2] - b[1] * A[2][0])) / det;
        x[2] = (A[0][0] * (A[1][1] * b[2] - b[1] * A[2][1]) - A[0][1] * (A[1][0] * b[2] - b[1] * A[2][0]) + b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])) / det;
    }else if(n == 4){
        double det = A[0][0] * (A[1][1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) + A[1][3] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])) - A[0][1] * (A[1][0] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])) + A[0][2] * (A[1][0] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) - A[1][1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])) - A[0][3] * (A[1][0] * (A[2][1] * A[3][2] - A[2][2] * A[3][1]) - A[1][1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0]) + A[1][2] * (A[2][0] * A[3][1] - A[2][1] * A[3][0]));
        if(det == 0){
            cout << "Determinant is zero, cannot solve the equations" << endl;
            return false;
        }

        x[0] = (b[0] * (A[1][1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) + A[1][3] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])) - A[0][1] * (b[1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (b[2] * A[3][3] - A[2][3] * b[3]) + A[1][3] * (b[2] * A[3][2] - A[2][2] * b[3])) + A[0][2] * (b[1] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) - A[1][1] * (b[2] * A[3][3] - A[2][3] * b[3]) + A[1][3] * (b[2] * A[3][1] - A[2][1] * b[3])) - A[0][3] * (b[1] * (A[2][1] * A[3][2] - A[2][2] * A[3][1]) - A[1][1] * (b[2] * A[3][2] - A[2][2] * b[3]) + A[1][2] * (b[2] * A[3][1] - A[2][1] * b[3]))) / det;
        x[1] = (A[0][0] * (b[1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (b[2] * A[3][3] - A[2][3] * b[3]) + A[1][3] * (b[2] * A[3][2] - A[2][2] * b[3])) - b[0] * (A[1][0] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])) + A[0][2] * (A[1][0] * (b[2] * A[3][3] - A[2][3] * b[3]) - b[1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * b[3] - b[2] * A[3][0])) - A[0][3] * (A[1][0] * (b[2] * A[3][2] - A[2][2] * b[3]) - A[1][2] * (A[2][0] * b[3] - b[2] * A[3][0]) + b[1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0]))) / det;
        x[2] = (A[0][0] * (A[1][1] * (b[2] * A[3][3] - A[2][3] * b[3]) - b[1] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) + A[1][3] * (A[2][1] * b[3] - b[2] * A[3][1])) - A[0][1] * (A[1][0] * (b[2] * A[3][3] - A[2][3] * b[3]) - b[1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * b[3] - b[2] * A[3][0])) + b[0] * (A[1][0] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) - A[1][1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])) - A[0][3] * (A[1][0] * (A[2][1] * b[3] - b[2] * A[3][1]) - A[1][1] * (A[2][0] * b[3] - b[2] * A[3][0]) + b[1] * (A[2][0] * A[3][1] - A[2][1] * A[3][0]))) / det;
        x[3] = (A[0][0] * (A[1][1] * (A[2][2] * b[3] - b[2] * A[3][2]) - A[1][2] * (A[2][1] * b[3] - b[2] * A[3][1]) + b[1] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])) - A[0][1] * (A[1][0] * (A[2][2] * b[3] - b[2] * A[3][2]) - A[1][2] * (A[2][0] * b[3] - b[2] * A[3][0]) + b[1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])) + A[0][2] * (A[1][0] * (A[2][1] * b[3] - b[2] * A[3][1]) - A[1][1] * (A[2][0] * b[3] - b[2] * A[3][0]) + b[1] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])) - b[0] * (A[1][0] * (A[2][1] * A[3][2] - A[2][2] * A[3][1]) - A[1][1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0]) + A[1][2] * (A[2][0] * A[3][1] - A[2][1] * A[3][0]))) / det;
    }
    return true;
};



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
                vector<double>x(A.size());
                if(luFactorization(LU_A,LU_B,x)){
                    print(x);
                }
                else{
                    cout<<"Determinant 0, can't solve."<<endl;
                }
                performOperation(A,B);
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
    copyData(A,B);
    performOperation(A,B);

    return 0;
}
