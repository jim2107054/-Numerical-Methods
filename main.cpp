#include<bits/stdc++.h>
#include"2107015.cpp"
#include"2107054.cpp"
#include"2107087.cpp"
using namespace std;


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

// JIM_2107054 function
void JIM_2107054() {
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
}

//ZISAN_2107015 functions.
void ZISAN_2107015(){
    int choice;
    do{
        cout << "Numerical Methods Application\n";
        cout << "1. Runge-Kutta Method\n";
        cout << "2. Matrix Inversion (n x n matrix)\n";
        cout << "0. Exit\n";
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice){
            case 1:{
                //equation is: 2*x+1
                double x0 = 0, y0 = 1, x = 2, h = 0.2;
                rungekutta(x0, y0, x, h);
                break;
            }
            case 2:{
                int n;
                cout << "Enter the dimension of the matrix (n x n): ";
                cin >> n;
                vector<vector<double>> A(n, vector<double>(n));

                cout << "Enter the elements of the matrix:\n";
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        cin >> A[i][j];
                    }
                }

                matrixinverse(A);
                break;
            }
            case 0:
                cout << "Exiting..." << endl;
                break;
            default:
                cout << "Invalid choice!" << endl;
        }
    }while(choice != 0);

}

// AFIFA_2107087 functions.
void AFIFA_2107087() {
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
                    return;
                default:
                    cout << "Invalid choice!";
            }
        } catch (const exception& e) {
            cout << e.what() << "\n";
        }
    }
}

int main()
{
    int option;
    cout << "\nChoose Application Structure:\n";
    cout << "1. Solution of Linear Equations \n";
    cout << "2. Solution of Non-linear Equations \n";
    cout << "3. Solution of Differential Equations or Matrix Inversion \n";
    cout << "4. Exit\n";
    cout << "Enter your choice: ";
    cin >>option;
;

    switch (option) {
        case 1:
            {
                JIM_2107054();
                main();
                break;
            }
        case 2:
            {
                AFIFA_2107087();
                main();
                break;
            }
        case 3:
            {
                ZISAN_2107015();
                main();
                break;
            }
        case 4:
            {
                return 0;
                break;
            }
        default:
            cout << "Invalid choice!";
            break;
    }
    return 0;
}