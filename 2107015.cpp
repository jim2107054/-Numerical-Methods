#include<bits/stdc++.h>
using namespace std;


//Runge Kutta method
double dydx(double x, double y){
    return (2*x+1); 
}

void rungekutta(double x0, double y0, double x, double h)
{
    int n = (int)((x - x0) / h);
    double k1, k2, k3, k4, y;
    y = y0;

    for(int i=0; i<n; i++)
     {
        k1 = h * dydx(x0, y);
        k2 = h * dydx(x0 + 0.5 * h, y + 0.5 * k1);
        k3 = h * dydx(x0 + 0.5 * h, y + 0.5 * k2);
        k4 = h * dydx(x0 + h, y + k3);
        y += (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        x0 += h;
    }
     cout<<"Solution at x = "<<x<< " is y = "<<y<< endl;

}
void printMatrix(const vector<vector<double>> &matrix){
    for (const auto &row : matrix){
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

//Inverse matrix
void matrixinverse(vector<vector<double>>&A){
    int n = A.size();
    vector<vector<double>> inv(n, vector<double>(n,0));

    
    for(int i=0; i < n; i++) inv[i][i] = 1;

    for(int i=0; i<n; i++){
        double pivot = A[i][i];
        if (pivot == 0) {
            cout<< "Matrix is singular and cannot be inverted."<<endl;
            return;
        }
        for(int j=0; j<n; j++){
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        
        for(int k=0; k<n;k++){
            if (i!=k) {
                double factor = A[k][i];
                for(int j=0; j<n; j++){
                    A[k][j] -= factor*A[i][j];
                    inv[k][j] -= factor*inv[i][j];
                }
            }
        }
    }

    cout << "Inverted Matrix:" << endl;
    printMatrix(inv);
}
