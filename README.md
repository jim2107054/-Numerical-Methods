# Console Application Development Using Numerical Methods

## Objective
This project aims to create a console application that implements various numerical methods to solve linear equations, non-linear equations, and differential equations, along with matrix inversion. 

## Application Structure
The application includes the following numerical methods:

### 1. Solution of Linear Equations
- **Jacobi iterative method**
- **Gauss-Seidel iterative method**
- **Gauss elimination**
- **Gauss-Jordan elimination**
- **LU factorization**

### 2. Solution of Non-linear Equations
- **Bisection method**
- **False position method**
- **Secant method**
- **Newton-Raphson method**

### 3. Solution of Differential Equations
- **Runge-Kutta method**

### 4. Matrix Inversion

The application is designed to solve a system of a minimum of 5 linear equations for the solution of linear equations.

## How to Run the Application
1. Clone this repository to your local machine using:
    ```bash
    git clone https://github.com/jim2107054/-Numerical-Methods.git
    ```
2. Navigate to the project directory.
3. Compile and run the application using your preferred programming environment.

## Equations: 
Example System of Equations:
 \(2x - y - 2z = -2\)
 \(-4x + 6y + 3z = 9\)
 \(-4x - 2y + 8z = -4\)

# Solving a System of Equations

## Example System of Equations

Consider the following system of equations:

1. 2x - y - 2z = -2
2. -4x + 6y + 3z = 9
3. -4x - 2y + 8z = -4

# Solution of Linear Equations

## 1. Jacobi Iterative Method

The **Jacobi Iterative Method** is an algorithm for determining the solutions of a diagonally dominant system of linear equations. 

### Initialization
- Start with an initial guess for the solution vector \( x^{(0)} \).
- Set the tolerance level for convergence.

### Iteration
1. For each equation, calculate the new value for each variable:
   \[
   x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \neq i} a_{ij} x_j^{(k)} \right)
   \]
   where \( k \) is the iteration index.
   
2. Repeat the above step for each variable until the changes between iterations are less than the specified tolerance.

### Convergence
- The method converges if the matrix \( A \) is strictly diagonally dominant or symmetric positive definite.

---

## 2. Gauss-Seidel Iterative Method

The **Gauss-Seidel Method** is an improvement over the Jacobi method that uses the most recently updated values.

### Initialization
- Start with an initial guess for the solution vector \( x^{(0)} \).
- Set the tolerance level for convergence.

### Iteration
1. For each equation, update the variable immediately after calculating its new value:
   \[
   x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)} \right)
   \]
   
2. Continue iterating until the solution converges.

### Convergence
- Similar to the Jacobi method, convergence is guaranteed under certain conditions, such as when the matrix is strictly diagonally dominant.

---

## 3. Gauss Elimination

**Gauss Elimination** transforms the system of equations into an upper triangular matrix form.

### Steps
1. Form the augmented matrix \([A | b]\).
2. Apply row operations to eliminate variables from the equations below.
3. Transform the matrix into upper triangular form.

### Back Substitution
- Once in upper triangular form, use back substitution to find the solution:
   \[
   x_n = \frac{b_n - \sum_{j=n+1}^{m} a_{nj} x_j}{a_{nn}}
   \]
   Repeat for \( x_{n-1}, x_{n-2}, \ldots, x_1 \).

---

## 4. Gauss-Jordan Elimination

**Gauss-Jordan Elimination** further reduces the upper triangular matrix to reduced row echelon form (RREF).

### Steps
1. Start with the augmented matrix \([A | b]\).
2. Use row operations to achieve RREF:
   - Normalize the leading coefficient to 1.
   - Eliminate all other entries in the leading coefficient's column.

### Solution
- The resulting matrix directly gives the solutions for each variable, as every leading variable corresponds to a column in the matrix.

---

## 5. LU Factorization

**LU Factorization** decomposes the matrix \( A \) into a lower triangular matrix \( L \) and an upper triangular matrix \( U \).

### Steps
1. Decompose \( A \) such that \( A = LU \) using Gaussian elimination without row swaps.
2. Solve \( Ly = b \) using forward substitution to find \( y \).
3. Solve \( Ux = y \) using back substitution to find the solution vector \( x \).

---




# Numerical Methods for Root Finding of Non-Linear Equations

## 1. Bisection Method
The **Bisection Method** is a root-finding technique that systematically narrows down the interval within which a root exists. It works as follows:

1. **Initialization**: Choose two points, \( a \) and \( b \), such that \( f(a) \) and \( f(b) \) have opposite signs. This indicates that there is at least one root in the interval \([a, b]\) due to the Intermediate Value Theorem.

2. **Iteration**:
   - Calculate the midpoint \( c = \frac{a + b}{2} \).
   - Evaluate the function at this midpoint, \( f(c) \).
   - Determine the next interval:
     - If \( f(c) \) is close to zero, then \( c \) is the root.
     - If \( f(a) \) and \( f(c) \) have opposite signs, set \( b = c \).
     - Otherwise, set \( a = c \).
   - Repeat this process until the interval is sufficiently small, indicating that \( c \) is an accurate approximation of the root.

3. **Convergence**: The method guarantees convergence, but it can be slow, especially if the root is located near one end of the interval.

## 2. False Position Method
The **False Position Method** (or Regula Falsi) is an improvement over the Bisection Method that uses a linear interpolation approach:

1. **Initialization**: Similar to the Bisection Method, start with two points \( a \) and \( b \) with \( f(a) \) and \( f(b) \) of opposite signs.

2. **Iteration**:
   - Instead of taking the midpoint, compute the intersection of the line connecting \( (a, f(a)) \) and \( (b, f(b)) \) with the x-axis to find a new point \( c \).
   - Evaluate \( f(c) \) and determine the next interval:
     - If \( f(c) \) is close to zero, then \( c \) is the root.
     - If \( f(a) \) and \( f(c) \) have opposite signs, set \( b = c \).
     - Otherwise, set \( a = c \).
   - Repeat until the desired accuracy is achieved.

3. **Convergence**: This method often converges faster than the Bisection Method, especially when the function is approximately linear in the vicinity of the root.

## 3. Secant Method
The **Secant Method** provides a way to find roots without requiring a bracketing interval:

1. **Initialization**: Choose two initial approximations \( x_0 \) and \( x_1 \) for the root.

2. **Iteration**:
   - Compute the next approximation using the secant line formed by the points \( (x_0, f(x_0)) \) and \( (x_1, f(x_1)) \).
   - Update the approximations:
     - \( x_{n+1} = x_n - f(x_n) \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})} \).
   - Repeat this process until convergence criteria are met (e.g., when the change in \( x \) is below a certain threshold).

3. **Convergence**: The method can converge rapidly, but it may fail if the function behaves poorly or if the initial guesses are not close enough to the actual root.

## 4. Newton-Raphson Method
The **Newton-Raphson Method** is a powerful technique based on the derivative of the function:

1. **Initialization**: Start with an initial guess \( x_0 \) for the root.

2. **Iteration**:
   - Compute the next approximation using:
     - \( x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)} \),
     where \( f'(x_n) \) is the derivative of \( f \) at \( x_n \).
   - Repeat the iteration until the change in \( x \) is sufficiently small or until \( f(x_n) \) is close to zero.

3. **Convergence**: The method generally exhibits quadratic convergence, meaning it can quickly approach the root if the initial guess is close enough and the function is well-behaved.

## Method Implementation Details

### 1. Bisection Method
- **Initialization**: 
  - Requires two initial guesses, \a\ and \b\, such that \f(\a\) and \f(\b \) have opposite signs.
  - Checks this condition at the start using an \if\ statement to confirm a root exists between them.
  
- **Iteration**: 
  - Calculate the midpoint \`mid\` of the interval \[\a\, \b\] using \`mid = (\a\ + \b\) / 2\`.
  - Evaluate the function value at the midpoint with \`f(mid)\`.
  - Update \`a\` or \`b\` based on the sign of \`f(mid)\`:
    - If \`f(mid)\` is close to zero (within tolerance), return \`mid\`.
    - If \`f(\a\) and \f(mid)\` have opposite signs, update \`b\` to \`mid\`.
    - Otherwise, update \`a\` to \`mid\`.

- **Termination**: 
  - Continue until the absolute difference between \`a\` and \`b\` is less than the specified tolerance level, ensuring convergence.

### 2. False Position Method
- **Initialization**: 
  - Begins with two initial points \`a\` and \`b\` that bracket a root, ensuring \`f(\a\)` and \`f(\b\)` have opposite signs.

- **Iteration**: 
  - Computes the x-intercept \`c\` of the line connecting \`(\a\, f(\a\))\` and \`(\b\, f(\b\))\`:
    
    \c \= \b\ - f(\b\) * (\a\ - \b\) / (f(\a\) - f(\b\));
   
  - Evaluate \f(c)\.

- **Interval Update**: 
  - Depending on the sign of \`f(c)\`, update \`a\` or \`b\`:
    - If \`f(c)\` is close to zero, consider \`c\` as the root.
    - If \`f(\a\) and \`f(c)\ have opposite signs, update \b\ to \c\.
    - Otherwise, update \`a\` to \`c\`.

- **Termination**: 
  - Repeat until the difference between \`a\` and \`b\` is less than the specified tolerance level.

### 3. Secant Method
- **Initialization**: 
  - Starts with two initial guesses, \`x0\` and \`x1\`, expected to bracket the root.

- **Iteration**: 
  - Calculate the next approximation \`x2\`:
    
   \ x2\ = \x1\ - f(\x1\) * (\x1\ - \x0\) / (f(\x1\) - f(\x0\));

- **Updating Values**: 
  - Update \`x0\` and \`x1\` to \`x1\` and \`x2\`, respectively.

- **Termination**: 
  - Continue until the absolute value of \`f(x2)\` is less than the specified tolerance or maximum iterations reached.

### 4. Newton-Raphson Method
- **Initialization**: 
  - Begins with a single initial guess \`x0\`.

- **Iteration**: 
  - Calculate the next approximation:
   
    \x1\ = \x0\ - f(\x0\) / f_prime(\x0\);
  

- **Updating Values**: 
  - Update \`x0\` to \`x1\`.

- **Termination**: 
  - Repeat until \`f(x1)\` is sufficiently close to zero or maximum iterations reached.

### 5. Finding Interval for Root Detection
- **Scanning Range**: 
  - Scans a specified range to find two consecutive points where function values have opposite signs.

- **Interval Return**: 
  - Returns the interval \`[\x_i\, \x_{i+1}\]` for use in root-finding methods.

- **Error Handling**: 
  - If no interval is found, throw an exception indicating no root exists.

### 6. Polynomial Evaluation Functions

#### Function: \`evaluatePoly\`
- **Purpose**: Evaluates a polynomial at a given \`x\` based on its coefficients.
- **Parameters**: 
  - \`const vector<double>& coeffs\`: Coefficients of the polynomial.
  - \`double x\`: Value at which the polynomial is evaluated.

- **Implementation**: 
  - Initialize \`result\` to zero.
  - Iterate over coefficients to calculate terms:
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * pow(x, coeffs.size() - i - 1);
    }

- **Return Value**: Returns the accumulated \`result\`.

#### Function: \`evaluatePolyDerivative\`
- **Purpose**: Computes the derivative of a polynomial at a given \`x\`.
- **Parameters**: 
  - \`const vector<double>& coeffs\`: Coefficients of the polynomial.
  - \`double x\`: Point at which the derivative is evaluated.

- **Implementation**: 
  - Initialize \`result\` to zero.
  - Iterate through coefficients (excluding constant term):
    for (size_t i = 0; i < coeffs.size() - 1; ++i) {
        result += coeffs[i] * (coeffs.size() - i - 1) * pow(x, coeffs.size() - i - 2);
    }

- **Return Value**: Returns the accumulated \`result\`.

## 5. Runge-Kutta Method
The **Runge-Kutta Method** is commonly used for solving ordinary differential equations numerically:

1. **Initialization**: Start with an initial condition \( (t_0, y_0) \) where \( y_0 \) is the solution at time \( t_0 \).

2. **Iteration**:
   - Calculate several intermediate values to improve accuracy:
     - Compute \( k_1, k_2, k_3, k_4 \) based on the derivative function \( f(t, y) \).
   - Update the value of \( y \):
     - Use a weighted average of the \( k \) values to estimate the new \( y \) value at the next time step \( t_1 = t_0 + h \).

3. **Convergence**: The fourth-order Runge-Kutta method (RK4) is particularly popular due to its balance between computational efficiency and accuracy.

## Video Presentation
A video demonstrating the application can be viewed [here](link_to_video).
