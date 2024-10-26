<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body>

<h1>Console Application Development Using Numerical Methods</h1>

<h2>Objective</h2>
<p>This project aims to create a console application that implements various numerical methods to solve linear equations, non-linear equations, and differential equations, along with matrix inversion.</p>

<h2>Application Structure</h2>
<p>The application includes the following numerical methods:</p>

<h3>1. Solution of Linear Equations</h3>
<ul>
    <li><b>Jacobi iterative method</b></li>
    <li><b>Gauss-Seidel iterative method</b></li>
    <li><b>Gauss elimination</b></li>
    <li><b>Gauss-Jordan elimination</b></li>
    <li><b>LU factorization</b></li>
</ul>

<h3>2. Solution of Non-linear Equations</h3>
<ul>
    <li><b>Bisection method</b></li>
    <li><b>False position method</b></li>
    <li><b>Secant method</b></li>
    <li><b>Newton-Raphson method</b></li>
</ul>

<h3>3. Solution of Differential Equations</h3>
<ul>
    <li><b>Runge-Kutta method</b></li>
</ul>

<h3>4. Matrix Inversion</h3>

<p>The application is designed to solve a system of a minimum of 5 linear equations for the solution of linear equations.</p>

<h2>How to Run the Application</h2>
<ol>
    <li>Clone this repository to your local machine using:
        <pre><code>git clone https://github.com/jim2107054/-Numerical-Methods.git</code></pre>
    </li>
    <li>Navigate to the project directory.</li>
    <li>Compile and run the application using your preferred programming environment.</li>
</ol>

<h2>Solving a System of Equations</h2>

<h3>Example System of Equations</h3>
<p>Consider the following system of equations:</p>
<ol>
    <li>2x - y - 2z = -2</li>
    <li>-4x + 6y + 3z = 9</li>
    <li>-4x - 2y + 8z = -4</li>
</ol>

<h2>Solution of Linear Equations</h2>

<h3>1. Jacobi Iterative Method</h3>
<p>The <b>Jacobi Iterative Method</b> is an algorithm for determining the solutions of a diagonally dominant system of linear equations.</p>

<h4>Initialization</h4>
<ul>
    <li>Start with an initial guess for the solution vector \( x^{(0)} \).</li>
    <li>Set the tolerance level for convergence.</li>
</ul>

<h4>Iteration</h4>
<ol>
    <li>For each equation, calculate the new value for each variable:
        <p>
        [
        x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \neq i} a_{ij} x_j^{(k)} \right)
        ]
        where ( k ) is the iteration index.
        </p>
    </li>
    <li>Repeat the above step for each variable until the changes between iterations are less than the specified tolerance.</li>
</ol>

<h4>Convergence</h4>
<p>The method converges if the matrix ( A ) is strictly diagonally dominant or symmetric positive definite.</p>

<hr>

<h3>2. Gauss-Seidel Iterative Method</h3>
<p>The <b>Gauss-Seidel Method</b> is an improvement over the Jacobi method that uses the most recently updated values.</p>

<h4>Initialization</h4>
<ul>
    <li>Start with an initial guess for the solution vector ( x^{(0)} ).</li>
    <li>Set the tolerance level for convergence.</li>
</ul>

<h4>Iteration</h4>
<ol>
    <li>For each equation, update the variable immediately after calculating its new value:
        <p>
        [
        x_i^{(k+1)} = frac{1}{a_{ii}} left( b_i - sum_{j < i} a_{ij} x_j^{(k+1)} - sum_{j > i} a_{ij} x_j^{(k)} right)
        ]
        </p>
    </li>
    <li>Continue iterating until the solution converges.</li>
</ol>

<h4>Convergence</h4>
<p>Similar to the Jacobi method, convergence is guaranteed under certain conditions, such as when the matrix is strictly diagonally dominant.</p>

<hr>

<h3>3. Gauss Elimination</h3>
<p><b>Gauss Elimination</b> transforms the system of equations into an upper triangular matrix form.</p>

<h4>Steps</h4>
<ol>
    <li>Form the augmented matrix ([A | b]).</li>
    <li>Apply row operations to eliminate variables from the equations below.</li>
    <li>Transform the matrix into upper triangular form.</li>
</ol>

<h4>Back Substitution</h4>
<p>Once in upper triangular form, use back substitution to find the solution:</p>
<p>
[
x_n = frac{b_n - sum_{j=n+1}^{m} a_{nj} x_j}{a_{nn}}
]
Repeat for ( x_{n-1}, x_{n-2}, ldots, x_1 ).</p>

<hr>

<h3>4. Gauss-Jordan Elimination</h3>
<p><b>Gauss-Jordan Elimination</b> further reduces the upper triangular matrix to reduced row echelon form (RREF).</p>

<h4>Steps</h4>
<ol>
    <li>Start with the augmented matrix ([A | b]).</li>
    <li>Use row operations to achieve RREF:
        <ul>
            <li>Normalize the leading coefficient to 1.</li>
            <li>Eliminate all other entries in the leading coefficient's column.</li>
        </ul>
    </li>
</ol>

<h4>Solution</h4>
<p>The resulting matrix directly gives the solutions for each variable, as every leading variable corresponds to a column in the matrix.</p>

<hr>

<h3>5. LU Factorization</h3>
<p><b>LU Factorization</b> decomposes the matrix ( A ) into a lower triangular matrix ( L ) and an upper triangular matrix ( U ).</p>

<h4>Steps</h4>
<ol>
    <li>Decompose ( A ) such that ( A = LU ) using Gaussian elimination without row swaps.</li>
    <li>Solve ( Ly = b ) using forward substitution to find ( y ).</li>
    <li>Solve ( Ux = y ) using back substitution to find the solution vector ( x ).</li>
</ol>

<h1>Numerical Methods for Root Finding of Non-Linear Equations</h1>

<h2>1. Bisection Method</h2>
<p>The <strong>Bisection Method</strong> is a root-finding technique that systematically narrows down the interval within which a root exists. It works as follows:</p>
<ol>
    <li><strong>Initialization</strong>: Choose two points, <em>a</em> and <em>b</em>, such that <em>f(a)</em> and <em>f(b)</em> have opposite signs. This indicates that there is at least one root in the interval <em>[a, b]</em> due to the Intermediate Value Theorem.</li>
    <li><strong>Iteration</strong>:
        <ul>
            <li>Calculate the midpoint <em>c = (a + b) / 2</em>.</li>
            <li>Evaluate the function at this midpoint, <em>f(c)</em>.</li>
            <li>Determine the next interval:
                <ul>
                    <li>If <em>f(c)</em> is close to zero, then <em>c</em> is the root.</li>
                    <li>If <em>f(a)</em> and <em>f(c)</em> have opposite signs, set <em>b = c</em>.</li>
                    <li>Otherwise, set <em>a = c</em>.</li>
                </ul>
            </li>
        </ul>
    </li>
    <li><strong>Convergence</strong>: The method guarantees convergence, but it can be slow, especially if the root is located near one end of the interval.</li>
</ol>

<h2>2. False Position Method</h2>
<p>The <strong>False Position Method</strong> (or Regula Falsi) is an improvement over the Bisection Method that uses a linear interpolation approach:</p>
<ol>
    <li><strong>Initialization</strong>: Similar to the Bisection Method, start with two points <em>a</em> and <em>b</em> with <em>f(a)</em> and <em>f(b)</em> of opposite signs.</li>
    <li><strong>Iteration</strong>:
        <ul>
            <li>Instead of taking the midpoint, compute the intersection of the line connecting <em>(a, f(a))</em> and <em>(b, f(b))</em> with the x-axis to find a new point <em>c</em>.</li>
            <li>Evaluate <em>f(c)</em> and determine the next interval:
                <ul>
                    <li>If <em>f(c)</em> is close to zero, then <em>c</em> is the root.</li>
                    <li>If <em>f(a)</em> and <em>f(c)</em> have opposite signs, set <em>b = c</em>.</li>
                    <li>Otherwise, set <em>a = c</em>.</li>
                </ul>
            </li>
        </ul>
    </li>
    <li><strong>Convergence</strong>: This method often converges faster than the Bisection Method, especially when the function is approximately linear in the vicinity of the root.</li>
</ol>

<h2>3. Secant Method</h2>
<p>The <strong>Secant Method</strong> provides a way to find roots without requiring a bracketing interval:</p>
<ol>
    <li><strong>Initialization</strong>: Choose two initial approximations <em>x<sub>0</sub></em> and <em>x<sub>1</sub></em> for the root.</li>
    <li><strong>Iteration</strong>:
        <ul>
            <li>Compute the next approximation using the secant line formed by the points <em>(x<sub>0</sub>, f(x<sub>0</sub>))</em> and <em>(x<sub>1</sub>, f(x<sub>1</sub>))</em>.</li>
            <li>Update the approximations:
                <em>x<sub>n+1</sub> = x<sub>n</sub> - f(x<sub>n</sub>) * (x<sub>n</sub> - x<sub>n-1</sub>) / (f(x<sub>n</sub>) - f(x<sub>n-1</sub>))</em>.</li>
        </ul>
    </li>
    <li><strong>Convergence</strong>: The method can converge rapidly, but it may fail if the function behaves poorly or if the initial guesses are not close enough to the actual root.</li>
</ol>

<h2>4. Newton-Raphson Method</h2>
<p>The <strong>Newton-Raphson Method</strong> is a powerful technique based on the derivative of the function:</p>
<ol>
    <li><strong>Initialization</strong>: Start with an initial guess <em>x<sub>0</sub></em> for the root.</li>
    <li><strong>Iteration</strong>:
        <ul>
            <li>Compute the next approximation using:
                <em>x<sub>n+1</sub> = x<sub>n</sub> - f(x<sub>n</sub>) / f'(x<sub>n</sub>)</em>,
                where <em>f'(x<sub>n</sub>)</em> is the derivative of <em>f</em> at <em>x<sub>n</sub></em>.</li>
        </ul>
    </li>
    <li><strong>Convergence</strong>: The method generally exhibits quadratic convergence, meaning it can quickly approach the root if the initial guess is close enough and the function is well-behaved.</li>
</ol>

<h2>Method Implementation Details</h2>

<h3>1. Bisection Method</h3>
<ul>
    <li><strong>Initialization</strong>: 
        <ul>
            <li>Requires two initial guesses, <em>a</em> and <em>b</em>, such that <em>f(a)</em> and <em>f(b)</em> have opposite signs.</li>
            <li>Checks this condition at the start using an <code>if</code> statement to confirm a root exists between them.</li>
        </ul>
    </li>
    <li><strong>Iteration</strong>: 
        <ul>
            <li>Calculate the midpoint <code>mid</code> of the interval <em>[a, b]</em> using <code>mid = (a + b) / 2</code>.</li>
            <li>Evaluate the function value at the midpoint with <code>f(mid)</code>.</li>
            <li>Update <code>a</code> or <code>b</code> based on the sign of <code>f(mid)</code>:
                <ul>
                    <li>If <code>f(mid)</code> is close to zero (within tolerance), return <code>mid</code>.</li>
                    <li>If <em>f(a)</em> and <code>f(mid)</code> have opposite signs, update <code>b</code> to <code>mid</code>.</li>
                    <li>Otherwise, update <code>a</code> to <code>mid</code>.</li>
                </ul>
            </li>
        </ul>
    </li>
    <li><strong>Termination</strong>: 
        <ul>
            <li>Continue until the absolute difference between <code>a</code> and <code>b</code> is less than the specified tolerance level, ensuring convergence.</li>
        </ul>
    </li>
</ul>

<h3>2. False Position Method</h3>
<ul>
    <li><strong>Initialization</strong>: 
        <ul>
            <li>Begins with two initial points <code>a</code> and <code>b</code> that bracket a root, ensuring <code>f(a)</code> and <code>f(b)</code> have opposite signs.</li>
        </ul>
    </li>
    <li><strong>Iteration</strong>: 
        <ul>
            <li>Computes the x-intercept <code>c</code> of the line connecting <em>(a, f(a))</em> and <em>(b, f(b))</em>:
                <code>c = b - f(b) * (a - b) / (f(a) - f(b))</code>;
            </li>
            <li>Evaluate <code>f(c)</code>.</li>
        </ul>
    </li>
    <li><strong>Interval Update</strong>: 
        <ul>
            <li>Depending on the sign of <code>f(c)</code>, update <code>a</code> or <code>b</code>:
                <ul>
                    <li>If <code>f(c)</code> is close to zero, consider <code>c</code> as the root.</li>
                    <li>If <em>f(a)</em> and <code>f(c)</code> have opposite signs, update <code>b</code> to <code>c</code>.</li>
                    <li>Otherwise, update <code>a</code> to <code>c</code>.</li>
                </ul>
            </li>
        </ul>
    </li>
    <li><strong>Termination</strong>: 
        <ul>
            <li>Repeat until the difference between <code>a</code> and <code>b</code> is less than the specified tolerance level.</li>
        </ul>
    </li>
</ul>

<h3>3. Secant Method</h3>
<ul>
    <li><strong>Initialization</strong>: 
        <ul>
            <li>Starts with two initial guesses, <code>x0</code> and <code>x1</code>, expected to bracket the root.</li>
        </ul>
    </li>
    <li><strong>Iteration</strong>: 
        <ul>
            <li>Calculate the next approximation <code>x2</code>:
                <code>x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));</code>
            </li>
        </ul>
    </li>
    <li><strong>Updating Values</strong>: 
        <ul>
            <li>Update <code>x0</code> and <code>x1</code> to <code>x1</code> and <code>x2</code>, respectively.</li>
        </ul>
    </li>
    <li><strong>Termination</strong>: 
        <ul>
            <li>Continue until the absolute value of <code>f(x2)</code> is less than the specified tolerance or maximum iterations reached.</li>
        </ul>
    </li>
</ul>

<h3>4. Newton-Raphson Method</h3>
<ul>
    <li><strong>Initialization</strong>: 
        <ul>
            <li>Begins with a single initial guess <code>x0</code>.</li>
        </ul>
    </li>
    <li><strong>Iteration</strong>: 
        <ul>
            <li>Calculate the next approximation:
                <code>x1 = x0 - f(x0) / f_prime(x0);</code>
            </li>
        </ul>
    </li>
    <li><strong>Updating Values</strong>: 
        <ul>
            <li>Update <code>x0</code> to <code>x1</code>.</li>
        </ul>
    </li>
    <li><strong>Termination</strong>: 
        <ul>
            <li>Repeat until <code>f(x1)</code> is sufficiently close to zero or maximum iterations reached.</li>
        </ul>
    </li>
</ul>

<h3>5. Finding Interval for Root Detection</h3>
<ul>
    <li><strong>Scanning Range</strong>: 
        <ul>
            <li>Scans a specified range to find two consecutive points where function values have opposite signs.</li>
        </ul>
    </li>
    <li><strong>Interval Return</strong>: 
        <ul>
            <li>Returns the interval <code>[x_i, x_{i+1}]</code> for use in root-finding methods.</li>
        </ul>
    </li>
    <li><strong>Error Handling</strong>: 
        <ul>
            <li>If no interval is found, throw an exception indicating no root exists.</li>
        </ul>
    </li>
</ul>

<h3>6. Polynomial Evaluation Functions</h3>

<h4>Function: <code>evaluatePoly</code></h4>
<ul>
    <li><strong>Purpose</strong>: Evaluates a polynomial at a given <code>x</code> based on its coefficients.</li>
    <li><strong>Parameters</strong>: 
        <ul>
            <li><code>const vector<double>& coeffs</code>: Coefficients of the polynomial.</li>
            <li><code>double x</code>: Value at which the polynomial is evaluated.</li>
        </ul>
    </li>
    <li><strong>Implementation</strong>: 
        <ul>
            <li>Initialize <code>result</code> to zero.</li>
            <li>Iterate over coefficients to calculate terms:
                <pre>for (size_t i = 0; i < coeffs.size(); ++i) {
                    result += coeffs[i] * pow(x, coeffs.size() - i - 1);
                }</pre>
            </li>
        </ul>
    </li>
    <li><strong>Return Value</strong>: Returns the accumulated <code>result</code>.</li>
</ul>

<h4>Function: <code>evaluatePolyDerivative</code></h4>
<ul>
    <li><strong>Purpose</strong>: Computes the derivative of a polynomial at a given <code>x</code>.</li>
    <li><strong>Parameters</strong>: 
        <ul>
            <li><code>const vector<double>& coeffs</code>: Coefficients of the polynomial.</li>
            <li><code>double x</code>: Point at which the derivative is evaluated.</li>
        </ul>
    </li>
    <li><strong>Implementation</strong>: 
        <ul>
            <li>Initialize <code>result</code> to zero.</li>
            <li>Iterate through coefficients (excluding constant term):
                <pre>for (size_t i = 0; i < coeffs.size() - 1; ++i) {
                    result += coeffs[i] * (coeffs.size() - i - 1) * pow(x, coeffs.size() - i - 2);
                }</pre>
            </li>
        </ul>
    </li>
    <li><strong>Return Value</strong>: Returns the accumulated <code>result</code>.</li>
</ul>

<h1>Solution of Differential Equations</h1>
<h2>1. Runge-Kutta Method</h2>
<p>The <strong>Runge-Kutta Method</strong> is a numerical technique used to solve ordinary differential equations (ODEs). In this implementation, the method is specifically designed to solve first-order ODEs of the form <code>dy/dx = f(x, y)</code>. The steps involved are as follows:</p>
<ul>
    <li><strong>Derivative Function:</strong> The function <code>dydx</code> defines the derivative <code>f(x, y)</code>, which in this case is <code>2*x + 1</code>.</li>
    <li><strong>Initialization:</strong> The method starts with initial conditions <code>(x0, y0)</code> and a target <code>x</code> value, along with a step size <code>h</code> to determine how far to progress in each iteration.</li>
    <li><strong>Iteration:</strong> The method calculates the solution by performing several iterations:
        <ul>
            <li>Calculate four slopes (<code>k1</code>, <code>k2</code>, <code>k3</code>, <code>k4</code>) based on the derivative function at various points in the interval.</li>
            <li>Update the value of <code>y</code> using a weighted average of these slopes, where the weights are derived from the Runge-Kutta method formulation.</li>
            <li>Increment <code>x0</code> by the step size <code>h</code> for the next iteration.</li>
        </ul>
    </li>
    <li><strong>Output:</strong> After completing the iterations, the method outputs the estimated value of <code>y</code> at the target <code>x</code>.</li>
</ul>

<h2>2. Matrix Inversion</h2>
<p>The <strong>Matrix Inversion</strong> method computes the inverse of a square matrix using Gaussian elimination. The key steps are as follows:</p>
<ul>
    <li><strong>Initialization:</strong> Create an identity matrix of the same dimensions as the input matrix <code>A</code>, which will be transformed into the inverse.</li>
    <li><strong>Pivoting:</strong> For each row, identify the pivot element (the diagonal element). If the pivot is zero, the matrix is singular, and inversion is not possible.</li>
    <li><strong>Normalization:</strong> Normalize the current row by dividing all its elements by the pivot to make the pivot element equal to 1.</li>
    <li><strong>Elimination:</strong> Use the normalized row to eliminate corresponding elements in all other rows, effectively transforming the input matrix into the identity matrix while simultaneously applying the same operations to the identity matrix, thereby obtaining the inverse.</li>
    <li><strong>Completion:</strong> After processing all rows, the identity matrix will represent the inverse of the original matrix if the process is successful.</li>
</ul>


<h2>Video Presentation</h2>
<p>A video demonstrating the application can be viewed <a href="link_to_video">here</a>.</p>

</body>
</html>
