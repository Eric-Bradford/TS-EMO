# Fast and Accurate Symmetric Positive Definite Matrix Inverse Using Cholesky Decomposition

version 1.11 (2.69 KB) by Eric Blake

use LAPACK Cholesky to invert real positive definite symmetric matrix; faster more accurate than inv

We can exploit the structure of a real, positive definite, symmetric matrix by using the Cholesky decomposition to compute the inverse. 
The standard MATLAB inv function uses LU decomposition which requires twice as many operations as the Cholesky decomposition and is less
accurate. Rarely does one need to compute the inverse of a matrix (e.g. when solving a linear system, we should use \), 
but when it is needed (e.g. least squares or Kalman Filtering applications), the matrix is positive definite and symmetric. 
Thus, we can speed up and improve the accuracy of the computation by using the Cholesky decomposition. Because the Cholesky 
decomposition takes half as many operations, the algorithm will decrease run-time by about 50% for large matrices. Accuracy is 
also improved with some examples showing orders of magnitude improvement (see Felix Govaers's comment).

For the original unmodified algorithm refer to 
https://se.mathworks.com/matlabcentral/fileexchange/34511-fast-and-accurate-symmetric-positive-definite-matrix-inverse-using-cholesky-decomposition
