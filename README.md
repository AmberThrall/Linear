# Linear
Linear is a C++14 linear algebra all-header library.

## Features
- Statically or dynamically sized matrices
- Multiple matrix decomposition methods, including: LUP, QR, Cholesky, Eigen, SVD and Schur
- Eigenvalue and eigenvector methods
- Equation solving of the form Ax=b
- Vectors with a distinction between row and column vectors

## Installation
Linear is an all-header library. Simply copy the src/ directory to your project's directory and include Linear.h.

## Example

```cpp
#include "Linear/Linear.h"
#include <iostream>

int main() {  
  // Get the eigenpairs of the matrix A
  Matrix3d A = { 
    {-4,14,0}, {-5,13,0}, {-1,0,2} 
  };
  std::vector<Eigenpair<double,3>> eigenpairs = Eigen(A);
  std::cout << "A = " << A << std::endl;
  for (size_t i = 0; i < eigenpairs.size(); ++i) {
      std::cout << "(" << eigenpairs[i].value << "," << Transpose(eigenpairs[i].vector) << ")" << std::endl;
  }
  
  // Compute the SVD of the matrix B
  Linear::Matrix<double,4,5> B = {
      {1,0,0,0,2}, {0,0,3,0,0}, {0,0,0,0,0}, {0,2,0,0,0}
  };
  Linear::SVD<double,4,5> svd(B);
  std::cout << "B = " << B << std::endl;
  std::cout << "USV* = " << svd.U*svd.S*svd.Vh << std::endl;
  
  // Solve the matrix equation Cx=b
  Matrix<double,3,4> C = {
      {1,2,2,2}, {2,4,6,8}, {3,6,8,10}
  };
  Vector3d b = {1,5,6};
  Vector4d x = Solve(C, b);
  std::cout << "x = " << Transpose(x) << std::endl;
}
```

## Documentation
A doxygen generated documentation is available in the doc/ directory. It is also avaiable online: http://amber.thrall.me/Linear/doc/html/index.html. The tests/ directory also contains various examples.
