# Getting Started

Matrices in Linear are comprised of the following things:
1. A type format to use, i.e., float, double, etc.
2. Two non-negative integers M and N representing the number of rows and number of columns respectively.
3. Whether the matrix is row-major or column-major for storage.
4. An array of Complex numbers.

Construction of a row-major matrix A then follows by the code
```cpp
    Matrix<double,3,3> A = {
        {1,2,3}, {4,5,6}, {7,8,9}
    };
```
This construction is a bit cumbersome, so the shortcuts Matrix{M}x{N}{T} have been added, where M=2,3,4; N=2,3,4 and T=f,d,ld for float, double, and long double respectively. For square matrices, the shortcuts Matrix{N}{T} have also been added alongside with SquareMatrix<typename T, size_t N, unsigned int Flags=0>. So we can concise our code as follows:
```cpp
    Matrix3d A = {
        {1,2,3}, {4,5,6}, {7,8,9}
    };
```
Note these two code snippets are identical.

## Operators

Matrices in Linear have all the operators one would except.
```cpp
    Matrix3d A = {
        {1,2,3}, {4,5,6}, {7,8,9}
    };
    Matrix3d B(1.0);
    std::cout << "A+B = " << A+B << std::endl;
    std::cout << "A-B = " << A-B << std::endl;
    std::cout << "AB = " << A*B << std::endl;
    std::cout << "BA = " << B*A << std::endl;
    std::cout << "A*2 = " << A*2.0 << std::endl;
    std::cout << "A/2 = " << A/2.0 << std::endl;
    std::cout << "A==B = " << A==B << std::endl;
    std::cout << "A!=B = " << A!=B << std::endl;
```
Note here that B is the 3x3 row-major matrix with every entry set to 1.

## Accessing

There are multiple ways to determine the size of a matrix A.
```cpp
    Matrix3d A = {
        {1,2,3}, {4,5,6}, {7,8,9}
    };
    std::cout << A.NumRows() << std::endl;
    std::cout << A.NumColumns() << std::endl;
    std::cout << A.Size() << std::endl;
    std::cout << A(1,0) << std::endl;
```
The above code prints out
```
3
3
(1x2)[3, 3]
4
```
A few things to point out here. First, indices for matrices start at 0 for both row and column. Second, the order for the accessor A(1,0)
is always (row,column) even if the matrix is column-major. Finally, A.Size() returns a 1x2 Matrix containing the dimensions of A.

## Dynamic Sizing

So far we have been working with statically sized matrices. Statically sized matrices have a set size, given by their template parameters,
and never change. Unfortunately, statically sized matrices must have the size determined at compile time. If we need to determine their size at runtime we can make use of dynamically sized matrices.
```cpp
    Matrix<double,Dynamic,Dynamic> A = {
        {1,2,3}, {4,5,6}, {7,8,9}
    };
```
The keyword 'Dynamic' replacing the previous 3's given as sizes indicate to the matrix that is dynamically sized. Construction of dynamically sized matrices has several methods. The first is with initializer lists, as seen above. A in the above code is a 3x3 matrix whose size was determined dynamically. A few other constructors are
```cpp
    MatrixXd A(3, 1.0);
    MatrixXd B(3, 3, 1.0);
```
Both A and B here are dynamically sized 3x3 matrices with every entry set to one (note MatrixXd is a shortcut for Matrix<double,Dynamic,Dynamic>).
One benefit of dynamically sized matrices, is that their size can change. There are two ways to do this, the assignment operator '=' and the Resize function.
```cpp
    Matrix4d A(1.0);
    MatrixXd B(3, 1.0);
    B = A;
    B.Resize(3, 3);
```
If a matrix resize is attempted on a statically sized matrix, an exception is thrown.

Finally, you can combine static and dynamic sizing.
```cpp
    Matrix<double,3,Dynamic> A(3, 1.0);
```
In the above code we create a 3x3 matrix where the number of rows is static at 3, and the number of columns is dynamically determined to be 3. The constructor A(*,3,1.0) would have also worked. Just the first parameter would be ignored.

## Vectors

Column and row vectors in Linear are just special cases of matrices. The shortcut Vector<typename T, size_t N> is identical to Matrix<T,N,1>. Likewise the shortcut RowVector<typename T, size_t N> is identical to Matrix<T,1,N>. Like before, there are a few more shortcuts. Namely Vector{N}{T} and RowVector{N}{T} where N=2,3,4 and T=f,d,ld. The code below is an example of a vector construction of length 3.
```cpp
    Vector3d x = {1,2,3};
    std::cout << x.Length() << std::endl;
    std::cout << x[1] << std::endl;
```
The above code prints out
```
3
2
```
Notice that a new size accessor pops up for vectors, Length. This does exactly what it sounds like, it gets the length of the vector.

While vectors can be accessed by x(i,0), it is suggested to use the square brackets x[i]. That way one doesn't have to worry about if x is a row or column vector.
