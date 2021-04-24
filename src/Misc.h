#pragma once
#include "Matrix.h"
#include "Basics.h"
#include "Vector.h"
#include "Types.h"
#include "Decomp.h"
#include "Global.h"

namespace Linear {
    /**
     * Computes a basis for the null space of A.
     * The null space of a matrix is the set of vectors v such that Av=0.
     * @param A MxN matrix
     * @return List of vectors v such that \f$Null(A)=span\{v[0],\dots,v[len(v)-1]\}\f$
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    std::vector<Vector<T,N>> NullSpace(const Matrix<T,M,N,Flags>& A) {
        std::vector<Vector<T,N>> basis;

        SquareMatrix<T,N,Flags> eye = Identity<T>(A.NumColumns());
        Matrix<T,(M==Dynamic||N==Dynamic?Dynamic:M+N),N,Flags> aug = RowAugmented(A, eye);
        aug = Transpose(RREF(Transpose(aug)));
        SquareMatrix<T,N,Flags> c = SubMatrix(aug, A.NumColumns(), A.NumColumns(), A.NumRows(), 0);
        for (size_t col = 0; col < A.NumColumns(); ++col) {
            bool add = true;
            for (size_t row = 0; row < A.NumRows(); ++row) {
                if (aug(row, col) != T(0)) {
                    add = false;
                    break;
                }
            }
            if (add)
                basis.push_back(c.GetColumn(col));
        }
        return basis;
    }
    /**
     * Computes the dimension of the A's null space.
     * @param A MxN matrix
     * @return dim(Null(A))
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    unsigned int Nullity(const Matrix<T,M,N,Flags>& A) {
        return NullSpace(A).size();
    }

    /**
     * Computes a basis for the column space of A.
     * The column space of a matrix is the set of vectors v such that Ax=v for some vector x.
     * @param A MxN matrix
     * @return List of vectors v such that \f$colsp(A)=span\{v[0],\dots,v[len(v)-1]\}\f$
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    std::vector<Vector<T,N>> ColumnSpace(const Matrix<T,M,N,Flags>& A) {
        std::vector<Vector<T,N>> basis;

        Matrix<T,M,N,Flags> B = RREF(A);
        size_t pivot = 0;
        for (size_t c = 0; c < A.NumColumns(); ++c) {
            bool add = true;
            if (B(pivot, c) == T(1)) {
                basis.push_back(A.GetColumn(c));
                pivot += 1;
            }
        }
        return basis;
    }
    /**
     * Computes the dimension of the A's column space.
     * @param A MxN matrix
     * @return dim(colsp(A))
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    unsigned int Rank(const Matrix<T,M,N,Flags>& A) {
        return ColumnSpace(A).size();
    }

    /**
     * Solves the matrix equation Ax=b for x.
     * This functions breaks into various cases.
     * 1. If A is square we have 2 special cases.
     *   a. If A is invertible, then \f$x=A^{-1}b\f$.
     *   b. If A is upper triangular or lower triangular, then we can use forward/backward subsitution to solve.
     * 2. If those two cases don't hold, we employ Guassian elimination.
     * If M != P, an exception is thrown. If there is no solution, an exception is thrown.
     * @param A MxN matrix
     * @param b Vector of length P
     * @return Vector of length N
     */
    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P>
    typename std::enable_if<(M==P||P==Dynamic||M==Dynamic), Vector<T,N>>::type Solve(const Matrix<T,M,N,Flags>& A, const Vector<T,P>& b) {
        if (A.NumRows() != b.Length())
            throw "Cannot solve matrix equation Ax=b, b is incorrect size.";

        Vector<T,N> x(A.NumColumns(), T(0));

        if (IsSquare(A)) { // Handle square matrices case.
            if (Abs(Determinant(A)) > T(Tol)) { // Simply invert, x = A^{-1}b
                return Inverse(A)*b;
            }
            else if (IsUpperTriangular(A)) { // Upper triangular, back substition.
                for (int i = A.NumColumns()-1; i >= 0; i--) {
                    // A(i,i)x(i) + A(i,i+1)x(i+1) + ... + A(i,N-1)x(N-1) = b(i)
                    Complex<T> sum;
                    for (size_t j = i+1; j < A.NumColumns(); ++j)
                        sum += A(i,j)*x[j];
                    if (Abs(A(i,i)) < T(Tol)) {
                        x[i] = T(0);
                        if (Abs(b[i]-sum) > T(Tol))
                            throw "No solution.";
                    }
                    else
                        x[i] = (b[i]-sum)/A(i,i);
                }
                return x;
            }
            else if (IsLowerTriangular(A)) { // Lower triangular, forward substition.
                for (size_t i = 0; i < A.NumColumns(); ++i) {
                    // A(i,i)x(i) + A(i,i+1)x(i+1) + ... + A(i,N-1)x(N-1) = b(i)
                    Complex<T> sum;
                    for (size_t j = i+1; j < A.NumColumns(); ++j)
                        sum += A(i,j)*x[j];
                    if (Abs(A(i,i)) < T(Tol)) {
                        x[i] = T(0);
                        if (Abs(b[i]-sum) > T(Tol))
                            throw "No solution.";
                    }
                    else
                        x[i] = (b[i]-sum)/A(i,i);
                }
                return x;
            }
        }

        Matrix<T,M,(N==Dynamic?Dynamic:N+1),Flags> Ab = Augmented(A, b);
        Ab = RREF(Ab);
        size_t next_pivot = A.NumColumns();
        for (int i = Ab.NumRows()-1; i >= 0; i--) {
            size_t pivot = 0;
            for (; pivot < A.NumColumns(); ++pivot) {
                if (Abs(Ab(i, pivot)) > T(Tol))
                    break;
            }
            if (pivot == A.NumColumns()) {
                if (Abs(Ab(i,pivot)) > T(Tol)) // 0x+0y+0z = 1
                    throw "No solution.";
            }
            else if (pivot >= next_pivot) { // Not row reduced.
                throw "No solution.";
            }
            else { // Ab(i,pivot)x(pivot) + Ab(i,pivot+1)x(pivot+1)+ ... + Ab(i,N-1)x(N-1) = Ab(i,N)
                Complex<T> sum;
                for (size_t j = pivot+1; j < A.NumColumns(); ++j)
                    sum += Ab(i,j)*x[j]; // Free variables are filled in with 0.
                x[pivot] = (Ab(i,A.NumColumns())-sum)/Ab(i,pivot);
            }
            next_pivot = pivot;
        }

        return x;
    }
}
