#pragma once
#include <vector>
#include "Matrix.h"
#include "Global.h"

namespace Linear {
    /**
     * Returns the PxQ submatrix located at (i,j)
     * If i+P>M or j+Q > N, an exception is raised.
     *
     * Example: B is the matrix {{2,3},{5,6}}
     *
     *     Matrix3d A = { {1,2,3}, {4,5,6}, {7,8,9} };
     *     Matrix2d B = SubMatrix<2,2>(A, 0, 1);
     *
     * @param A MxN matrix
     * @param i Row offset (default = 0)
     * @param j Column offset (default = 0)
     * @return PxQ submatrix of A
     */
    template <size_t P, size_t Q, typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(P>0&&Q>0),Matrix<T,P,Q,Flags>>::type SubMatrix(Matrix<T,M,N,Flags> A, size_t i = 0, size_t j = 0) {
        if (i+P > A.NumRows() || j+Q > A.NumColumns())
            throw "Cannot create submatrix, indices out of bounds.";
        Matrix<T,P,Q,Flags> ret(P,Q,T(0));
        for (size_t r = 0; r < P; ++r) {
            for (size_t c = 0; c < Q; ++c) {
                ret(r,c) = A(i+r,j+c);
            }
        }
        return ret;
    }
    /**
     * Returns the nrowsxncols submatrix located at (i,j)
     * If i+nrows>M or j+ncols > N, an exception is raised.
     *
     * Example: B is the matrix {{2,3},{5,6}}
     *
     *     Matrix3d A = { {1,2,3}, {4,5,6}, {7,8,9} };
     *     MatrixXd B = SubMatrix(A, 2, 2, 0, 1);
     *
     * @param A MxN matrix
     * @param nrows Number of rows in submatrix
     * @param ncols Number of columns in submatrix
     * @param i Row offset (default = 0)
     * @param j Column offset (default = 0)
     * @return nrowsxncols submatrix of A
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,Dynamic,Dynamic,Flags> SubMatrix(Matrix<T,M,N,Flags> A, size_t nrows, size_t ncols, size_t i = 0, size_t j = 0) {
        if (i+nrows > A.NumRows() || j+ncols > A.NumColumns())
            throw "Cannot create submatrix, indices out of bounds.";
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows,ncols,T(0));
        for (size_t r = 0; r < nrows; ++r) {
            for (size_t c = 0; c < ncols; ++c) {
                ret(r,c) = A(i+r,j+c);
            }
        }
        return ret;
    }
    /**
     * Returns the (M-1)xN submatrix formed by removing the ith row.
     * @param A MxN matrix
     * @param i Row to remove
     * @return (M-1)xN submatrix of A
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,(M==Dynamic?Dynamic:M-1),N,Flags> RemoveRow(Matrix<T,M,N,Flags> A, size_t i) {
        if (A.NumRows() == 1)
            throw "Cannot create a 0xN matrix.";
        Matrix<T,(M==Dynamic?Dynamic:M-1),N,Flags> ret(A.NumRows()-1,A.NumColumns(),T(0));
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                if (r < i)
                    ret(r,c) = A(r,c);
                if (r > i)
                    ret(r-1,c) = A(r,c);
            }
        }
        return ret;
    }
    /**
     * Returns the Mx(N-1) submatrix formed by removing the ith column.
     * @param A MxN matrix
     * @param i Column to remove
     * @return Mx(N-1) submatrix of A
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,(N==Dynamic?Dynamic:N-1),Flags> RemoveColumn(Matrix<T,M,N,Flags> A, size_t i) {
        if (A.NumColumns() == 1)
            throw "Cannot create a Mx0 matrix.";
        Matrix<T,M,(N==Dynamic?Dynamic:N-1),Flags> ret(A.NumRows(),A.NumColumns()-1,T(0));
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                if (c < i)
                    ret(r,c) = A(r,c);
                if (c > i)
                    ret(r,c-1) = A(r,c);
            }
        }
        return ret;
    }
    /**
     * Returns the (M-1)x(N-1) submatrix formed by removing the ith row and jth column.
     * @param A MxN matrix
     * @param i Row to remove
     * @param j Column to remove
     * @return (M-1)x(N-1) submatrix of A
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,(M==Dynamic?Dynamic:M-1),(N==Dynamic?Dynamic:N-1),Flags> RemoveRowAndColumn(Matrix<T,M,N,Flags> A, size_t i, size_t j) {
        if (A.NumRows() == 1 || A.NumColumns() == 1)
            throw "Cannot create a 0x0 matrix.";
        Matrix<T,(M==Dynamic?Dynamic:M-1),(N==Dynamic?Dynamic:N-1),Flags> ret(A.NumRows()-1,A.NumColumns()-1,T(0));
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                if (r < i && c < j)
                    ret(r,c) = A(r,c);
                if (r > i && c < j)
                    ret(r-1,c) = A(r,c);
                if (r < i && c > j)
                    ret(r,c-1) = A(r,c);
                if (r > i && c > j)
                    ret(r-1,c-1) = A(r,c);
            }
        }
        return ret;
    }

    /**
     * Returns the NxM matrix \f$B=A^\top\f$ defined by \f$b_{ij}=a_{ji}\f$.
     * @param A MxN matrix
     * @return NxM matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,N,M,Flags> Transpose(const Matrix<T,M,N,Flags>& A) {
        Matrix<T,N,M,Flags> ret(A.NumColumns(), A.NumRows(), T(0));
        for (size_t r = 0; r < A.NumColumns(); ++r) {
            for (size_t c = 0; c < A.NumRows(); ++c) {
                ret(r,c) = A(c,r);
            }
        }
        return ret;
    }

    /**
     * Returns the NxM matrix \f$B=A^*\f$ defined by \f$b_{ij}=\overline{a_{ji}}\f$.
     * @param A MxN matrix
     * @return NxM matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,N,M,Flags> ConjugateTranspose(Matrix<T,M,N,Flags> A) {
        Matrix<T,N,M,Flags> ret(A.NumColumns(), A.NumRows(), T(0));
        for (size_t r = 0; r < A.NumColumns(); ++r) {
            for (size_t c = 0; c < A.NumRows(); ++c) {
                ret(r,c) = Conjugate(A(c,r));
            }
        }
        return ret;
    }

    /**
     * Computes A's reduced row echelon form.
     * @param A MxN matrix
     * @return MxN matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> RREF(Matrix<T,M,N,Flags> A) {
        size_t lead = 0;
        for (size_t r = 0; r < A.NumRows(); ++r) {
            if (lead >= A.NumColumns())
                return A;

            // Find the pivot.
            size_t i = r;
            while (A(i, lead) == 0) {
                i += 1;
                if (i == A.NumRows()) {
                    i = r;
                    lead += 1;
                    if (lead == A.NumColumns())
                        return A;
                }
            }

            // Swap rows i and r
            A.SwapRows(i, r);
            // R_r / pivot -> R_r
            A.ScaleRow(r, 1/A(r,lead));
            // For each row i!=r, R_i - R_r*m(i,lead) -> R_i
            for (size_t i = 0; i < A.NumRows(); ++i) {
                if (i == r) continue;

                A.AddRows(i, r, -A(i,lead));
            }

            lead += 1;
        }
        return A;
    }

    /**
     * Computes the trace of a square matrix \f$\sum_{i=0}^{N-1}a_{ii}\f$. If A is not square, an exception is thrown.
     * @param A MxN matrix
     * @return Complex number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Complex<T>>::type Trace(const Matrix<T,M,N,Flags>& A) {
        if (A.NumRows() != A.NumColumns())
            throw "Cannot take the trace of a non-square matrix.";
        Complex<T> ret;
        for (size_t r = 0; r < A.NumRows(); ++r) {
            ret += A(r,r);
        }
        return ret;
    }

    /**
     * Computes the determinant of a square matrix defined by \f$\sum_{i=0}^{N-1}(-1)^ia_{0i}\det(A_{0i})\f$ where \f$A_{0i}\f$ is the (N-1)x(N-1)
     * matrix obtained by deleting the 0th row and ith column. If A is not square, an exception is thrown.
     * @param A MxN matrix
     * @return Complex number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Complex<T>>::type Determinant(const Matrix<T,M,N,Flags>& A) {
        if (A.NumRows() != A.NumColumns())
            throw "Cannot take the determinant of a non-square matrix.";

        if (A.NumRows() == 1)
            return A(0,0);
        if (A.NumRows() == 2)
            return (A(0,0)*A(1,1) - A(0,1)*A(1,0));

        Complex<T> ret;
        for (size_t i = 0; i < A.NumRows(); ++i) {
            if (i % 2 == 0)
                ret += A(0,i)*Determinant(RemoveRowAndColumn(A, 0, i));
            else
                ret += A(0,i)*Determinant(RemoveRowAndColumn(A, 0, i));
        }
        return ret;
    }

    /**
     * Computes the (i,j)-minor of a square matrix defined by \f$\det(A_{ij})\f$ where \f$A_{ij}\f$ is the (N-1)x(N-1)
     * matrix obtained by deleting the ith row and jth column. If A is not square, an exception is thrown. Likewise if i or j is out-of-bounds
     * an exception is thrown.
     * @param A MxN matrix
     * @param i Row index
     * @param j Column index
     * @return Complex number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Complex<T>>::type Minor(const Matrix<T,M,N,Flags>& A, size_t i, size_t j) {
        if (A.NumRows() != A.NumColumns())
            throw "Cannot take the minor of a non-square matrix.";
        if (i >= A.NumRows() || j >= A.NumRows())
            throw "Minor coordinates out of bounds.";
        if (A.NumRows() == 1)
            return A(0,0);

        return Determinant(RemoveRowAndColumn(A, i, j));
    }

    /**
     * Computes the (i,j)-cofactor of a square matrix defined by \f$(-1)^{i+j}\det(A_{ij})\f$ where \f$A_{ij}\f$ is the (N-1)x(N-1)
     * matrix obtained by deleting the ith row and jth column. If A is not square, an exception is thrown. Likewise if i or j is out-of-bounds
     * an exception is thrown.
     * @param A MxN matrix
     * @param i Row index
     * @param j Column index
     * @return Complex number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Complex<T>>::type Cofactor(const Matrix<T,M,N,Flags>& A, size_t i, size_t j) {
        if (A.NumRows() != A.NumColumns())
            throw "Cannot take the cofactor of a non-square matrix.";
        if (i >= A.NumRows() || j >= A.NumRows())
            throw "Cofactor coordinates out of bounds.";
        if (i+j % 2 == 0)
            return Determinant(RemoveRowAndColumn(A, i, j));
        return -Determinant(RemoveRowAndColumn(A, i, j));
    }

    /**
     * Computes the NxN adjugate B of a square matrix A defined by \f$b_{ij}=c_{ji}\f$ where \f$c_{ij}\f$ represents the (i,j)-cofactor of A.
     * If A is not square, an exception is thrown.
     * @param A MxN matrix
     * @return NxN matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Matrix<T,M,N,Flags>>::type Adjugate(const Matrix<T,M,N,Flags>& A) {
        if (A.NumRows() != A.NumColumns())
            throw "Cannot take the adjugate of a non-square matrix.";
        Matrix<T,M,N,Flags> ret(A.NumRows(), A.NumColumns(), T(0));
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                ret(c,r) = Cofactor(A, r, c);
            }
        }
        return ret;
    }

    /**
     * Computes the NxN inverse \f$A^{-1}=\frac{1}{\det A}Adj\f$ of a square matrix A where Adj is the adjugate of A.
     * If A is not square or if A is singular, an exception is thrown.
     * @param A MxN matrix
     * @return NxN matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Matrix<T,M,N,Flags>>::type Inverse(const Matrix<T,M,N,Flags>& A) {
        if (A.NumRows() != A.NumColumns())
            throw "Cannot take the inverse of a non-square matrix.";
        Complex<T> det = Determinant(A);
        if (det == 0)
            throw "Cannot take the inverse of a singular matrix.";
        return Adjugate(A)/det;
    }
}
