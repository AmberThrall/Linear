#pragma once
#include "Matrix.h"
#include "Global.h"

namespace Linear {
    /**
     * Checks if a matrix is a row or column vector.
     * @param A MxN matrix
     * @return True if M=1 or N=1
     */
    template<typename T, size_t M, size_t N, unsigned int Flags>
    bool IsVector(const Matrix<T,M,N,Flags>& A) {
        return (A.NumColumns() == 1 || A.NumRows() == 1);
    }

    /**
     * Checks if the matrix is squre.
     * @param A MxN matrix
     * @return True if M=N.
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsSquare(Matrix<T,M,N,Flags> A) {
        return (A.NumRows() == A.NumColumns());
    }

    /**
     * Checks if A is upper triangular form.
     * @param A MxN matrix
     * @return True if A is upper triangular form
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsUpperTriangular(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            return false;
        for (size_t c = 0; c < A.NumRows()-1; ++c) {
            for (size_t r = c+1; r < A.NumRows(); ++r) {
                if (A(r,c) != T(0))
                    return false;
            }
        }
        return true;
    }
    /**
     * Checks if A is lower triangular form.
     * @param A MxN matrix
     * @return True if A is lower triangular form
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsLowerTriangular(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            return false;
        for (size_t r = 0; r < A.NumRows()-1; ++r) {
            for (size_t c = r+1; c < A.NumRows(); ++c) {
                if (A(r,c) != T(0))
                    return false;
            }
        }
        return true;
    }
    /**
     * Checks if A is lower triangular form or upper triangular form.
     * @param A MxN matrix
     * @return True if A is upper triangular form or lower triangular form
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsTriangular(const Matrix<T,M,N,Flags>& A) {
        return (IsUpperTriangular(A) || IsLowerTriangular(A));
    }
    /**
     * Checks if A is a diagonal matrix.
     * @param A MxN matrix
     * @return True if A is diagonal
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsDiagonal(const Matrix<T,M,N,Flags>& A) {
        return (IsUpperTriangular(A) && IsLowerTriangular(A));
    }
    /**
     * Checks if A is the identity matrix.
     * @param A MxN matrix
     * @return True if A is the identity
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsIdentity(const Matrix<T,M,N,Flags>& A) {
        if (!IsDiagonal(A))
            return false;
        for (size_t i = 0; i < A.NumRows(); ++i) {
            if (Abs(A(i,i)-T(1)) > T(Tol))
                return false;
        }
        return true;
    }
    /**
     * Checks if A is a companion matrix
     * @param A MxN matrix
     * @return True if A is a companion matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsCompanion(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            return false;

        for (size_t i = 0; i < A.NumColumns()-1; ++i) {
            if (Abs(A(0,i)) > T(Tol))
                return false;
        }
        for (size_t r = 1; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns()-1; ++c) {
                if (r-1 == c && Abs(A(r,c)-T(1)) > T(Tol))
                    return false;
                if (r-1 != c && Abs(A(r,c)) > T(Tol))
                    return false;
            }
        }
        return true;
    }

    /**
     * Checks if a matrix is symmetric.
     * @param A MxN matrix
     * @return True if A is square and \f$a_{ij}=a_{ji}\f$ for all i and j.
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsSymmetric(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            return false;
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < r+1; ++c) {
                if (Abs(A(r,c)-A(c,r)) > T(Tol))
                    return false;
            }
        }
        return true;
    }
    /**
     * Checks if a matrix is Hermitian.
     * @param A MxN matrix
     * @return True if A is square and \f$a_{ij}=\overline{a_{ji}}\f$ for all i and j.
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsHermitian(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            return false;
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < r+1; ++c) {
                if (Abs(A(r,c)-Conjugate(A(c,r))) > T(Tol))
                    return false;
            }
        }
        return true;
    }

    /**
     * Checks if a matrix is upper Hessenberg.
     * @param A MxN matrix
     * @return True if A is square and \f$a_{ij}=0\f$ for all i,j such that \f$i>j+1\f$
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsUpperHessenberg(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            return false;
        if (A.NumColumns() == 1 || A.NumColumns() == 2)
            return true;
        for (size_t c = 0; c < A.NumColumns(); ++c) {
            for (size_t r = c+2; r < A.NumRows(); ++r) {
                if (Abs(A(r,c)) > T(Tol))
                    return false;
            }
        }
        return true;
    }
    /**
     * Checks if a matrix is lower Hessenberg.
     * @param A MxN matrix
     * @return True if A is square and \f$a_{ij}=0\f$ for all i,j such that \f$j>i+1\f$
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsLowerHessenberg(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            return false;
        if (A.NumColumns() == 1 || A.NumColumns() == 2)
            return true;
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = r+2; c < A.NumColumns(); ++c) {
                if (Abs(A(r,c)) > T(Tol))
                    return false;
            }
        }
        return true;
    }
    /**
     * Checks if a matrix is Hessenberg.
     * @param A MxN matrix
     * @return True if A is either lower or upper Hessenberg
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsHessenberg(const Matrix<T,M,N,Flags>& a) {
        return (IsLowerHessenberg(a) || IsUpperHessenberg(a));
    }
    /**
     * Checks if a matrix is Tridiagonal.
     * @param A MxN matrix
     * @return True if A is both lower and upper Hessenberg
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsTridiagonal(const Matrix<T,M,N,Flags>& a) {
        return (IsLowerHessenberg(a) && IsUpperHessenberg(a));
    }
}
