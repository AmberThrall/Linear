#pragma once
#include <vector>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"
#include "Construction.h"
#include "Global.h"

namespace Linear {
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
}
