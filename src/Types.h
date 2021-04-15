#pragma once
#include <vector>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"
#include "Construction.h"
#include "Global.h"

namespace Linear {
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsUpperTriangular(const Matrix<T,M,N,Flags>& a) {
        if (!IsSquare(a))
            return false;
        for (size_t c = 0; c < a.NumRows()-1; ++c) {
            for (size_t r = c+1; r < a.NumRows(); ++r) {
                if (a(r,c) != T(0))
                    return false;
            }
        }
        return true;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsLowerTriangular(const Matrix<T,M,N,Flags>& a) {
        if (!IsSquare(a))
            return false;
        for (size_t r = 0; r < a.NumRows()-1; ++r) {
            for (size_t c = r+1; c < a.NumRows(); ++c) {
                if (a(r,c) != T(0))
                    return false;
            }
        }
        return true;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsTriangular(const Matrix<T,M,N,Flags>& a) {
        return (IsUpperTriangular(a) || IsLowerTriangular(a));
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsDiagonal(const Matrix<T,M,N,Flags>& a) {
        return (IsUpperTriangular(a) && IsLowerTriangular(a));
    }
}
