#pragma once
#include "Matrix.h"
#include "Basics.h"
#include "Vector.h"
#include "Types.h"
#include "Global.h"

namespace Linear {
    /**
     * Computes a basis for the nullspace of A.
     * The nullspace of a matrix is the set of vectors v such that Av=0.
     * @param A MxN matrix
     * @return List of vectors v such that \f$Null(A)=span\{v[0],\dots,v[len(v)-1]\}\f$
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    std::vector<Vector<T,N>> Nullspace(const Matrix<T,M,N,Flags>& A) {
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
     * Computes the dimension of the A's nullspace.
     * @param A MxN matrix
     * @return dim(Null(A))
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    unsigned int Nullity(const Matrix<T,M,N,Flags>& A) {
        return Nullspace(A).size();
    }
}
