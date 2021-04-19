#pragma once
#include "Matrix.h"
#include "Basics.h"
#include "Vector.h"
#include "Types.h"
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
        std::cout << "B = " << B << std::endl;
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
}
