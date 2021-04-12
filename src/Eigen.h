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

    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)&&(P==N||P==Dynamic||N==Dynamic)), std::pair<Complex<T>,Vector<T,P>>>::type
    PowerIteration(const Matrix<T,M,N,Flags>& a, Vector<T,P> b, unsigned int num_simulations) {
        if (!IsSquare(a))
            throw "Eigenvalues are only defined for square matrices.";
        if (a.NumRows() != b.NumRows())
            throw "Size mismatch in PowerIteration().";

        Complex<T> lambda;
        for (unsigned int i = 0; i < num_simulations; ++i) {
            Vector<T,P> b_next = a*b;
            lambda = Dot(b,b_next)/Dot(b,b);
            b = b_next/Length(b_next);
        }
        return std::make_pair(lambda, b);
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    std::vector<Vector<T,N>> Nullspace(const Matrix<T,M,N,Flags>& m) {
        std::vector<Vector<T,N>> basis;

        SquareMatrix<T,N,Flags> eye = Identity<T>(m.NumColumns());
        Matrix<T,(M==Dynamic||N==Dynamic?Dynamic:M+N),N,Flags> aug = RowAugmented(m, eye);
        aug = Transpose(RREF(Transpose(aug)));
        SquareMatrix<T,N,Flags> c = SubMatrix(aug, m.NumColumns(), m.NumColumns(), m.NumRows(), 0);
        for (size_t col = 0; col < m.NumColumns(); ++col) {
            bool add = true;
            for (size_t row = 0; row < m.NumRows(); ++row) {
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

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)), std::vector<std::pair<Complex<T>,Vector<T,N>>>>::type
    WielandtDeflationAlgorithm(const Matrix<T,M,N,Flags> a, unsigned int num_simulations = 15) {
        if (!IsSquare(a))
            throw "Eigenvalues are only defined for square matrices.";

        // Final step: 1x1 matrix.
        if (a.NumRows() == 1) {
            std::vector<std::pair<Complex<T>,Vector<T,N>>> eigenpairs;
            eigenpairs.push_back(std::make_pair(a(0,0), Vector<T,1>(T(1))));
            return eigenpairs;
        }

        // Step 1: Computer the dominant eigenvalue and vector.
        std::vector<std::pair<Complex<T>,Vector<T,N>>> eigenpairs;
        Vector<T,N> b = Random<T>(a.NumRows(),1);
        std::pair<Complex<T>,Vector<T,N>> pair = PowerIteration(a, b, num_simulations);
        Complex<T> lambda1 = pair.first;
        Vector<T,N> x1 = pair.second;

        // Step 2: Find the p such x1[p] is maximal.
        size_t p = 0;
        T max = T(0);
        for (size_t i = 0; i < a.NumRows(); ++i) {
            if (max < x1[i].Abs()) {
                p = i;
                max = x1[i].Abs();
            }
        }
        x1 = x1 / x1[p];

        // Step 3/4: Compute Ap and remove row p and column p
        RowVector<T,N> ap = a.GetRow(p);
        SquareMatrix<T,Dynamic,Flags> Ap = RemoveRowAndColumn(a - x1*ap, p, p);

        // Step 5: Repeat.
        eigenpairs.push_back(std::make_pair(lambda1, x1));
        std::vector<std::pair<Complex<T>,Vector<T,Dynamic>>> res = WielandtDeflationAlgorithm(Ap, num_simulations);
        for (size_t i = 0; i < res.size(); ++i) {
            Vector<T,N> yi(a.NumRows(), T(0));
            for (size_t j = 0; j < yi.Size(); ++j) {
                if (j == p)
                    yi[j] = T(0);
                else if (j > p)
                    yi[j] = res[i].second[j-1];
                else
                    yi[j] = res[i].second[j];
            }
            yi = yi + (Dot(ap, yi)/(res[i].first-lambda1))*x1;
            eigenpairs.push_back(std::make_pair(res[i].first, yi));
        }
        return eigenpairs;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)), std::vector<std::pair<Complex<T>,Vector<T,N>>>>::type
    Eigen(Matrix<T,M,N,Flags> a, unsigned int num_simulations = 15) {
        if (!IsSquare(a))
            throw "Eigenvalues are only defined for square matrices.";

        std::vector<std::pair<Complex<T>,Vector<T,N>>> eigenpairs;
        if (IsTriangular(a)) {
            SquareMatrix<T,Dynamic,Flags> eye = Identity<T>(a.NumColumns());
            for (size_t i = 0; i < a.NumRows(); ++i) {
                Complex<T> lambda = a(i,i);
                std::vector<Vector<T,N>> basis = Nullspace(a - lambda*eye);
                eigenpairs.push_back(std::make_pair(lambda, Normalize(basis[0])));
            }
            return eigenpairs;
        }

        // Generic case via Deflation and power iteration.
        return WielandtDeflationAlgorithm(a, num_simulations);
    }
}
