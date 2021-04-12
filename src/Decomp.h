#pragma once
#include <vector>
#include <tuple>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"
#include "Construction.h"
#include "Global.h"

namespace Linear {
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic),
        std::tuple<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>>>::type LUP(Matrix<T,M,N,Flags> a) {
        if (!IsSquare(a))
            throw "LUP Decomposition is only defined for square matrices.";

        std::vector<size_t> perms;
        for (size_t i = 0; i < a.NumRows(); ++i)
            perms.push_back(i);
        for (size_t i = 0; i < a.NumRows(); ++i) {
            T maxA = T(0);
            size_t imax = i;
            for (size_t j = i; j < a.NumRows(); ++j) {
                if (a(j,i).Abs() > maxA) {
                    maxA = a(j,i).Abs();
                    imax = j;
                }
            }

            if (imax != i) {
                size_t k = perms[i];
                perms[i] = perms[imax];
                perms[imax] = k;
                a.SwapRows(i, imax);
            }

            for (size_t j = i+1; j < a.NumRows(); ++j) {
                a(j,i) /= a(i,i);
                for (size_t k = i+1; k < a.NumRows(); ++k)
                    a(j,k) -= a(j,i) * a(i,k);
            }
        }

        // Create the permutation matrix.
        SquareMatrix<T,N,Flags> p(a.NumRows(), T(0));
        for (size_t i = 0; i < perms.size(); ++i)
            p.SetColumn(i, Basis<T>(a.NumRows(), perms[i]));
        p = Transpose(p);
        // Get the U matrix.
        SquareMatrix<T,N,Flags> u = a;
        for (size_t c = 0; c < a.NumColumns(); ++c) {
            for (size_t r=c+1; r < a.NumRows(); ++r)
                u(r,c) = T(0);
        }
        // Get the L matrix.
        SquareMatrix<T,N,Flags> l = a - u + Identity<T>(a.NumRows());
        return std::make_tuple(l, u, p);
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M>=N||M==Dynamic||N==Dynamic), std::pair<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>>>::type QR(const Matrix<T,M,N,Flags>& a) {
        if (a.NumRows() < a.NumColumns())
            throw "QR Decomposition is only defined for m-by-n matrices where m>=n.";

        // Perform Householder reflections.
        std::vector<SquareMatrix<T,M,Flags>> qs;
        SquareMatrix<T,Dynamic,Flags> aprime = a;
        size_t niters = std::min(a.NumRows()-1,a.NumColumns());
        for (size_t k = 0; k < niters; ++k) {
            SquareMatrix<T,Dynamic,Flags> eye = Identity<T,Flags>(aprime.NumRows());
            Vector<T,Dynamic> e1 = Basis<T>(aprime.NumRows(), 0);
            Vector<T,Dynamic> x = aprime.GetColumn(0);
            Vector<T,Dynamic> u = x - Length(x)*e1;
            T len = Length(u);
            if (len >= T(Tol))
                u = u/len;

            SquareMatrix<T,Dynamic,Flags> qprime = eye - 2*u*ConjugateTranspose(u);
            SquareMatrix<T,M,Flags> qk = (k > 0 ? Diag(Identity<T,Flags>(k), qprime) : qprime);
            aprime = RemoveRowAndColumn(qprime*aprime, 0, 0);
            qs.push_back(qk);
        }

        SquareMatrix<T,M,Flags> q = qs[0];
        for (size_t k = 1; k < qs.size(); ++k) {
            q = q*Transpose(qs[k]);
        }
        Matrix<T,M,N,Flags> r = Transpose(q)*a;

        return std::make_pair(q, r);
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic), SquareMatrix<T,N,Flags>>::type LL(const Matrix<T,M,N,Flags>& a) {
        if (!IsSquare(a))
            throw "LL* Decomposition is only defined for square matrices.";

        SquareMatrix<T,N,Flags> ret(a.NumRows(), T(0));
        for (size_t j = 0; j < a.NumRows(); ++j) {
            ret(j,j) = a(j,j);
            for (size_t k = 0; k < j; ++k)
                ret(j,j) -= ret(j,k)*ret(j,k).Conjugate();
            ret(j,j) = Complex<T>::Sqrt(ret(j,j));

            for (size_t i = j+1; i < a.NumRows(); ++i) {
                ret(i,j) = a(i,j);
                for (size_t k = 0; k < j; ++k)
                    ret(i,j) -= ret(i,k)*ret(j,k).Conjugate();
                ret(i,j) /= ret(j,j);
            }
        }

        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic), std::pair<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>>>::type LDL(const Matrix<T,M,N,Flags>& a) {
        if (!IsSquare(a))
            throw "LDL* Decomposition is only defined for square matrices.";

        SquareMatrix<T,N,Flags> l(a.NumRows(), T(0));
        SquareMatrix<T,N,Flags> d(a.NumRows(), T(0));
        for (size_t j = 0; j < a.NumRows(); ++j) {
            d(j,j) = a(j,j);
            for (size_t k = 0; k < j; ++k)
                d(j,j) -= l(j,k)*l(j,k).Conjugate()*d(k,k);

            l(j,j) = 1;
            for (size_t i = j+1; i < a.NumRows(); ++i) {
                l(i,j) = a(i,j);
                for (size_t k = 0; k < j; ++k)
                    l(i,j) -= l(i,k)*l(j,k).Conjugate()*d(k,k);
                l(i,j) /= d(j,j);
            }
        }

        return std::make_pair(l, d);
    }
}
