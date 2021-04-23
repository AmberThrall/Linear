#pragma once
#include <vector>
#include <tuple>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"
#include "Construction.h"
#include "Eigen.h"
#include "Types.h"
#include "Global.h"

namespace Linear {
    /**
     * Computes the LUP decomposition of A.
     * This finds NxN matrix L, NxN matrix U and NxN matrix P such that \f$PA=LU\f$, L is lower triangular, U is upper triangular
     * and P is a permutation matrix. If A is not square an exception is thrown.
     * @param A MxN matrix
     * @return Tuple (L,U,P)
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    std::tuple<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>> LUP(Matrix<T,M,N,Flags> A) {
        if (!IsSquare(A))
            throw "LUP Decomposition is only defined for square matrices.";
        std::cout << "LUP(A). A = " << A <<std::endl;

        SquareMatrix<T,N,Flags> P(A.NumRows(), T(0));
        for (size_t i = 0; i < A.NumRows(); ++i)
            P(i,i) = T(1);

        for (size_t i = 0; i < A.NumRows(); ++i) {
            // Step 1. Select the pivot
            T Amax = T(0);
            size_t imax = i;
            for (size_t j = i; j < A.NumRows(); ++j) {
                Complex<T> Aii=A(j,i);
                for (size_t q = 0; q < i; ++q) {
                    Aii -= A(j, q)*A(q, j);
                }
                if (Abs(Aii) > Amax) {
                    Amax = Abs(Aii);
                    imax = j;
                }

            }

            // Step 2. Swap rows.
            if (imax != i) {
                P.SwapRows(i, imax);
                A.SwapRows(i, imax);
            }

            for (size_t j = i; j < A.NumRows(); ++j) {
                for (size_t q = 0; q < i; ++q)
                    A(i,j) -= A(i,q)*A(q,j);
            }
            for (size_t j = i+1; j < A.NumRows(); ++j) {
                for (size_t q = 0; q < i; ++q)
                    A(j,i) -= A(j,q)*A(q,i);
                if (Abs(A(i,i)) > T(Tol))
                    A(j,i) = A(j,i)/A(i,i);
            }
        }

        // Get the U matrix.
        SquareMatrix<T,N,Flags> U = A;
        for (size_t c = 0; c < A.NumColumns(); ++c) {
            for (size_t r=c+1; r < A.NumRows(); ++r)
                U(r,c) = T(0);
        }
        // Get the L matrix.
        SquareMatrix<T,N,Flags> L = A - U + Identity<T,Flags>(A.NumRows());
        return std::make_tuple(L, U, P);
    }

    /**
     * Computes the LUP decomposition of A.
     * This finds MxN matrix Q and NxN matrix R such that \f$A=QR\f$ and R is upper triangular. Moreover, if A is square then Q is unitary.
     * If M < N, then an exception is thrown.
     * @param A MxN matrix
     * @return Pair (Q,R)
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M>=N||M==Dynamic||N==Dynamic), std::pair<SquareMatrix<T,N,Flags>,Matrix<T,M,N,Flags>>>::type QR(const Matrix<T,M,N,Flags>& A) {
        if (A.NumRows() < A.NumColumns())
            throw "QR Decomposition is defined for m-by-n matrices where m>=n.";

        // Perform Householder reflections.
        std::vector<SquareMatrix<T,M,Flags>> qs;
        SquareMatrix<T,Dynamic,Flags> aprime = A;
        size_t niters = std::min(A.NumRows()-1,A.NumColumns());
        for (size_t k = 0; k < niters; ++k) {
            SquareMatrix<T,Dynamic,Flags> eye = Identity<T,Flags>(aprime.NumRows());
            Vector<T,Dynamic> e1 = Basis<T>(aprime.NumRows(), 1, 0, 0);
            Vector<T,Dynamic> x = aprime.GetColumn(0);
            Vector<T,Dynamic> u = x - Norm(x)*e1;
            u = Normalize(u);

            SquareMatrix<T,Dynamic,Flags> qprime = eye - 2*u*ConjugateTranspose(u);
            SquareMatrix<T,M,Flags> qk = (k > 0 ? Diag(Identity<T,Flags>(k), qprime) : qprime);
            aprime = RemoveRowAndColumn(qprime*aprime, 0, 0);
            qs.push_back(qk);
        }

        SquareMatrix<T,M,Flags> Q = qs[0];
        for (size_t k = 1; k < qs.size(); ++k) {
            Q = Q*Transpose(qs[k]);
        }
        Matrix<T,M,N,Flags> R = Transpose(Q)*A;

        return std::make_pair(Q, R);
    }

    /**
     * Computes the LL* decomposition of A.
     * This finds NxN matrix L such that \f$A=LL^*\f$ and L is lower triangular with real and positive diagonal entries.
     * If A is not square then an exception is thrown. If A is not a Hermitian positive-definite matrix, the decomposition may fail.
     * @param A MxN matrix
     * @return L
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic), SquareMatrix<T,N,Flags>>::type LL(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            throw "LL* Decomposition is only defined for square matrices.";

        SquareMatrix<T,N,Flags> ret(A.NumRows(), T(0));
        for (size_t j = 0; j < A.NumRows(); ++j) {
            ret(j,j) = A(j,j);
            for (size_t k = 0; k < j; ++k)
                ret(j,j) -= ret(j,k)*Conjugate(ret(j,k));
            ret(j,j) = Sqrt(ret(j,j));

            for (size_t i = j+1; i < A.NumRows(); ++i) {
                ret(i,j) = A(i,j);
                for (size_t k = 0; k < j; ++k)
                    ret(i,j) -= ret(i,k)*Conjugate(ret(j,k));
                ret(i,j) /= ret(j,j);
            }
        }

        return ret;
    }

    /**
     * Computes the LDL* decomposition of A.
     * This finds NxN matrix L and NxN matrix D such that \f$A=LDL^*\f$, L is lower unit triangular and D is diagonal.
     * If A is not square then an exception is thrown.  If A is not a Hermitian positive-definite matrix, the decomposition may fail.
     * @param A MxN matrix
     * @return Pair (L, D)
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic), std::pair<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>>>::type LDL(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            throw "LDL* Decomposition is only defined for square matrices.";

        SquareMatrix<T,N,Flags> L(A.NumRows(), T(0));
        SquareMatrix<T,N,Flags> D(A.NumRows(), T(0));
        for (size_t j = 0; j < A.NumRows(); ++j) {
            D(j,j) = A(j,j);
            for (size_t k = 0; k < j; ++k)
                D(j,j) -= L(j,k)*Conjugate(L(j,k))*D(k,k);

            L(j,j) = 1;
            for (size_t i = j+1; i < A.NumRows(); ++i) {
                L(i,j) = A(i,j);
                for (size_t k = 0; k < j; ++k)
                    L(i,j) -= L(i,k)*Conjugate(L(j,k))*D(k,k);
                L(i,j) /= D(j,j);
            }
        }

        return std::make_pair(L, D);
    }

    /**
     * Computes the Eigendecomposition of A.
     * This finds NxN matrix V and NxN matrix D such that \f$A=VDV^{-1}\f$, D is diagonal with the eigenvalues of A along the diagonal
     * and V is the eigenvectors of A.
     * If A is not square or V is not invertible, then an exception is thrown.
     * @param A MxN matrix
     * @return Tuple \f$(V, D, V^{-1})\f$
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic), std::tuple<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>>>::type
    Eigendecomposition(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            throw "Eigendecomposition is only defined for square matrices.";

        std::vector<std::pair<Complex<T>,Vector<T,N>>> eigens = Eigen(A);

        SquareMatrix<T,N,Flags> V(A.NumRows(), T(0));
        SquareMatrix<T,N,Flags> D(A.NumRows(), T(0));
        for (size_t i = 0; i < A.NumRows(); ++i) {
            V.SetColumn(i, eigens[i].second);
            D(i,i) = eigens[i].first;
        }

        if (Determinant(V) == T(0))
            throw "The matrix is not diagonalizable.";

        return std::make_tuple(V, D, Inverse(V));
    }

    /**
     * Computes the SVD decomposition of A.
     * This finds MxM matrix U, MxN matrix S and NxN matrix V such that \f$A=USV^*\f$ where U and V are both unitary, and S is a rectangular diagonal
     * matrix containing the real non-negative singular values in decreasing order.
     * @param A MxN matrix
     * @return Tuple (U,S,V)
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    std::tuple<SquareMatrix<T,M,Flags>,Matrix<T,M,N,Flags>,SquareMatrix<T,N,Flags>> SVD(const Matrix<T,M,N,Flags>& A) {
        std::vector<std::pair<Complex<T>,Vector<T,M>>> left = Eigen(A*ConjugateTranspose(A));
        std::vector<std::pair<Complex<T>,Vector<T,N>>> right = Eigen(ConjugateTranspose(A)*A);

        SquareMatrix<T,M,Flags> U(A.NumRows(), T(0));
        Matrix<T,M,N,Flags> S(A.NumRows(), A.NumColumns(), T(0));
        SquareMatrix<T,N,Flags> V(A.NumColumns(), T(0));

        // Sort singular values by largest.
        size_t i = 0;
        while (left.size() > 0 && right.size() > 0) {
            std::cout << std::endl;

            size_t ileft = 0, iright = 0;
            Complex<T> sval;
            T abs = T(0);
            for (size_t j = 0; j < std::min(left.size(), right.size()); ++j) {
                if (left.size() <= right.size() && Abs(left[j].first) > abs) {
                    ileft = j;
                    sval = left[j].first;
                    abs = Abs(sval);
                }
                else if (left.size() > right.size() && Abs(right[j].first) > abs) {
                    iright = j;
                    sval = right[j].first;
                    abs = Abs(sval);
                    for (size_t k = 0; k < left.size(); ++k) {
                        if (Abs(sval-left[k].first) < T(Tol)) {
                            ileft = k;
                            break;
                        }
                    }
                }
            }
            if (left.size() <= right.size()) {
                // Find iright by taking the closest to the value.
                abs = T(-1);
                for (size_t k = 0; k < right.size(); ++k) {
                    T diff = Abs(sval - right[k].first);
                    if (diff < abs || abs < 0) {
                        iright = k;
                        abs = diff;
                    }
                }
            }
            else  {
                // Find ileft by taking the closest to the value.
                abs = T(-1);
                for (size_t k = 0; k < left.size(); ++k) {
                    T diff = Abs(sval - left[k].first);
                    if (diff < abs || abs < 0) {
                        ileft = k;
                        abs = diff;
                    }
                }
            }

            // Add it.
            U.SetColumn(i, left[ileft].second);
            S(i,i) = Sqrt(sval);
            V.SetColumn(i, right[iright].second);
            i += 1;

            // Remove entries
            left.erase(left.begin()+ileft);
            right.erase(right.begin()+iright);
        }

        return std::make_tuple(U,S,V);
    }

    /**
     * Reduces the matrix A to Hessenberg.
     * This finds NxN matrix Q and NxN matrix H such that \f$A=QHQ^*\f$ where Q is unitary and H is upper Hessenberg.
     * If A is not square, an exception is thrown.
     * @param A MxN matrix
     * @return Pair (Q,H)
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic), std::pair<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>>>::type
    Hessenberg(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            throw "Hessenberg decomposition requires a squire matrix.";

        SquareMatrix<T,N,Flags> Q = Identity<T>(A.NumRows());
        SquareMatrix<T,N,Flags> H = A;
        SquareMatrix<T,Dynamic,Flags> Aprime = A;

        for (size_t j = 0; j < A.NumRows()-2; ++j) {
            SquareMatrix<T,Dynamic,Flags> P = Householder(Aprime.GetColumn(0), Aprime.NumRows()-2);
            SquareMatrix<T,N,Flags> P2 = (j > 0 ? Diag(Identity<T,Flags>(j), P) : P);
            Aprime = RemoveRowAndColumn(ConjugateTranspose(P)*Aprime*P, 0, 0);
            H = ConjugateTranspose(P2)*H*P2;
            Q = Q*P2;
        }

        return std::make_pair(Q, H);
    }

    /**
     * Computes the real Schur decomposition of A
     * This finds NxN matrix Q and NxN matrix U such that \f$A=QUQ^*\f$ where Q is unitary and U is block upper triangular with 1x1 and 2x2 blocks.
     * The matrix U and A have the same eigenvalues.
     * If A is not square, an exception is thrown.
     * @param A MxN matrix
     * @return Pair (Q,U)
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic), std::pair<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>>>::type
    Schur(const Matrix<T,M,N,Flags>& A, size_t max_iterations = 100) {
        if (!IsSquare(A))
            throw "Schur decomposition requires a squire matrix.";

        SquareMatrix<T,N,Flags> eye = Identity<T>(A.NumRows());
        SquareMatrix<T,N,Flags> Q = Identity<T>(A.NumRows());
        SquareMatrix<T,N,Flags> U = A;
        if (!IsHessenberg(A)) {
            std::pair<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>> pair = Hessenberg(A);
            Q = pair.first;
            U = pair.second;
        }

        for (size_t i = A.NumRows()-1; i >= 1; i--) {
            size_t k = 0;
            while (k < max_iterations && Abs(U(i,i-1)) > T(Tol)) {
                Complex<T> sigma = U(i,i);
                std::pair<SquareMatrix<T,N,Flags>,Matrix<T,M,N,Flags>> qr = QR(U-sigma*eye);
                U = qr.second*qr.first + sigma*eye;
                Q = Q*qr.first;
                k += 1;
            }
        }

        return std::make_pair(Q, U);
    }
}
