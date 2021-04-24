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
    /// \cond DO_NOT_DOCUMENT
    template <typename T, size_t N>
    struct Eigenpair;
    /// \endcond

    /**
     * Struct for LUP decomposition.
     * This struct finds NxN matrix L, NxN matrix U and NxN matrix P such that \f$PA=LU\f$, L is lower triangular, U is upper triangular
     * and P is a permutation matrix.
     * @param T Type to store matrix entries as.
     * @param N Size of the L,U,P matrices (all three are square). Dynamic is allowed for N.
     * @param Flags Flags to pass to the matrices (default = row major).
     */
    template <typename T, size_t N, unsigned int Flags = 0>
    struct LUP {
        SquareMatrix<T,N,Flags> L; /*!< Lower triangular matrix with ones along the diagonal */
        SquareMatrix<T,N,Flags> U; /*!< Upper triangular matrix */
        SquareMatrix<T,N,Flags> P; /*!< Permutation matrix */

        /**
         * Constructor. Just calls Compute.
         * @param A PxQ Matrix
         */
        template <size_t P, size_t Q, unsigned int Flags2>
        LUP(const Matrix<T,P,Q,Flags2>& A) {
            Compute(A);
        }
        /**
         * Computes the LUP decomposition. If A is not square or Q != N, then an exception is thrown.
         * @param A PxQ Matrix
         */
        template <size_t P, size_t Q, unsigned int Flags2>
        void Compute(Matrix<T,P,Q,Flags2> A) {
            if (!IsSquare(A))
                throw "LUP Decomposition is only defined for square matrices.";
            if (A.NumRows() != N && N != Dynamic)
                throw "Cannot perform LUP Decomposition; size mismatch.";

            this->P = Identity<T,Flags>(A.NumRows());
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
                    this->P.SwapRows(i, imax);
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
            this->U = A + Zero<T,Flags>(A.NumRows(), A.NumRows()); // to avoid move operator.
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                for (size_t r=c+1; r < A.NumRows(); ++r)
                    this->U(r,c) = T(0);
            }
            // Get the L matrix.
            this->L = A - this->U + Identity<T,Flags>(A.NumRows());
        }
    };

    /**
     * Struct for QR decomposition.
     * This struct finds MxN matrix Q, NxN matrix R such that \f$A=QR\f$ and R is upper triangular. Moreover, if Q is square then Q is unitary.
     * @param T Type to store matrix entries as.
     * @param M Number of rows for Q. Dynamic is allowed for M.
     * @param N Number of columns for Q and R (R is square). Dynamic is allowed for N.
     * @param Flags Flags to pass to the matrices (default = row major).
    */
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    struct QR {
        Matrix<T,M,N,Flags> Q; /*!< Unitary when square */
        SquareMatrix<T,N,Flags> R; /*!< Upper triangular matrix */

        /**
        * Constructor. Just calls Compute.
        * @param A PxQ Matrix
        */
        template <size_t P, size_t Q, unsigned int Flags2>
        QR(const Matrix<T,P,Q,Flags2>& A) {
            Compute(A);
        }
        /**
        * Computes the QR decomposition. If P != M or Q != N, then an exception is thrown.
        * @param A PxQ Matrix
        */
        template <size_t P, size_t Q, unsigned int Flags2>
        void Compute(Matrix<T,P,Q,Flags2> A) {
            if (A.NumRows() < A.NumColumns())
                throw "QR Decomposition is defined for m-by-n matrices where m>=n.";
            if ((A.NumRows() != M && M != Dynamic) || (A.NumColumns() != N && N != Dynamic))
                throw "Cannot perform QR decomposition; size mismatch.";

            // Perform Householder reflections.
            std::vector<SquareMatrix<T,M,Flags>> qs;
            SquareMatrix<T,Dynamic,Flags> aprime = A;
            size_t niters = std::min(A.NumRows()-1,A.NumColumns());
            for (size_t k = 0; k < niters; ++k) {
                SquareMatrix<T,Dynamic,Flags> qprime = Householder(aprime.GetColumn(0), aprime.NumRows()-1);
                SquareMatrix<T,M,Flags> qk = (k > 0 ? Diag(Identity<T,Flags>(k), qprime) : qprime);
                aprime = RemoveRowAndColumn(qprime*aprime, 0, 0);
                qs.push_back(qk);
            }

            this->Q = qs[0];
            for (size_t k = 1; k < qs.size(); ++k) {
                this->Q = this->Q*Transpose(qs[k]);
            }
            this->R = Transpose(this->Q)*A;
        }
    };

    /**
     * Struct for Cholesky decomposition.
     * This struct finds NxN matrix L and NxN matrix D such that \f$A=LDL^*\f$, L is lower unit triangular and D is diagonal.
     * @param T Type to store matrix entries as.
     * @param N Number of rows/columns for L and D (both square). Dynamic is allowed for N.
     * @param Flags Flags to pass to the matrices (default = row major).
    */
    template <typename T, size_t N, unsigned int Flags = 0>
    struct Cholesky {
        SquareMatrix<T,N,Flags> L; /*!< Lower unit triangular */
        SquareMatrix<T,N,Flags> D; /*!< Diagonal */
        SquareMatrix<T,N,Flags> Lh; /*!< Lh=ConjugateTranspose(L) */

        /**
        * Constructor. Just calls Compute.
        * @param A PxQ Matrix
        */
        template <size_t P, size_t Q, unsigned int Flags2>
        Cholesky(const Matrix<T,P,Q,Flags2>& A) {
            Compute(A);
        }
        /**
        * Computes the Cholesky decomposition. If A is not square or Q != N, then an exception is thrown. This process
        * may fail if A is not a Hermitian positive-definite matrix.
        * @param A PxQ Matrix
        */
        template <size_t P, size_t Q, unsigned int Flags2>
        void Compute(const Matrix<T,P,Q,Flags2>& A) {
            if (!IsSquare(A))
                throw "Cholesky decomposition is defined for Hermitian positive-definite matrices.";
            if (A.NumColumns() != N && N != Dynamic)
                throw "Cannot perform Cholesky decomposition; size mismatch.";

            this->L = Zero<T,Flags>(A.NumRows(), A.NumRows());
            this->D = Zero<T,Flags>(A.NumRows(), A.NumRows());
            for (size_t j = 0; j < A.NumRows(); ++j) {
                this->D(j,j) = A(j,j);
                for (size_t k = 0; k < j; ++k)
                    this->D(j,j) -= this->L(j,k)*Conjugate(this->L(j,k))*this->D(k,k);

                this->L(j,j) = 1;
                for (size_t i = j+1; i < A.NumRows(); ++i) {
                    this->L(i,j) = A(i,j);
                    for (size_t k = 0; k < j; ++k)
                        this->L(i,j) -= this->L(i,k)*Conjugate(this->L(j,k))*this->D(k,k);
                    this->L(i,j) /= this->D(j,j);
                }
            }
            this->Lh = ConjugateTranspose(this->L);
        }
    };

    /**
     * Struct for Eigendecomposition.
     * This struct finds NxN matrix Q and NxN matrix D such that \f$A=QDQ^{-1}\f$, D is diagonal with entries equaling the eigenvalues of A
     * and V's columns are the corresponding eigenvectors of A. Only diagonalizable matrices can be decomposed in this fashion.
     * @param T Type to store matrix entries as.
     * @param N Number of rows/columns for L and D (both square). Dynamic is allowed for N.
     * @param Flags Flags to pass to the matrices (default = row major).
    */
    template <typename T, size_t N, unsigned int Flags = 0>
    struct Eigendecomposition {
        SquareMatrix<T,N,Flags> Q; /*!< Comprised of the eigenvectors of A */
        SquareMatrix<T,N,Flags> D; /*!< Diagonal matrix containing the eigenvalues of A */
        SquareMatrix<T,N,Flags> Qinv; /*!< Qinv=Inverse(Q) */

        /**
        * Constructor. Just calls Compute.
        * @param A M2xN2 Matrix
        */
        template <size_t M2, size_t N2, unsigned int Flags2>
        Eigendecomposition(const Matrix<T,M2,N2,Flags2>& A) {
            Compute(A);
        }
        /**
        * Computes the Eigendecomposition. If A is not square or Q != N, then an exception is thrown. If the matrix is not
        * diagonalizable, an exception will be thrown.
        * @param A M2xN2 Matrix
        */
        template <size_t M2, size_t N2, unsigned int Flags2>
        void Compute(const Matrix<T,M2,N2,Flags2>& A) {
            if (!IsSquare(A))
                throw "Eigendecomposition is defined for square matrices.";
            if (A.NumColumns() != N && N != Dynamic)
                throw "Cannot perform Cholesky decomposition; size mismatch.";

            std::vector<Eigenpair<T,N>> eigens = Eigen(A);

            this->Q = Zero<T,Flags>(A.NumRows(), A.NumRows());
            this->D = Zero<T,Flags>(A.NumRows(), A.NumRows());
            for (size_t i = 0; i < A.NumRows(); ++i) {
                this->Q.SetColumn(i, eigens[i].vector);
                this->D(i,i) = eigens[i].value;
            }

            if (Determinant(this->Q) == T(0))
                throw "The matrix is not diagonalizable.";
            this->Qinv = Inverse(this->Q);
        }
    };

    enum SVDType {
        FULL_SVD,
        THIN_SVD
    };
    /**
     * Struct for SVD decomposition.
     * This struct finds MxM matrix U, MxN matrix S and NxN matrix V such that \f$A=USV^*\f$ where U and V are both unitary, and S is a rectangular
     * diagonal matrix containing the real non-negative singular values of A in decreasing order.
     * If THIN_SVD is passed, then U is MxK, S is KxK, and Vh is KxN where K=min(M,N).
     * @param T Type to store matrix entries as.
     * @param M Number of rows for U and S. Dynamic is allowed for M.
     * @param N Number of columns for S and Vh. Dynamic is allowed for N.
     * @param Flags Flags to pass to the matrices (default = row major).
     * @param Type Either FULL_SVD or THIN_SVD (default = FULL_SVD).
    */
    template <typename T, size_t M, size_t N, unsigned int Flags = 0, SVDType Type = FULL_SVD>
    struct SVD {
        std::conditional_t<Type==FULL_SVD, SquareMatrix<T,M,Flags>, Matrix<T,M,(M>N?N:M),Flags>> U; /*!< Unitary matrix containing the left singular vectors */
        std::conditional_t<Type==FULL_SVD, Matrix<T,M,N,Flags>, SquareMatrix<T,(M>N?N:M),Flags>> S; /*!< Rectangular diagonal matrix containing the singular values of A in decreasing order */
        std::conditional_t<Type==FULL_SVD, SquareMatrix<T,N,Flags>, Matrix<T,(M>N?N:M),N,Flags>> Vh; /*!< Vh=V* */

        /**
        * Constructor. Just calls Compute.
        * @param A PxQ Matrix
        */
        template <size_t P, size_t Q, unsigned int Flags2>
        SVD(const Matrix<T,P,Q,Flags2>& A) {
            Compute(A);
        }
        /**
        * Computes the SVD decomposition. If P != M or Q != N, then an exception is thrown. If the matrix is not
        * diagonalizable, an exception will be thrown.
        * @param A PxQ Matrix
        */
        template <size_t P, size_t Q, unsigned int Flags2>
        void Compute(const Matrix<T,P,Q,Flags2>& A) {
            if ((A.NumRows() != M && M != Dynamic) || (A.NumColumns() != N && N != Dynamic))
                throw "Cannot perform SVD decomposition; size mismatch.";

            std::vector<Eigenpair<T,M>> left = Eigen(A*ConjugateTranspose(A));
            std::vector<Eigenpair<T,N>> right = Eigen(ConjugateTranspose(A)*A);

            size_t k = (A.NumRows() > A.NumColumns() ? A.NumColumns() : A.NumRows());
            std::conditional_t<Type==FULL_SVD, SquareMatrix<T,N,Flags>, Matrix<T,N,(M>N?N:M),Flags>> V = Zero<T,Flags>(A.NumColumns(), (Type==FULL_SVD?A.NumColumns():k));
            if (Type == FULL_SVD) {
                this->U = Zero<T,Flags>(A.NumRows(), A.NumRows());
                this->S = Zero<T,Flags>(A.NumRows(), A.NumColumns());
            }
            else {
                this->U = Zero<T,Flags>(A.NumRows(), k);
                this->S = Zero<T,Flags>(k,k);
            }

            // Sort singular values by largest.
            size_t i = 0;
            while (i < k && left.size() > 0 && right.size() > 0) {
                size_t ileft = 0, iright = 0;
                Complex<T> sval;
                T abs = T(0);
                for (size_t j = 0; j < std::min(left.size(), right.size()); ++j) {
                    if (left.size() <= right.size() && Abs(left[j].value) > abs) {
                        ileft = j;
                        sval = left[j].value;
                        abs = Abs(sval);
                    }
                    else if (left.size() > right.size() && Abs(right[j].value) > abs) {
                        iright = j;
                        sval = right[j].value;
                        abs = Abs(sval);
                        for (size_t k = 0; k < left.size(); ++k) {
                            if (Abs(sval-left[k].value) < T(Tol)) {
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
                        T diff = Abs(sval - right[k].value);
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
                        T diff = Abs(sval - left[k].value);
                        if (diff < abs || abs < 0) {
                            ileft = k;
                            abs = diff;
                        }
                    }
                }

                // Add it.
                U.SetColumn(i, left[ileft].vector);
                S(i,i) = Sqrt(sval);
                V.SetColumn(i, right[iright].vector);
                i += 1;

                // Remove entries
                left.erase(left.begin()+ileft);
                right.erase(right.begin()+iright);
            }

            this->Vh = ConjugateTranspose(V);
        }
    };

    /**
     * Struct for Hessenberg decomposition.
     * This struct NxN matrix Q and NxN matrix H such that \f$A=QHQ^*\f$ where Q is unitary and H is upper Hessenberg.
     * @param T Type to store matrix entries as.
     * @param N Number of rows/columns for Q and H (both square). Dynamic is allowed for N.
     * @param Flags Flags to pass to the matrices (default = row major).
    */
    template <typename T, size_t N, unsigned int Flags = 0>
    struct Hessenberg {
        SquareMatrix<T,N,Flags> Q; /*!< Unitary matrix */
        SquareMatrix<T,N,Flags> H; /*!< Upper Hessenberg matrix */
        SquareMatrix<T,N,Flags> Qh; /*!< Qh=ConjugateTranspose(Q) */

        /**
        * Constructor. Just calls Compute.
        * @param A M2xN2 Matrix
        */
        template <size_t M2, size_t N2, unsigned int Flags2>
        Hessenberg(const Matrix<T,M2,N2,Flags2>& A) {
            Compute(A);
        }
        /**
        * Computes the Hessenberg decomposition. If A is not square or Q != N, then an exception is thrown.
        * @param A M2xN2 Matrix
        */
        template <size_t M2, size_t N2, unsigned int Flags2>
        void Compute(const Matrix<T,M2,N2,Flags2>& A) {
            if (!IsSquare(A))
                throw "Hessenberg decomposition requires a square matrix.";
            if (A.NumColumns() != N && N != Dynamic)
                throw "Cannot perform Hessenberg decomposition; size mismatch.";

            this->Q = Identity<T>(A.NumRows());
            this->H = A;
            SquareMatrix<T,Dynamic,Flags> Aprime = A;

            for (size_t j = 0; j < A.NumRows()-2; ++j) {
                SquareMatrix<T,Dynamic,Flags> P = Householder(Aprime.GetColumn(0), Aprime.NumRows()-2);
                SquareMatrix<T,N,Flags> P2 = (j > 0 ? Diag(Identity<T,Flags>(j), P) : P);
                Aprime = RemoveRowAndColumn(ConjugateTranspose(P)*Aprime*P, 0, 0);
                this->H = ConjugateTranspose(P2)*this->H*P2;
                this->Q = this->Q*P2;
            }
            this->Qh = ConjugateTranspose(this->Q);
        }
    };

    /**
     * Struct for real Schur decomposition.
     * This struct NxN matrix Q and NxN matrix U such that \f$A=QUQ^*\f$ where Q is unitary and U is block upper triangular with 1x1 and 2x2 blocks.
     * The matrix U and A have the same eigenvalues.
     * @param T Type to store matrix entries as.
     * @param N Number of rows/columns for Q and U (both square). Dynamic is allowed for N.
     * @param Flags Flags to pass to the matrices (default = row major).
    */
    template <typename T, size_t N, unsigned int Flags = 0>
    struct Schur {
        SquareMatrix<T,N,Flags> Q; /*!< Unitary matrix */
        SquareMatrix<T,N,Flags> U; /*!< Block upper triangular with 1x1 and 2x2 blocks.*/
        SquareMatrix<T,N,Flags> Qh; /*!< Qh=ConjugateTranspose(Q) */

        /**
        * Constructor. Just calls Compute.
        * @param A M2xN2 Matrix
        */
        template <size_t M2, size_t N2, unsigned int Flags2>
        Schur(const Matrix<T,M2,N2,Flags2>& A) {
            Compute(A);
        }
        /**
        * Computes the Hessenberg decomposition. If A is not square or Q != N, then an exception is thrown.
        * @param A M2xN2 Matrix
        */
        template <size_t M2, size_t N2, unsigned int Flags2>
        void Compute(const Matrix<T,M2,N2,Flags2>& A, unsigned int max_iterations = 25) {
            if (!IsSquare(A))
                throw "Schur decomposition requires a square matrix.";
            if (A.NumColumns() != N && N != Dynamic)
                throw "Cannot perform Schur decomposition; size mismatch.";

            SquareMatrix<T,N,Flags> eye = Identity<T>(A.NumRows());
            this->Q = Identity<T>(A.NumRows());
            this->U = A;
            if (!IsHessenberg(A)) {
                Hessenberg<T,N,Flags> hess(A);
                this->Q = hess.Q;
                this->U = hess.H;
            }

            for (size_t i = A.NumRows()-1; i >= 1; i--) {
                size_t k = 0;
                while (k < max_iterations && Abs(this->U(i,i-1)) > T(Tol)) {
                    Complex<T> sigma = U(i,i);
                    QR<T,N,N,Flags> qr(this->U-sigma*eye);
                    this->U = qr.R*qr.Q + sigma*eye;
                    this->Q = this->Q*qr.Q;
                    k += 1;
                }
            }

            this->Qh = ConjugateTranspose(Q);
        }
    };
}
