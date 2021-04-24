#pragma once
#include <vector>
#include <array>
#include <random>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"
#include "Global.h"

namespace Linear {
    /**
     * Creates the NxN Identity matrix.
     *
     * Example: I is the 3x3 matrix {{1,0,0}, {0,1,0}, {0,0,1}}
     *
     *     Matrix3d I = Identity<double, 3>();
     *
     * @param T Type
     * @param N Size of matrix.
     * @param Flags Flags to pass to the matrix (default = row major)
     * @return NxN Matrix
     */
    template <typename T, size_t N, unsigned int Flags = 0>
    SquareMatrix<T,N,Flags> Identity() {
        SquareMatrix<T,N,Flags> ret(T(0));
        for (size_t i = 0; i < N; ++i)
            ret(i,i) = T(1);
        return ret;
    }
    /**
     * Dynamically creates the nxn Identity matrix.
     *
     * Example: I is the 3x3 matrix {{1,0,0}, {0,1,0}, {0,0,1}}
     *
     *     MatrixXd I = Identity<double>(3);
     *
     * @param T Type
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param n Size of matrix
     * @return nxn Matrix
     */
    template<typename T, unsigned int Flags = 0>
    SquareMatrix<T,Dynamic,Flags> Identity(size_t n) {
        SquareMatrix<T,Dynamic,Flags> ret(n, n, T(0));
        for (size_t i = 0; i < n; ++i)
            ret(i,i) = T(1);
        return ret;
    }

    /**
     * Creates the MxN all zeros matrix.
     *
     * Example: A is the 2x3 matrix {{0,0,0}, {0,0,0}}
     *
     *     Matrix<double,2,3> A = Zero<double, 2, 3>();
     *
     * @param T Type
     * @param M Number of rows
     * @param N Number of columns
     * @param Flags Flags to pass to the matrix (default = row major)
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    Matrix<T,M,N,Flags> Zero() {
        Matrix<T,M,N,Flags> ret(T(0));
        return ret;
    }
    /**
     * Dynamically creates the nrowsxncols all zeros matrix.
     *
     * Example: A is the 2x3 matrix {{0,0,0}, {0,0,0}}
     *
     *     MatrixXd A = Zero<double>(2, 3);
     *
     * @param T Type
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @return nrowsxncols Matrix
     */
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> Zero(size_t nrows, size_t ncols) {
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, T(0));
        return ret;
    }

    /**
     * Creates the MxN all ones matrix.
     *
     * Example: A is the 2x3 matrix {{1,1,1}, {1,1,1}}
     *
     *     Matrix<double,2,3> A = One<double, 2, 3>();
     *
     * @param T Type
     * @param M Number of rows
     * @param N Number of columns
     * @param Flags Flags to pass to the matrix (default = row major)
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    Matrix<T,M,N,Flags> One() {
        Matrix<T,M,N,Flags> ret(T(1));
        return ret;
    }
    /**
     * Dynamically creates the nrowsxncols all ones matrix.
     *
     * Example: A is the 2x3 matrix {{1,1,1}, {1,1,1}}
     *
     *     MatrixXd A = One<double>(2, 3);
     *
     * @param T Type
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @return nrowsxncols Matrix
     */
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> One(size_t nrows, size_t ncols) {
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, T(1));
        return ret;
    }

    /**
     * Creates an MxN matrix with every entry set to z.
     *
     * Example: A is the 2x3 matrix {{1+1i,1+1i,1+1i}, {1+1i,1+1i,1+1i}}
     *
     *     Matrix<double,2,3> A = Constant<double, 2, 3>(Complexd(1,1));
     *
     * @param T Type
     * @param M Number of rows
     * @param N Number of columns
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param z Complex number
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    Matrix<T,M,N,Flags> Constant(Complex<T> z) {
        Matrix<T,M,N,Flags> ret(z);
        return ret;
    }
    /**
     * Dynamically creates an nrowsxncols matrix with every entry set to v.
     *
     * Example: A is the 2x3 matrix {{1+1i,1+1i,1+1i}, {1+1i,1+1i,1+1i}}
     *
     *     MatrixXd A = Constant<double>(2, 3, Complexd(1,1));
     *
     * @param T Type
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param z Complex number
     * @return nrowsxncols Matrix
     */
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> Constant(size_t nrows, size_t ncols, Complex<T> z) {
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, z);
        return ret;
    }

    /**
     * Creates an MxN basis matrix \f$A\f$ where all entries are set to 0 except \f$a_{ij}=1\f$. If i or j are
     * out of bounds, an exception is thrown.
     *
     * Example: A is the 2x2 matrix {{0,1}, {0,0}}
     *
     *     Matrix2d A = Basis<double, 2, 2>(0, 1);
     *
     * @param T Type
     * @param M Number of rows
     * @param N Number of columns
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param i Row index
     * @param j Column index
     * @return MxN Matrix
     */
    template<typename T,size_t M,size_t N, unsigned int Flags = 0>
    typename std::enable_if<(M>0&&N>0), Matrix<T,M,N,Flags>>::type Basis(size_t i, size_t j) {
        if (i >= M && j >= N)
            throw "Cannot create basis matrix, indices out of bounds.";
        Matrix<T,M,N,Flags> ret(T(0));
        ret(i,j) = T(1);
        return ret;
    }
    /**
     * Dynamically creates an nrowsxncols basis matrix \f$A\f$ where all entries are set to 0 except \f$a_{ij}=1\f$. If i or j are
     * out of bounds, an exception is thrown.
     *
     * Example: A is the 2x2 matrix {{0,1}, {0,0}}
     *
     *     MatrixXd A = Basis<double>(2, 2, 0, 1);
     *
     * @param T Type
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param i Row index
     * @param j Column index
     * @return nrowsxncols Matrix
     */
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> Basis(size_t nrows, size_t ncols, size_t i, size_t j) {
        if (i >= nrows && j >= ncols)
            throw "Cannot create basis matrix, indices out of bounds.";
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, T(0));
        ret(i,j) = T(1);
        return ret;
    }

    std::mt19937 random_number_generator;
    /**
     * Seeds the random number generator used in Random.
     * @param seed Seed to feed random_number_generator
     */
    void SeedRandom(unsigned int seed) {
        random_number_generator.seed(seed);
    }
    /**
     * Creates an MxN random matrix \f$A\f$ where all entries are randomly determined between min and max. If \f$min=a+bi\f$ and \f$max=c+di\f$,
     * then each entry of \f$A\f$ will be of the form \f$x+iy\f$ where \f$a\le x\le c\f$ and \f$b\le y\le d\f$.
     * @param T Type
     * @param M Number of rows
     * @param N Number of columns
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param min Complex number (default = 0)
     * @param max Complex number (default = 1+0i)
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    Matrix<T,M,N,Flags> Random(Complex<T> min = Complex<T>(0.0), Complex<T> max = Complex<T>(1.0)) {
        std::uniform_real_distribution<T> distrRe(min.Re, max.Re);
        std::uniform_real_distribution<T> distrIm(min.Im, max.Im);
        Matrix<T,M,N,Flags> ret(T(0));
        for (size_t r = 0; r < ret.NumRows(); ++r) {
            for (size_t c = 0; c < ret.NumColumns(); ++c) {
                ret(r,c) = Complex<T>(distrRe(random_number_generator), distrIm(random_number_generator));
            }
        }
        return ret;
    }
    /**
     * Dynamically creates an nrowsxncols random matrix \f$A\f$ where all entries are randomly determined between min and max. If \f$min=a+bi\f$ and \f$max=c+di\f$,
     * then each entry of \f$A\f$ will be of the form \f$x+iy\f$ where \f$a\le x\le c\f$ and \f$b\le y\le d\f$.
     * @param T Type
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param min Complex number (default = 0)
     * @param max Complex number (default = 1+0i)
     * @return nrowsxncols Matrix
     */
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> Random(size_t nrows, size_t ncols, Complex<T> min = Complex<T>(0.0), Complex<T> max = Complex<T>(1.0)) {
        std::uniform_real_distribution<T> distrRe(min.Re, max.Re);
        std::uniform_real_distribution<T> distrIm(min.Im, max.Im);
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, T(0));
        for (size_t r = 0; r < ret.NumRows(); ++r) {
            for (size_t c = 0; c < ret.NumColumns(); ++c) {
                ret(r,c) = Complex<T>(distrRe(random_number_generator), distrIm(random_number_generator));
            }
        }
        return ret;
    }

    /**
     * Creates an NxN Householder transformation \f$H\f$ defined by \f$H = I - 2vv^*.\f$
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param v Column vector of length N
     * @return NxN Matrix
     */
    template <typename T, size_t N, unsigned int Flags = 0>
    SquareMatrix<T,N,Flags> Householder(Vector<T,N> v) {
        v = Normalize(v);
        SquareMatrix<T,N,Flags> ret = Identity<T,Flags>(v.Length()) - T(2)*v*ConjugateTranspose(v);
        return ret;
    }
    /**
     * Creates an NxN Householder transformation \f$H\f$ such that \f$Hx\f$ eliminates the last k elements.
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param x Column vector of length N
     * @param k Number of elements to zero out
     * @return NxN Matrix
     */
    template <typename T, size_t N, unsigned int Flags = 0>
    SquareMatrix<T,N,Flags> Householder(const Vector<T,N>& x, size_t k) {
        if (k == 0 || k >= x.Length())
            throw "Cannot create Householder matrix, k must be such that 1<k<x.Length().";

        size_t index = x.Length()-k-1;
        Vector<T,Dynamic> x2 = SubVector(x,k+1,index);
        Complex<T> alpha = -Sign(x[index])*Norm(x2);
        if (Abs(alpha) < T(Tol))
            alpha = Sqrt(2);
        Vector<T,N> v(x.Length(),T(0));
        for (size_t i = 0; i < index; ++i)
            v[i] = T(0);
        for (size_t i = 0; i < x2.Length(); ++i)
            v[index+i] = x2[i];
        v[index] = x2[0] - alpha;
        return Householder<T,Flags>(v);
    }


    /**
     * Creates the M1x(N1+N2) augmented matrix
     * \f\[ \begin{bmatrix} left & \mid & right \end{bmatrix}. \f\]
     * If the two matrices have differing numbers of rows, an exception is thrown.
     * @param left M1xN1 Matrix
     * @param right M2xN2 Matrix
     * @return M1x(N1+N2) Matrix
     */
    template<typename T, size_t M1, size_t N1, unsigned int Flags1, size_t M2, size_t N2, unsigned int Flags2>
    typename std::enable_if<(M1==M2||M1==Dynamic||M2==Dynamic), Matrix<T,M1,(N1==Dynamic||N2==Dynamic?Dynamic:N1+N2),Flags1>>::type
    Augmented(const Matrix<T,M1,N1,Flags1>& left, const Matrix<T,M2,N2,Flags2>& right) {
        if (left.NumRows() != right.NumRows())
            throw "Malformed augmented matrix. Differing number of rows.";
        Matrix<T,M1,(N1==Dynamic||N2==Dynamic?Dynamic:N1+N2),Flags1> ret(left.NumRows(),left.NumColumns()+right.NumColumns(),T(0));
        for (size_t r = 0; r < left.NumRows(); ++r) {
            for (size_t c = 0; c < left.NumColumns(); ++c) {
                ret(r,c) = left(r,c);
            }
        }
        for (size_t r = 0; r < right.NumRows(); ++r) {
            for (size_t c = 0; c < right.NumColumns(); ++c) {
                ret(r,c+left.NumColumns()) = right(r,c);
            }
        }
        return ret;
    }
    /**
     * Creates the (M1+M2)xN1 row augmented matrix
     * \f\[ \begin{bmatrix} top \\ \hline bottom \end{bmatrix}. \f\]
     * If the two matrices have differing numbers of columns, an exception is thrown.
     * @param top M1xN1 Matrix
     * @param bottom M2xN2 Matrix
     * @return (M1+M2)xN1 Matrix
     */
    template<typename T, size_t M1, size_t N1, unsigned int Flags1, size_t M2, size_t N2, unsigned int Flags2>
    typename std::enable_if<(N1==N2||N1==Dynamic||N2==Dynamic), Matrix<T,(M1==Dynamic||M2==Dynamic?Dynamic:M1+M2),N1,Flags1>>::type
    RowAugmented(const Matrix<T,M1,N1,Flags1>& top, const Matrix<T,M2,N2,Flags2>& bottom) {
        if (top.NumColumns() != bottom.NumColumns())
            throw "Malformed row augmented matrix. Differing number of columns.";
        Matrix<T,(M1==Dynamic||M2==Dynamic?Dynamic:M1+M2),N1,Flags1> ret(top.NumRows()+bottom.NumRows(),top.NumColumns(),T(0));
        for (size_t r = 0; r < top.NumRows(); ++r) {
            for (size_t c = 0; c < top.NumColumns(); ++c) {
                ret(r,c) = top(r,c);
            }
        }
        for (size_t r = 0; r < bottom.NumRows(); ++r) {
            for (size_t c = 0; c < bottom.NumColumns(); ++c) {
                ret(r+top.NumRows(),c) = bottom(r,c);
            }
        }
        return ret;
    }

    /**
     * Creates the (M1+M3)x(N1+N2) block matrix
     * \f\[ \begin{bmatrix} tl & tr \\ bl & br \end{bmatrix}. \f\]
     * If \f$M1+M3\ne M2+M4\f$ or \f$N1+N2\ne N3+N4\f$, then an exception is thrown.
     * @param tl M1xN1 Matrix
     * @param tr M2xN2 Matrix
     * @param bl M3xN3 Matrix
     * @param br M4xN4 Matrix
     * @return (M1+M3)x(N1+N2) Matrix
     */
    template<typename T, size_t M1, size_t N1, unsigned int Flags1,
            size_t M2, size_t N2, unsigned int Flags2,
            size_t M3, size_t N3, unsigned int Flags3,
            size_t M4, size_t N4, unsigned int Flags4>
    typename std::enable_if<((M1+M3==M2+M4||M1==Dynamic||M2==Dynamic||M3==Dynamic||M4==Dynamic)&&(N1+N2==N3+N4||N1==Dynamic||N2==Dynamic||N3==Dynamic||N4==Dynamic)),
        Matrix<T,(M1==Dynamic||M2==Dynamic||M3==Dynamic||M4==Dynamic?Dynamic:M1+M3),(N1==Dynamic||N2==Dynamic||N3==Dynamic||N4==Dynamic?Dynamic:N1+N2),Flags1>>::type
    Block(const Matrix<T,M1,N1,Flags1>& tl, const Matrix<T,M2,N2,Flags2>& tr, const Matrix<T,M3,N3,Flags3>& bl, const Matrix<T,M4,N4,Flags4>& br) {
        if (tl.NumRows()+bl.NumRows() != tr.NumRows()+br.NumRows() || tl.NumColumns()+tr.NumColumns() != bl.NumColumns()+br.NumColumns())
            throw "Malformed block matrix. Sizes don't add up equally.";
        Matrix<T,(M1==Dynamic||M2==Dynamic||M3==Dynamic||M4==Dynamic?Dynamic:M1+M3),(N1==Dynamic||N2==Dynamic||N3==Dynamic||N4==Dynamic?Dynamic:N1+N2),Flags1> ret(tl.NumRows()+bl.NumRows(),tl.NumColumns()+tr.NumColumns());
        for (size_t r = 0; r < tl.NumRows(); ++r) {
            for (size_t c = 0; c < tl.NumColumns(); ++c) {
                ret(r,c) = tl(r,c);
            }
        }

        for (size_t r = 0; r < tr.NumRows(); ++r) {
            for (size_t c = 0; c < tr.NumColumns(); ++c) {
                ret(r,c+tl.NumColumns()) = tr(r,c);
            }
        }

        for (size_t r = 0; r < bl.NumRows(); ++r) {
            for (size_t c = 0; c < bl.NumColumns(); ++c) {
                ret(r+tl.NumRows(),c) = bl(r,c);
            }
        }

        for (size_t r = 0; r < br.NumRows(); ++r) {
            for (size_t c = 0; c < br.NumColumns(); ++c) {
                ret(r+tl.NumRows(),c+bl.NumColumns()) = br(r,c);
            }
        }
        return ret;
    }
    /**
     * Creates the (M+1)x(N+1) block matrix
     * \f\[ \begin{bmatrix} tl & tr \\ bl & br \end{bmatrix}. \f\]
     * If \f$P\ne M\f$ or \f$Q\ne N\f$, then an exception is thrown.
     * @param tl MxN Matrix
     * @param tr Vector of length P
     * @param bl Row vector of length Q
     * @param br Complex number
     * @return (M+1)x(N+1) Matrix
     */
    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q>
    typename std::enable_if<((M==P||M==Dynamic||P==Dynamic)&&(N==Q||N==Dynamic||Q==Dynamic)),
        Matrix<T,(M==Dynamic?Dynamic:M+1),(N==Dynamic?Dynamic:N+1),Flags>>::type
    Block(const Matrix<T,M,N,Flags>& tl, const Vector<T,P>& tr, const RowVector<T,Q>& bl, Complex<T> br) {
        return Block(tl, tr, bl, Matrix<T,1,1,Flags>(br));
    }
    /**
     * Creates the (M+1)x(N+1) block matrix
     * \f\[ \begin{bmatrix} tl & tr \\ bl & br \end{bmatrix}. \f\]
     * If \f$P\ne M\f$ or \f$Q\ne N\f$, then an exception is thrown.
     * @param tl Complex number
     * @param tr Row vector of length Q
     * @param bl Vector of length P
     * @param br MxN Matrix
     * @return (M+1)x(N+1) Matrix
     */
    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q>
    typename std::enable_if<((M==P||M==Dynamic||P==Dynamic)&&(N==Q||N==Dynamic||Q==Dynamic)),
        Matrix<T,(M==Dynamic?Dynamic:M+1),(N==Dynamic?Dynamic:N+1),Flags>>::type
    Block(Complex<T> tl, const RowVector<T,Q>& tr, const Vector<T,P>& bl, const Matrix<T,M,N,Flags>& br) {
        return Block(Matrix<T,1,1,Flags>(tl), tr, bl, br);
    }

    /**
     * Computes the (MP)x(NQ) block matrix
     * \f\[
     *    A\otimes B = \begin{bmatrix}
     *    a_{0,0}B & \dots & a_{0,N-1}B \\
     *    \vdots & \ddots & \vdots \\
     *    a_{M-1,0}B & \dots & a_{M-1,N-1}B
     *    \end{bmatrix}.
     * \f\]
     * @param A MxN Matrix
     * @param B PxQ Matrix
     * @return (MP)x(NQ) Matrix
     */
    template<typename T, size_t M, size_t N, unsigned int Flags1, size_t P, size_t Q, unsigned int Flags2>
    Matrix<T,M*P,Q*N,Flags1> Kronecker(const Matrix<T,M,N,Flags1>& A, const Matrix<T,P,Q,Flags2>& B) {
        Matrix<T,M*P,Q*N,Flags1> ret(A.NumRows()*B.NumRows(), A.NumColumns()*B.NumColumns(), T(0));
        for (size_t r = 0; r < ret.NumRows(); ++r) {
            for (size_t c = 0; c < ret.NumColumns(); ++c) {
                ret(r,c) = A(r/B.NumRows(),c/B.NumColumns())*B(r%B.NumRows(),c%B.NumColumns());
            }
        }
        return ret;
    }
    /**
     * Computes the (NQ)x(NQ) Kronecker sum defined by
     * \f\[
     *    A \oplus B = A\otimes I_Q + I_N\otimes B
     * \f\]
     * where \f$I_N\f$ represents the NxN identity matrix.
     * @param A NxN Matrix
     * @param B QxQ Matrix
     * @return (NQ)x(NQ) Matrix
     */
    template<typename T, size_t M, size_t N, unsigned int Flags1, size_t P, size_t Q, unsigned int Flags2>
    SquareMatrix<T,N*Q,Flags1> KroneckerSum(const Matrix<T,M,N,Flags1>& A, const Matrix<T,P,Q,Flags2>& B) {
        if (!IsSquare(A) || !IsSquare(B))
            throw "The Kronecker sum is only defined for two square matrices.";
        if (A.NumEntries() == 0 || B.NumEntries() == 0) {
            SquareMatrix<T,N*Q,Flags1> ret(T(0));
            return ret;
        }
        SquareMatrix<T,Dynamic,Flags1> In = Identity<T,Flags1>(A.NumRows());
        SquareMatrix<T,Dynamic,Flags1> Iq = Identity<T,Flags1>(B.NumRows());
        SquareMatrix<T,N*Q,Flags1> first = Kronecker(A, Iq);
        SquareMatrix<T,N*Q,Flags2> second = Kronecker(In, B);
        return first+second;
    }

    /**
     * Computes the (M1+M2)x(N1+N2) block diagonal matrix
     * \f\[
     *    \begin{bmatrix} A & 0 \\ 0 & B \end{bmatrix}.
     * \f\]
     * @param A M1xN1 Matrix
     * @param B M2xN2 Matrix
     * @return (M1+M2)x(N1+N2) Matrix
     */
    template<typename T, size_t M1, size_t N1, unsigned int Flags1, size_t M2, size_t N2, unsigned int Flags2>
    Matrix<T,(M1==Dynamic||M2==Dynamic?Dynamic:M1+M2),(N1==Dynamic||N2==Dynamic?Dynamic:N1+N2),Flags1> Diag(const Matrix<T,M1,N1,Flags1>& A, const Matrix<T,M2,N2,Flags2>& B) {
        Matrix<T,(M1==Dynamic||M2==Dynamic?Dynamic:M1+M2),(N1==Dynamic||N2==Dynamic?Dynamic:N1+N2),Flags1> ret(A.NumRows()+B.NumRows(),A.NumColumns()+B.NumColumns(),T(0));
        for (size_t r = 0; r < ret.NumRows(); ++r) {
            for (size_t c = 0; c < ret.NumColumns(); ++c) {
                if (r < A.NumRows() && c < A.NumColumns())
                    ret(r,c) = A(r,c);
                else if (r >= A.NumRows() && c >= A.NumColumns())
                    ret(r,c) = B(r-A.NumRows(),c-A.NumColumns());
                else
                    ret(r,c) = T(0);
            }
        }
        return ret;
    }
    /**
     * Computes the (M+1)x(N+1) block diagonal matrix
     * \f\[
     *    \begin{bmatrix} A & 0 \\ 0 & z \end{bmatrix}.
     * \f\]
     * @param A MxN Matrix
     * @param z Complex number
     * @return (M+1)x(N+1) Matrix
     */
    template<typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,(M==Dynamic?Dynamic:M+1),(N==Dynamic?Dynamic:N+1),Flags> Diag(const Matrix<T,M,N,Flags>& A, Complex<T> z) {
        Matrix<T,(M==Dynamic?Dynamic:M+1),(N==Dynamic?Dynamic:N+1),Flags> ret(A.NumRows()+1,A.NumColumns()+1,T(0));
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                ret(r,c) = A(r,c);
            }
        }
        ret(ret.NumRows()-1,ret.NumColumns()-1) = z;
        return ret;
    }
    /**
     * Computes the 2x2 diagonal matrix
     * \f\[
     *    \begin{bmatrix} z & 0 \\ 0 & w \end{bmatrix}.
     * \f\]
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param z Complex number
     * @param w Complex number
     * @return 2x2 Matrix
     */
    template<typename T, unsigned int Flags=0>
    Matrix<T,2,2,Flags> Diag(Complex<T> z, Complex<T> w) {
        Matrix<T,2,2,Flags> ret(T(0));
        ret(0,0) = z;
        ret(1,1) = w;
        return ret;
    }
    template<typename T, size_t M1, size_t N1, unsigned int Flags1, size_t M2, size_t N2, unsigned int Flags2, typename... Ts>
    Matrix<T,Dynamic,Dynamic,Flags1> Diag(const Matrix<T,M1,N1,Flags1>& a, const Matrix<T,M2,N2,Flags2>& b, Ts&&... ts) {
        Matrix<T,Dynamic,Dynamic,Flags1> m = Diag(a,b);
        return Diag(m, Matrix<T,Dynamic,Dynamic,Flags1>(std::forward<Ts>(ts))...);
    }
    template<typename T, typename... Ts>
    Matrix<T,Dynamic,Dynamic> Diag(Complex<T> a, Complex<T> b, Ts&&... ts) {
        Matrix<T,Dynamic,Dynamic> m = Diag(a,b);
        return Diag(m, (Complex<T>)std::forward<Ts>(ts)...);
    }
    /**
     * Dynamically computes the NxN diagonal matrix
     * \f\[
     *    diag(v[0],\dots,v[N-1])=\begin{bmatrix} v[0] & & \\ & \ddots & \\ & & v[N-1] \end{bmatrix}.
     * \f\]
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param v List of N Complex numbers
     * @return NxN Matrix
     */
    template<typename T, unsigned int Flags=0>
    SquareMatrix<T,Dynamic> Diag(std::vector<Complex<T>> v) {
        if (v.size() == 0)
            throw "Cannot create a diagonal matrix based on a empty list.";
        SquareMatrix<T,Dynamic,Flags> ret(v.size(), v.size(), T(0));
        for (size_t i = 0; i < v.size(); ++i)
            ret(i,i) = v[i];
        return ret;
    }
    /**
     * Computes the PxQ diagonal matrix, where \f$P=\sum_{i=0}^{N-1}NumRows(As[i])\f$ and \f$Q=\sum_{i=0}^{N-1}NumColumns(As[i])\f$,
     * \f\[
     *    diag(As[0],\dots,As[N-1])=\begin{bmatrix} As[0] & & \\ & \ddots & \\ & & As[N-1] \end{bmatrix}.
     * \f\]
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param As List of N matrices
     * @return PxQ Matrix
     */
    template<typename T, unsigned int Flags>
    Matrix<T,Dynamic,Dynamic,Flags> Diag(std::vector<Matrix<T,Dynamic,Dynamic,Flags>> As) {
        if (As.size() == 0)
            throw "Cannot create a diagonal matrix based on a empty list.";
        Matrix<T,Dynamic,Dynamic,Flags> ret = As[0];
        for (size_t i = 0; i < As.size(); ++i)
            ret = Diag(ret, As[i]);
        return ret;
    }
    /**
     * Computes the NxN diagonal matrix
     * \f\[
     *    diag(v_0,\dots,v_{N-1}) = \begin{bmatrix} v_0 & & \\ & \ddots & \\ & & v_{N-1} \end{bmatrix}.
     * \f\]
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param v Vector of size N
     * @return NxN Matrix
     */
    template<typename T, size_t N, unsigned int Flags=0>
    SquareMatrix<T,N> Diag(Vector<T,N> v) {
        SquareMatrix<T,N,Flags> ret(v.NumRows(), v.NumRows(), T(0));
        for (size_t i = 0; i < v.Length(); ++i)
            ret(i,i) = v[i];
        return ret;
    }
    /**
     * Computes the NxN diagonal matrix
     * \f\[
     *    diag(v_0,\dots,v_{N-1})=\begin{bmatrix} v_0 & & \\ & \ddots & \\ & & v_{N-1} \end{bmatrix}.
     * \f\]
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param v Row vector of size N
     * @return NxN Matrix
     */
    template<typename T, size_t N, unsigned int Flags=0>
    SquareMatrix<T,N> Diag(RowVector<T,N> v) {
        return Diag<T,N,Flags>(Transpose(v));
    }

    /**
     * Computes the (N-1)x(N-1) companion matrix
     * \f\[
     *    C = \begin{bmatrix}
     *      0 & 0 & \dots & 0 & -c_0/c_{N-1} \\
     *      1 & 0 & \dots & 0 & -c_1/c_{N-1}  \\
     *      0 & 1 & \dots & 0 & -c_2/c_{N-1} \\
     *      \vdots & \vdots & \ddots & \vdots & \vdots \\
     *      0 & 0 & \dots & 1 & -c_{N-2}/c_{N-1}
     *    \end{bmatrix}
     * \f\]
     * such that \f$\det(C-tI)=c_0/c_{N-1}+c_1/c_{N-1}t+\dots+c_{N-2}/c_{N-1}t^{N-2}+t^{N-1}\f$.
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param c Vector of size N
     * @return (N-1)x(N-1) Matrix
     */
    template <typename T, size_t N, unsigned int Flags=0>
    typename std::enable_if<(N>1||N==Dynamic), SquareMatrix<T,(N==Dynamic?Dynamic:N-1),Flags>>::type Companion(const Vector<T,N>& c) {
        if (c.Length() <= 1)
            throw "Cannot create a 0x0 companion matrix.";
        if (Abs(c[c.Length()-1]) < T(Tol))
            throw "Leading coefficient cannot be zero.";
        // Create the matrix.
        SquareMatrix<T,(N==Dynamic?Dynamic:N-1),Flags> ret(c.Length()-1,T(0));
        for (size_t i = 0; i < ret.NumColumns()-1; ++i)
            ret(i+1,i) = T(1);
        for (size_t i = 0; i < ret.NumColumns(); ++i)
            ret(i,ret.NumColumns()-1) = -c[i]/c[c.Length()-1];
        return ret;
    }
    /**
     * Computes the (N-1)x(N-1) companion matrix
     * \f\[
     *    C = \begin{bmatrix}
     *      0 & 0 & \dots & 0 & -c_0/c_{N-1} \\
     *      1 & 0 & \dots & 0 & -c_1/c_{N-1}  \\
     *      0 & 1 & \dots & 0 & -c_2/c_{N-1} \\
     *      \vdots & \vdots & \ddots & \vdots & \vdots \\
     *      0 & 0 & \dots & 1 & -c_{N-2}/c_{N-1}
     *    \end{bmatrix}
     * \f\]
     * such that \f$\det(C-tI)=c_0/c_{N-1}+c_1/c_{N-1}t+\dots+c_{N-2}/c_{N-1}t^{N-2}+t^{N-1}\f$.
     * @param Flags Flags to pass to the matrix (default = row major)
     * @param c Row vector of size N
     * @return (N-1)x(N-1) Matrix
     */
    template <typename T, size_t N, unsigned int Flags=0>
    typename std::enable_if<(N>1||N==Dynamic), SquareMatrix<T,(N==Dynamic?Dynamic:N-1),Flags>>::type Companion(const RowVector<T,N>& c) {
        return Companion<Flags>(Transpose(c));
    }
}
