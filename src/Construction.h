#pragma once
#include <vector>
#include <array>
#include <random>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"
#include "Global.h"

namespace Linear {
    template <typename T, size_t N, unsigned int Flags = 0>
    typename std::enable_if<(N>0), SquareMatrix<T,N,Flags>>::type Identity() {
        SquareMatrix<T,N,Flags> ret(T(0));
        for (size_t i = 0; i < N; ++i)
            ret(i,i) = T(1);
        return ret;
    }
    template<typename T, unsigned int Flags = 0>
    SquareMatrix<T,Dynamic,Flags> Identity(size_t n) {
        SquareMatrix<T,Dynamic,Flags> ret(n, n, T(0));
        for (size_t i = 0; i < n; ++i)
            ret(i,i) = T(1);
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    typename std::enable_if<(M>0&&N>0), Matrix<T,M,N,Flags>>::type Zero() {
        Matrix<T,M,N,Flags> ret(T(0));
        return ret;
    }
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> Zero(size_t nrows, size_t ncols) {
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, T(0));
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    typename std::enable_if<(M>0&&N>0), Matrix<T,M,N,Flags>>::type One() {
        Matrix<T,M,N,Flags> ret(T(1));
        return ret;
    }
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> One(size_t nrows, size_t ncols) {
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, T(1));
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    typename std::enable_if<(M>0&&N>0), Matrix<T,M,N,Flags>>::type Constant(Complex<T> v) {
        Matrix<T,M,N,Flags> ret(v);
        return ret;
    }
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> Constant(size_t nrows, size_t ncols, Complex<T> v) {
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, v);
        return ret;
    }
    template<typename T,size_t M,size_t N, unsigned int Flags = 0>
    typename std::enable_if<(M>0&&N>0), Matrix<T,M,N,Flags>>::type Basis(size_t i, size_t j) {
        Matrix<T,M,N,Flags> ret(T(0));
        ret(i,j) = T(1);
        return ret;
    }
    template<typename T, unsigned int Flags = 0>
    Matrix<T,Dynamic,Dynamic,Flags> Basis(size_t nrows, size_t ncols, size_t i, size_t j) {
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows, ncols, T(0));
        ret(i,j) = T(1);
        return ret;
    }
    std::mt19937 random_number_generator;
    void SeedRandom(unsigned int seed) {
        random_number_generator.seed(seed);
    }
    template <typename T, size_t M, size_t N, unsigned int Flags = 0>
    typename std::enable_if<(M>0&&N>0), Matrix<T,M,N,Flags>>::type Random(Complex<T> min = Complex<T>(0.0), Complex<T> max = Complex<T>(1.0)) {
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

    template <size_t P, typename T, size_t N, unsigned int Flags = 0>
    typename std::enable_if<(P>0), SquareMatrix<T,P,Flags>>::type Givens(size_t i, size_t j, Complex<T> theta) {
        if (i >= P || j >= P)
            throw "Cannot create Givens rotation, indices are out of bounds.";
        SquareMatrix<T,P,Flags> ret = Identity<T,P,Flags>();
        ret(i,i) = Cos(theta);
        ret(i,j) = Sin(theta);
        ret(j,i) = Sin(theta);
        ret(j,j) = Cos(theta);
        return ret;
    }
    template <typename T, size_t N, unsigned int Flags = 0>
    SquareMatrix<T,Dynamic,Flags> Givens(size_t size, size_t i, size_t j, Complex<T> theta) {
        if (i >= size || j >= size)
            throw "Cannot create Givens rotation, indices are out of bounds.";
        SquareMatrix<T,Dynamic,Flags> ret = Identity<T,Flags>(size);
        ret(i,i) = Cos(theta);
        ret(i,j) = Sin(theta);
        ret(j,i) = Sin(theta);
        ret(j,j) = Cos(theta);
        return ret;
    }

    template <typename T, size_t N, unsigned int Flags = 0>
    SquareMatrix<T,N,Flags> Householder(Vector<T,N> v) {
        T len = Length(v);
        if (len > T(Tol))
            v = v/len;
        SquareMatrix<T,N,Flags> ret = Identity<T,Flags>(v.Size()) - T(2)*v*ConjugateTranspose(v);
        return ret;
    }
    template <typename T, size_t N, unsigned int Flags = 0>
    SquareMatrix<T,N,Flags> Householder(const Vector<T,N>& x, size_t k) {
        if (k == 0 || k >= x.Size())
            throw "Cannot create Householder matrix, k must be such that 1<k<x.Size().";

        size_t index = x.Size()-k-1;
        Vector<T,Dynamic> x2 = SubVector(x,k+1,index);
        Complex<T> alpha = -Sign(x[index])*Length(x2);
        if (Abs(alpha) < T(Tol))
            alpha = Sqrt(2);
        Vector<T,N> v(x.Size(),T(0));
        for (size_t i = 0; i < index; ++i)
            v[i] = T(0);
        for (size_t i = 0; i < x2.Size(); ++i)
            v[index+i] = x2[i];
        v[index] = x2[0] - alpha;
        return Householder<T,Flags>(v);
    }


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
                ret(r+top.NumColumns(),c) = bottom(r,c);
            }
        }
        return ret;
    }

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
    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q>
    typename std::enable_if<((M==P||M==Dynamic||P==Dynamic)&&(N==Q||N==Dynamic||Q==Dynamic)),
        Matrix<T,(M==Dynamic?Dynamic:M+1),(N==Dynamic?Dynamic:N+1),Flags>>::type
    Block(const Matrix<T,M,N,Flags>& tl, const Vector<T,P>& tr, const RowVector<T,Q>& bl, Complex<T> br) {
        return Block(tl, tr, bl, Matrix<T,1,1,Flags>(br));
    }
    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q>
    typename std::enable_if<((M==P||M==Dynamic||P==Dynamic)&&(N==Q||N==Dynamic||Q==Dynamic)),
        Matrix<T,(M==Dynamic?Dynamic:M+1),(N==Dynamic?Dynamic:N+1),Flags>>::type
    Block(Complex<T> tl, const RowVector<T,Q>& tr, const Vector<T,P>& bl, const Matrix<T,M,N,Flags>& br) {
        return Block(Matrix<T,1,1,Flags>(tl), tr, bl, br);
    }

    template<typename T, size_t M, size_t N, unsigned int Flags1, size_t P, size_t Q, unsigned int Flags2>
    Matrix<T,M*P,Q*N,Flags1> Kronecker(const Matrix<T,M,N,Flags1>& a, const Matrix<T,P,Q,Flags2>& b) {
        Matrix<T,M*P,Q*N,Flags1> ret(a.NumRows()*b.NumRows(), a.NumColumns()*b.NumColumns(), T(0));
        for (size_t r = 0; r < ret.NumRows(); ++r) {
            for (size_t c = 0; c < ret.NumColumns(); ++c) {
                ret(r,c) = a(r/b.NumRows(),c/b.NumColumns())*b(r%b.NumRows(),c%b.NumColumns());
            }
        }
        return ret;
    }
    template<typename T, size_t M, size_t N, unsigned int Flags1, size_t P, size_t Q, unsigned int Flags2>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)&&(P==Q||P==Dynamic||Q==Dynamic)), Matrix<T,N*Q,N*Q,Flags1>>::type
    KroneckerSum(const Matrix<T,M,N,Flags1>& a, const Matrix<T,P,Q,Flags2>& b) {
        if (!IsSquare(a) || !IsSquare(b))
            throw "The Kronecker sum is only defined for two square matrices.";
        SquareMatrix<T,Dynamic,Flags1> In = Identity<T>(a.NumRows());
        SquareMatrix<T,Dynamic,Flags1> Iq = Identity<T>(b.NumRows());
        SquareMatrix<T,N*Q,Flags1> first = Kronecker(a, Iq);
        SquareMatrix<T,N*Q,Flags2> second = Kronecker(In, b);
        return first+second;
    }

    template<typename T, size_t M1, size_t N1, unsigned int Flags1, size_t M2, size_t N2, unsigned int Flags2>
    Matrix<T,(M1==Dynamic||M2==Dynamic?Dynamic:M1+M2),(N1==Dynamic||N2==Dynamic?Dynamic:N1+N2),Flags1> Diag(const Matrix<T,M1,N1,Flags1>& a, const Matrix<T,M2,N2,Flags2>& b) {
        Matrix<T,(M1==Dynamic||M2==Dynamic?Dynamic:M1+M2),(N1==Dynamic||N2==Dynamic?Dynamic:N1+N2),Flags1> ret(a.NumRows()+b.NumRows(),a.NumColumns()+b.NumColumns(),T(0));
        for (size_t r = 0; r < ret.NumRows(); ++r) {
            for (size_t c = 0; c < ret.NumColumns(); ++c) {
                if (r < a.NumRows() && c < a.NumColumns())
                    ret(r,c) = a(r,c);
                else if (r >= a.NumRows() && c >= a.NumColumns())
                    ret(r,c) = b(r-a.NumRows(),c-a.NumColumns());
                else
                    ret(r,c) = T(0);
            }
        }
        return ret;
    }
    template<typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,(M==Dynamic?Dynamic:M+1),(N==Dynamic?Dynamic:N+1),Flags> Diag(const Matrix<T,M,N,Flags>& a, Complex<T> b) {
        Matrix<T,(M==Dynamic?Dynamic:M+1),(N==Dynamic?Dynamic:N+1),Flags> ret(a.NumRows()+1,a.NumColumns()+1,T(0));
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                ret(r,c) = a(r,c);
            }
        }
        ret(ret.NumRows()-1,ret.NumColumns()-1) = b;
        return ret;
    }
    template<typename T, unsigned int Flags=0>
    Matrix<T,2,2,Flags> Diag(Complex<T> a, Complex<T> b) {
        Matrix<T,2,2,Flags> ret(T(0));
        ret(0,0) = a;
        ret(1,1) = b;
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
    template<typename T>
    SquareMatrix<T,Dynamic> Diag(std::vector<Complex<T>> v) {
        if (v.size() == 0)
            throw "Cannot create a diagonal matrix based on a empty list.";
        SquareMatrix<T,Dynamic> ret(v.size(), v.size(), T(0));
        for (size_t i = 0; i < v.size(); ++i)
            ret(i,i) = v[i];
        return ret;
    }
    template<typename T, unsigned int Flags>
    Matrix<T,Dynamic,Dynamic,Flags> Diag(std::vector<Matrix<T,Dynamic,Dynamic,Flags>> v) {
        if (v.size() == 0)
            throw "Cannot create a diagonal matrix based on a empty list.";
        Matrix<T,Dynamic,Dynamic,Flags> ret = v[0];
        for (size_t i = 0; i < v.size(); ++i)
            ret = Diag(ret, v[i]);
        return ret;
    }
    template<typename T, size_t N>
    SquareMatrix<T,N> Diag(Vector<T,N> v) {
        SquareMatrix<T,N> ret(v.NumRows(), v.NumRows(), T(0));
        for (size_t i = 0; i < v.NumRows(); ++i)
            ret(i,i) = v(i,0);
        return ret;
    }
    template<typename T, size_t N>
    SquareMatrix<T,N> Diag(RowVector<T,N> v) {
        return Diag(Transpose(v));
    }
}
