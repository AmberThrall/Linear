#pragma once
#include "Matrix.h"
#include "Vector.h"
#include "Eigen.h"
#include "Global.h"
#include "Types.h"

namespace Linear {
    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    typename std::enable_if<((M==P||M==Dynamic||P==Dynamic)&&(N==Q||N==Dynamic||Q==Dynamic)), Matrix<T,M,N,Flags>>::type
    EntrywiseProduct(const Matrix<T,M,N,Flags>& a, const Matrix<T,P,Q,Flags2>& b) {
        if (a.NumRows() != b.NumRows() || a.NumColumns() != b.NumColumns())
            throw "Cannot perform entrywise product when matrices have varying sizes.";

        Matrix<T,M,N,Flags> ret(a.NumRows(),a.NumColumns(),T(0));
        for (size_t i = 0; i < a.NumRows(); ++i) {
            for (size_t j = 0; j < a.NumColumns(); ++j)
                ret(i,j) = a(i,j)*b(i,j);
        }
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    typename std::enable_if<((M==P||M==Dynamic||P==Dynamic)&&(N==Q||N==Dynamic||Q==Dynamic)), Matrix<T,M,N,Flags>>::type
    EntrywiseDivision(const Matrix<T,M,N,Flags>& a, const Matrix<T,P,Q,Flags2>& b) {
        if (a.NumRows() != b.NumRows() || a.NumColumns() != b.NumColumns())
            throw "Cannot perform entrywise division when matrices have varying sizes.";

        Matrix<T,M,N,Flags> ret(a.NumRows(),a.NumColumns(),T(0));
        for (size_t i = 0; i < a.NumRows(); ++i) {
            for (size_t j = 0; j < a.NumColumns(); ++j)
                ret(i,j) = a(i,j)/b(i,j);
        }
        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    T Norm(const Matrix<T,M,N,Flags>& a, T p = T(2)) {
        if (IsVector(a)) {
            T outerSum = T(0);
            for (size_t j = 0; j < a.NumColumns(); ++j) {
                for (size_t i = 0; i < a.NumRows(); ++i)
                    outerSum += Pow(Abs(a(i,j)), p);
            }
            return Pow(outerSum, 1/p);
        }
        if (p != T(1) && p != T(2))
            throw "Only p=1 and p=2 are supported for matrix norm.";
        if (p == T(1)) {
            T ret = T(0);
            for (size_t j = 0; j < a.NumColumns(); ++j) {
                T sum = T(0);
                for (size_t i = 0; i < a.NumRows(); ++i) {
                    sum += Abs(a(i,j));
                }
                if (sum > ret)
                    ret = sum;
            }
            return ret;
        }
        else if (p == T(2)) {
            Vector<T,N> b = Random<T,N,1>(a.NumRows(),1);
            std::pair<Complex<T>,Vector<T,N>> pair = PowerIteration(ConjugateTranspose(a)*a, b, 25);
            return Sqrt(pair.first.Re);
        }
        else
            throw "Only p=1 and p=2 are supported for matrix norm.";
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T EntrywiseNorm(const Matrix<T,M,N,Flags>& a, T p = T(2)) {
        T outerSum = T(0);
        for (size_t j = 0; j < a.NumColumns(); ++j) {
            for (size_t i = 0; i < a.NumRows(); ++i) {
                outerSum += Pow(Abs(a(i,j)), T(p));
            }
        }
        return Pow(outerSum, 1/T(p));
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T FrobeniusNorm(const Matrix<T,M,N,Flags>& a) {
        return EntrywiseNorm(a, 2,  2);
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T MaxNorm(const Matrix<T,M,N,Flags>& a) {
        T ret = T(0);
        for (size_t j = 0; j < a.NumColumns(); ++j) {
            for (size_t i = 0; i < a.NumRows(); ++i) {
                if (Abs(a(i,j)) > ret)
                    ret = Abs(a(i,j));
            }
        }
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T InfinityNorm(const Matrix<T,M,N,Flags>& a) {
        T ret = T(0);
        for (size_t i = 0; i < a.NumRows(); ++i) {
            T sum = T(0);
            for (size_t j = 0; j < a.NumColumns(); ++j) {
                sum += Abs(a(i,j));
            }
            if (sum > ret)
                ret = sum;
        }
        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsReal(const Matrix<T,M,N,Flags>& a) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                if (!IsReal(a(r,c)))
                    return false;
            }
        }
        return true;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Abs(Matrix<T,M,N,Flags> a) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Abs(a(r,c));
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Sign(Matrix<T,M,N,Flags> a) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Sign(a(r,c));
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Arg(Matrix<T,M,N,Flags> a) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Arg(a(r,c));
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Conjugate(Matrix<T,M,N,Flags> a) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Conjugate(a(r,c));
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Sqrt(Matrix<T,M,N,Flags> a) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Sqrt(a(r,c));
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Exp(Matrix<T,M,N,Flags> a) {
        if (!IsSquare(a)) {
            for (size_t r = 0; r < a.NumRows(); ++r) {
                for (size_t c = 0; c < a.NumColumns(); ++c) {
                    a(r,c) = Exp(a(r,c));
                }
            }
            return a;
        }
        if (IsDiagonal(a)) {
            for (size_t i = 0; i < a.NumRows(); ++i) {
                a(i,i) = Exp(a(i,i));
            }
            return a;
        }

        try {
            std::tuple<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>> decomp = Eigendecomposition(a);
            SquareMatrix<T,N,Flags> d = std::get<1>(decomp);
            for (size_t i = 0; i < a.NumRows(); ++i) {
                d(i,i) = Exp(d(i,i));
            }
            return std::get<0>(decomp)*d*std::get<2>(decomp);
        }
        catch (...) {}

        T factorial = 1;
        Matrix<T,M,N,Flags> ret = Identity<T,Flags>(a.NumRows());
        for (size_t k = 1; k < 10; ++k) {
            factorial *= k;
            ret += Pow(a, T(k)) / factorial;
        }
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Log(Matrix<T,M,N,Flags> a) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Log(a(r,c));
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Log(Matrix<T,M,N,Flags> a, T base) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Log(a(r,c), base);
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Log(Matrix<T,M,N,Flags> a, Complex<T> base) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Log(a(r,c), base);
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Pow(Matrix<T,M,N,Flags> a, Complex<T> power) {
        if (!IsSquare(a)) {
            for (size_t r = 0; r < a.NumRows(); ++r) {
                for (size_t c = 0; c < a.NumColumns(); ++c) {
                    a(r,c) = Pow(a(r,c), power);
                }
            }
            return a;
        }
        if (power == T(0)) {
            return Identity<T,Flags>(a.NumRows());
        }

        if (IsDiagonal(a)) {
            for (size_t i = 0; i < a.NumRows(); ++i) {
                a(i,i) = Pow(a(i,i), power);
            }
            return a;
        }

        try {
            std::tuple<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>> decomp = Eigendecomposition(a);
            SquareMatrix<T,N,Flags> d = std::get<1>(decomp);
            for (size_t i = 0; i < a.NumRows(); ++i) {
                d(i,i) = Pow(d(i,i), power);
            }
            return std::get<0>(decomp)*d*std::get<2>(decomp);
        }
        catch (...) {}

        if (IsReal(power) && Floor(power.Re) == power.Re) {
            for (size_t k = 0; k < Abs(power.Re); ++k) {
                a = a*a;
            }

            if (power.Re < T(0))
                return Inverse(a);
            return a;
        }

        throw "Couldn't perform matrix power.";
    }
    template <typename T, size_t M, size_t N, unsigned int Flags, typename U>
    Matrix<T,M,N,Flags> Pow(const Matrix<T,M,N,Flags>& a, U power) {
        return Pow(a, Complex<T>(T(power), 0));
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Mod(Matrix<T,M,N,Flags> a, Complex<T> y) {
        for (size_t r = 0; r < a.NumRows(); ++r) {
            for (size_t c = 0; c < a.NumColumns(); ++c) {
                a(r,c) = Mod(a(r,c), y);
            }
        }
        return a;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags, typename U>
    Matrix<T,M,N,Flags> Mod(Matrix<T,M,N,Flags> a, U y) {
        return Mod(a, Complex<T>(T(y), 0));
    }
}
