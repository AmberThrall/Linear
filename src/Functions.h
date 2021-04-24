#pragma once
#include "Matrix.h"
#include "Vector.h"
#include "Eigen.h"
#include "Global.h"
#include "Types.h"

namespace Linear {
    /**
     * Computes the matrix \f$C\f$ defined by \f$c_{ij}=a_{ij}b_{ij}\f$. If \f$M\ne P\f$ or \f$N\ne Q\f$ an exception is thrown.
     * @param A MxN Matrix
     * @param B PxQ Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    Matrix<T,M,N,Flags> EntrywiseProduct(const Matrix<T,M,N,Flags>& A, const Matrix<T,P,Q,Flags2>& B) {
        if (A.NumRows() != B.NumRows() || A.NumColumns() != B.NumColumns())
            throw "Cannot perform entrywise product when matrices have varying sizes.";

        Matrix<T,M,N,Flags> ret(A.NumRows(),A.NumColumns(),T(0));
        for (size_t i = 0; i < A.NumRows(); ++i) {
            for (size_t j = 0; j < A.NumColumns(); ++j)
                ret(i,j) = A(i,j)*B(i,j);
        }
        return ret;
    }
    /**
     * Computes the matrix \f$C\f$ defined by \f$c_{ij}=a_{ij}/b_{ij}\f$. If \f$M\ne P\f$ or \f$N\ne Q\f$ an exception is thrown.
     * @param A MxN Matrix
     * @param B PxQ Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    Matrix<T,M,N,Flags> EntrywiseDivision(const Matrix<T,M,N,Flags>& A, const Matrix<T,P,Q,Flags2>& B) {
        if (A.NumRows() != B.NumRows() || A.NumColumns() != B.NumColumns())
            throw "Cannot perform entrywise division when matrices have varying sizes.";

        Matrix<T,M,N,Flags> ret(A.NumRows(),A.NumColumns(),T(0));
        for (size_t i = 0; i < A.NumRows(); ++i) {
            for (size_t j = 0; j < A.NumColumns(); ++j)
                ret(i,j) = A(i,j)/B(i,j);
        }
        return ret;
    }

    /**
     * If \f$A\f$ is a vector, it computes the entrywise vector p-norm \f$\left(\sum_{i=0}^{N-1}|a_i|^p\right)^{1/p}\f$.
     * Otherwise it computes the matrix norm \f$\|A\|_p\f$.
     * If p=1, it returns \f$\max_{j=0,\dots,N-1}\sum_{i=0}^{M-1}|a_{ij}|\f$.
     * If p=2, it returns the largest singular value of \f$A\f$.
     * Otherwise an error is returned.
     * @param A MxN Matrix
     * @param p Real number (default = 2)
     * @return Real number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T Norm(const Matrix<T,M,N,Flags>& A, size_t p = 2) {
        if (A.NumEntries() == 0)
            return T(0);

        if (IsVector(A)) {
            return EntrywiseNorm(A, p);
        }
        if (p == 1) {
            T ret = T(0);
            for (size_t j = 0; j < A.NumColumns(); ++j) {
                T sum = T(0);
                for (size_t i = 0; i < A.NumRows(); ++i) {
                    sum += Abs(A(i,j));
                }
                if (sum > ret)
                    ret = sum;
            }
            return ret;
        }
        else if (p == 2) {
            Vector<T,N> b = Random<T>(A.NumRows(),1);
            Eigenpair<T,N> pair = PowerIteration(ConjugateTranspose(A)*A, b, 25);
            return Sqrt(pair.value).Re;
        }
        else
            throw "Only p=1 and p=2 are supported for matrix norm.";
    }
    /**
     * Computes the entrywise p-norm \f$\left(\sum_{j=0}^{N-1}\sum_{i=0}^{M-1}|a_{ij}|^p\right)^{1/p}\f$.
     * @param A MxN Matrix
     * @param p Real number (default = 2)
     * @return Real number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T EntrywiseNorm(const Matrix<T,M,N,Flags>& A, size_t p = 2) {
        T sum = T(0);
        for (size_t j = 0; j < A.NumColumns(); ++j) {
            for (size_t i = 0; i < A.NumRows(); ++i) {
                sum += Pow(Abs(A(i,j)), T(p));
            }
        }
        return Pow(sum, 1/T(p));
    }
    /**
     * Computes the Frobenius norm. This is identical to EntrywiseNorm(A,2).
     * @param A MxN Matrix
     * @return Real number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T FrobeniusNorm(const Matrix<T,M,N,Flags>& A) {
        return EntrywiseNorm(A, 2);
    }
    /**
     * Computes the max norm \f$\max_{ij}|a_{ij}|\f$.
     * @param A MxN Matrix
     * @return Real number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T MaxNorm(const Matrix<T,M,N,Flags>& A) {
        T ret = T(0);
        for (size_t j = 0; j < A.NumColumns(); ++j) {
            for (size_t i = 0; i < A.NumRows(); ++i) {
                if (Abs(A(i,j)) > ret)
                    ret = Abs(A(i,j));
            }
        }
        return ret;
    }
    /**
     * Computes the infinity norm \f$\max_{i=0,...M-1}\sum_{j=0}^{N-1}|a_{ij}|\f$.
     * @param A MxN Matrix
     * @return Real number
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    T InfinityNorm(const Matrix<T,M,N,Flags>& A) {
        T ret = T(0);
        for (size_t i = 0; i < A.NumRows(); ++i) {
            T sum = T(0);
            for (size_t j = 0; j < A.NumColumns(); ++j) {
                sum += Abs(A(i,j));
            }
            if (sum > ret)
                ret = sum;
        }
        return ret;
    }

    /**
     * Determines if a matrix is real or complex.
     * @param A MxN Matrix
     * @return Returns true if every entry of A is real.
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsReal(const Matrix<T,M,N,Flags>& A) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                if (!IsReal(A(r,c)))
                    return false;
            }
        }
        return true;
    }

    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=|a_{ij}|\f$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Abs(Matrix<T,M,N,Flags> A) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Abs(A(r,c));
            }
        }
        return A;
    }
    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=sgn(a_{ij})\f$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Sign(Matrix<T,M,N,Flags> A) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Sign(A(r,c));
            }
        }
        return A;
    }
    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=arg(a_{ij})\f$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Arg(Matrix<T,M,N,Flags> A) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Arg(A(r,c));
            }
        }
        return A;
    }
    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=\overline{a_{ij}}\f$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Conjugate(Matrix<T,M,N,Flags> A) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Conjugate(A(r,c));
            }
        }
        return A;
    }
    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=\sqrt{a_{ij}}\f$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Sqrt(Matrix<T,M,N,Flags> A) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Sqrt(A(r,c));
            }
        }
        return A;
    }
    /**
     * If \f$A\f$ is not square, it computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=e^{a_{ij}}\f$. Otherwise it computes the matrix
     * exponential \f$e^A=\sum_{k=0}^\infty\frac{1}{k!}A^k\f$. If \f$A=VDV^{-1}\f$ with \f$D=diag(d_1,\dots,d_N)\f$ a diagonal matrix,
     * then \f$e^{A}=Ve^{D}V^{-1}\f$ and \f$e^D=diag(e^{d_1},\dots,e^{d_N})\f$. If no such decomposition can be found, it simply estimates
     * the series upto the first 10 terms.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Exp(Matrix<T,M,N,Flags> A) {
        if (!IsSquare(A)) {
            for (size_t r = 0; r < A.NumRows(); ++r) {
                for (size_t c = 0; c < A.NumColumns(); ++c) {
                    A(r,c) = Exp(A(r,c));
                }
            }
            return A;
        }
        if (IsDiagonal(A)) {
            for (size_t i = 0; i < A.NumRows(); ++i) {
                A(i,i) = Exp(A(i,i));
            }
            return A;
        }

        try {
            Eigendecomposition<T,N,Flags> eigen(A);
            for (size_t i = 0; i < eigen.D.NumRows(); ++i) {
                eigen.D(i,i) = Exp(eigen.D(i,i));
            }
            return eigen.Q*eigen.D*eigen.Qinv;
        }
        catch (...) {}

        T factorial = 1;
        Matrix<T,M,N,Flags> ret = Identity<T,Flags>(A.NumRows());
        for (size_t k = 1; k < 10; ++k) {
            factorial *= k;
            ret += Pow(A, T(k)) / factorial;
        }
        return ret;
    }
    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=\log a_{ij}\f$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Log(Matrix<T,M,N,Flags> A) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Log(A(r,c));
            }
        }
        return A;
    }
    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=\log_{base} a_{ij}\f$.
     * @param A MxN Matrix
     * @param base Real number
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Log(Matrix<T,M,N,Flags> A, T base) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Log(A(r,c), base);
            }
        }
        return A;
    }
    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=\log_{base} a_{ij}\f$.
     * @param A MxN Matrix
     * @param base Complex number
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Log(Matrix<T,M,N,Flags> A, Complex<T> base) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Log(A(r,c), base);
            }
        }
        return A;
    }
    /**
     * If \f$A\f$ is not square, it computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=a_{ij}^{power}\f$. Otherwise it attempts to computes the matrix
     * power \f$A^{power}\f$. If \f$A=VDV^{-1}\f$ with \f$D=diag(d_1,\dots,d_N)\f$ a diagonal matrix,
     * then \f$A^{power}=VD^{power}V^{-1}\f$ and \f$D^{power}=diag(d_1^{power},\dots,d_N^{power})\f$.
     * If no such decomposition can be found and power is an integer, it multiplies A by itself |power|-times then taking the inverse if power<0.
     * Otherwise it gives up and throws an exception.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Pow(Matrix<T,M,N,Flags> A, Complex<T> power) {
        if (!IsSquare(A)) {
            for (size_t r = 0; r < A.NumRows(); ++r) {
                for (size_t c = 0; c < A.NumColumns(); ++c) {
                    A(r,c) = Pow(A(r,c), power);
                }
            }
            return A;
        }
        if (power == T(0))
            return Identity<T,Flags>(A.NumRows());
        if (power == T(1))
            return A;

        if (IsDiagonal(A)) {
            for (size_t i = 0; i < A.NumRows(); ++i) {
                A(i,i) = Pow(A(i,i), power);
            }
            return A;
        }

        try {
            Eigendecomposition<T,N,Flags> eigen(A);
            for (size_t i = 0; i < eigen.D.NumRows(); ++i) {
                eigen.D(i,i) = Pow(eigen.D(i,i), power);
            }
            return eigen.Q*eigen.D*eigen.Qinv;
        }
        catch (...) {}

        if (IsReal(power) && Floor(power.Re) == power.Re) {
            SquareMatrix<T,N,Flags> B = A;
            for (size_t k = 0; k < Abs(power.Re)-1; ++k) {
                B = B*A;
            }

            if (power.Re < T(0))
                return Inverse(B);
            return B;
        }

        throw "Couldn't perform matrix power.";
    }
    template <typename T, size_t M, size_t N, unsigned int Flags, typename U>
    Matrix<T,M,N,Flags> Pow(const Matrix<T,M,N,Flags>& A, U power) {
        return Pow(A, Complex<T>(T(power), 0));
    }

    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=a_{ij}\bmod z\f$.
     * @param A MxN Matrix
     * @param z Complex number
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Mod(Matrix<T,M,N,Flags> A, Complex<T> z) {
        for (size_t r = 0; r < A.NumRows(); ++r) {
            for (size_t c = 0; c < A.NumColumns(); ++c) {
                A(r,c) = Mod(A(r,c), z);
            }
        }
        return A;
    }
    /**
     * Computes the MxN Matrix \f$B\f$ defined by \f$b_{ij}=a_{ij}\bmod y\f$.
     * @param A MxN Matrix
     * @param y Real number
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags, typename U>
    Matrix<T,M,N,Flags> Mod(Matrix<T,M,N,Flags> A, U y) {
        return Mod(A, Complex<T>(T(y), 0));
    }

    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\sin a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Sin(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Sin(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\cos a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Cos(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Cos(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\tan a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Tan(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Tan(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\csc a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Csc(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Csc(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\sec a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Sec(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Sec(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\cot a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Cot(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Cot(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\sin^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ASin(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ASin(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\cos^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ACos(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ACos(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\tan^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ATan(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ATan(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\csc^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ACsc(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ACsc(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\sec^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ASec(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ASec(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\cot^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ACot(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ACot(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\sinh a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Sinh(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Sinh(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\cosh a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Cosh(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Cosh(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\tanh a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Tanh(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Tanh(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=csch a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Csch(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Csch(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=sech a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Sech(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Sech(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\coth a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Coth(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = Coth(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\sinh^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ASinh(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ASinh(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\cosh^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ACosh(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ACosh(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\tanh^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ATanh(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ATanh(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=csch^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ACsch(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ACsch(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=sech^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ASech(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ASech(A(r,c));
              }
         }
         return A;
    }
    /**
     * Computes the MxN Matrix$B$ defined by $b_{ij}=\coth^{-1} a_{ij}$.
     * @param A MxN Matrix
     * @return MxN Matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> ACoth(Matrix<T,M,N,Flags> A) {
         for (size_t r = 0; r < A.NumRows(); ++r) {
              for (size_t c = 0; c < A.NumColumns(); ++c) {
                   A(r,c) = ACoth(A(r,c));
              }
         }
         return A;
    }
}
