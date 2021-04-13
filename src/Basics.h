#pragma once
#include "Matrix.h"
#include "Vector.h"
#include "Global.h"

namespace Linear {
    template <size_t P, size_t Q, typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(P>0&&Q>0),Matrix<T,P,Q,Flags>>::type SubMatrix(Matrix<T,M,N,Flags> m, size_t i = 0, size_t j = 0) {
        if (i+P > m.NumRows() || j+Q > m.NumColumns())
            throw "Cannot create submatrix, indices out of bounds.";
        Matrix<T,P,Q,Flags> ret(P,Q,T(0));
        for (size_t r = 0; r < P; ++r) {
            for (size_t c = 0; c < Q; ++c) {
                ret(r,c) = m(i+r,j+c);
            }
        }
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,Dynamic,Dynamic,Flags> SubMatrix(Matrix<T,M,N,Flags> m, size_t nrows, size_t ncols, size_t i = 0, size_t j = 0) {
        if (i+nrows > m.NumRows() || j+ncols > m.NumColumns())
            throw "Cannot create submatrix, indices out of bounds.";
        Matrix<T,Dynamic,Dynamic,Flags> ret(nrows,ncols,T(0));
        for (size_t r = 0; r < nrows; ++r) {
            for (size_t c = 0; c < ncols; ++c) {
                ret(r,c) = m(i+r,j+c);
            }
        }
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,(M==Dynamic?Dynamic:M-1),N,Flags> RemoveRow(Matrix<T,M,N,Flags> m, size_t i) {
        if (m.NumRows() == 1)
            throw "Cannot create a 0xN matrix.";
        Matrix<T,(M==Dynamic?Dynamic:M-1),N,Flags> ret(m.NumRows()-1,m.NumColumns(),T(0));
        for (size_t r = 0; r < m.NumRows(); ++r) {
            for (size_t c = 0; c < m.NumColumns(); ++c) {
                if (r < i)
                    ret(r,c) = m(r,c);
                if (r > i)
                    ret(r-1,c) = m(r,c);
            }
        }
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,(N==Dynamic?Dynamic:N-1),Flags> RemoveColumn(Matrix<T,M,N,Flags> m, size_t i) {
        if (m.NumColumns() == 1)
            throw "Cannot create a Mx0 matrix.";
        Matrix<T,M,(N==Dynamic?Dynamic:N-1),Flags> ret(m.NumRows(),m.NumColumns()-1,T(0));
        for (size_t r = 0; r < m.NumRows(); ++r) {
            for (size_t c = 0; c < m.NumColumns(); ++c) {
                if (c < i)
                    ret(r,c) = m(r,c);
                if (c > i)
                    ret(r,c-1) = m(r,c);
            }
        }
        return ret;
    }
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,(M==Dynamic?Dynamic:M-1),(N==Dynamic?Dynamic:N-1),Flags> RemoveRowAndColumn(Matrix<T,M,N,Flags> m, size_t i, size_t j) {
        if (m.NumRows() == 1 || m.NumColumns() == 1)
            throw "Cannot create a 0x0 matrix.";
        Matrix<T,(M==Dynamic?Dynamic:M-1),(N==Dynamic?Dynamic:N-1),Flags> ret(m.NumRows()-1,m.NumColumns()-1,T(0));
        for (size_t r = 0; r < m.NumRows(); ++r) {
            for (size_t c = 0; c < m.NumColumns(); ++c) {
                if (r < i && c < j)
                    ret(r,c) = m(r,c);
                if (r > i && c < j)
                    ret(r-1,c) = m(r,c);
                if (r < i && c > j)
                    ret(r,c-1) = m(r,c);
                if (r > i && c > j)
                    ret(r-1,c-1) = m(r,c);
            }
        }
        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,N,M,Flags> Transpose(const Matrix<T,M,N,Flags>& m) {
        Matrix<T,N,M,Flags> ret(m.NumColumns(), m.NumRows(), T(0));
        for (size_t r = 0; r < m.NumColumns(); ++r) {
            for (size_t c = 0; c < m.NumRows(); ++c) {
                ret(r,c) = m(c,r);
            }
        }
        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,N,M,Flags> ConjugateTranspose(Matrix<T,M,N,Flags> m) {
        Matrix<T,N,M,Flags> ret(m.NumColumns(), m.NumRows(), T(0));
        for (size_t r = 0; r < m.NumColumns(); ++r) {
            for (size_t c = 0; c < m.NumRows(); ++c) {
                ret(r,c) = Conjugate(m(c,r));
            }
        }
        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> RREF(Matrix<T,M,N,Flags> m) {
        size_t lead = 0;
        for (size_t r = 0; r < m.NumRows(); ++r) {
            if (lead >= m.NumColumns())
                return m;

            // Find the pivot.
            size_t i = r;
            while (m(i, lead) == 0) {
                i += 1;
                if (i == m.NumRows()) {
                    i = r;
                    lead += 1;
                    if (lead == m.NumColumns())
                        return m;
                }
            }

            // Swap rows i and r
            m.SwapRows(i, r);
            // R_r / pivot -> R_r
            m.ScaleRow(r, 1/m(r,lead));
            // For each row i!=r, R_i - R_r*m(i,lead) -> R_i
            for (size_t i = 0; i < m.NumRows(); ++i) {
                if (i == r) continue;

                m.AddRows(i, r, -m(i,lead));
            }

            lead += 1;
        }
        return m;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    bool IsSquare(Matrix<T,M,N,Flags> m) {
        return (m.NumRows() == m.NumColumns());
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Complex<T>>::type Trace(Matrix<T,M,N,Flags> m) {
        if (!IsSquare(m))
            throw "Cannot take the trace of a non-square matrix.";
        Complex<T> ret;
        for (size_t r = 0; r < m.NumRows(); ++r) {
            ret += m(r,r);
        }
        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Complex<T>>::type Determinant(Matrix<T,M,N,Flags> m) {
        if (!IsSquare(m))
            throw "Cannot take the determinant of a non-square matrix.";

        if (m.NumRows() == 1)
            return m(0,0);
        if (m.NumRows() == 2)
            return (m(0,0)*m(1,1) - m(0,1)*m(1,0));

        Complex<T> ret;
        for (size_t i = 0; i < m.NumRows(); ++i) {
            ret += m(0,i)*Cofactor(m, 0, i);
        }
        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Complex<T>>::type Minor(Matrix<T,M,N,Flags> m, size_t i, size_t j) {
        if (!IsSquare(m))
            throw "Cannot take the minor of a non-square matrix.";
        if (i >= m.NumRows() || j >= m.NumRows())
            throw "Minor coordinates out of bounds.";
        if (m.NumRows() == 1)
            return m(0,0);

        Complex<T> ret;
        return Determinant(RemoveRowAndColumn(m, i, j));
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Complex<T>>::type Cofactor(Matrix<T,M,N,Flags> m, size_t i, size_t j) {
        if (!IsSquare(m))
            throw "Cannot take the cofactor of a non-square matrix.";
        if (i >= m.NumRows() || j >= m.NumRows())
            throw "Cofactor coordinates out of bounds.";
        return T(std::pow(-1.0,i+j))*Minor(m, i,j);
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Matrix<T,M,N,Flags>>::type Adjugate(Matrix<T,M,N,Flags> m) {
        if (!IsSquare(m))
            throw "Cannot take the adjugate of a non-square matrix.";
        Matrix<T,M,N,Flags> ret(m.NumRows(), m.NumColumns(), T(0));
        for (size_t r = 0; r < m.NumRows(); ++r) {
            for (size_t c = 0; c < m.NumColumns(); ++c) {
                ret(c,r) = Cofactor(m, r,c);
            }
        }
        return ret;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N)||M==Dynamic||N==Dynamic), Matrix<T,M,N,Flags>>::type Inverse(Matrix<T,M,N,Flags> m) {
        if (!IsSquare(m))
            throw "Cannot take the inverse of a non-square matrix.";
        Complex<T> det = Determinant(m);
        if (det == 0)
            throw "Cannot take the inverse of a singular matrix.";
        return Adjugate(m)/det;
    }
}
