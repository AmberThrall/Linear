#pragma once
#include <sstream>
#include "Complex.h"
#include "Global.h"

namespace Linear {
    const unsigned int Dynamic = 0;
    const unsigned int RowMajor = 0x0000;
    const unsigned int ColumnMajor = 0x0001;

    template<typename T, size_t M, size_t N, unsigned int Flags = 0>
    class Matrix {
    public:
        Matrix(T val = T(0)) {
            if (M == Dynamic || N == Dynamic)
                throw "Cannot create a 0x0 matrix.";
            this->m = M;
            this->n = N;
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < Size(); ++i) {
                this->data[i] = val;
            }
        }
        Matrix(Complex<T> val) {
            if (M == Dynamic || N == Dynamic)
                throw "Cannot create a 0x0 matrix.";
            this->m = M;
            this->n = N;
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < Size(); ++i) {
                this->data[i] = val;
            }
        }
        Matrix(size_t size, T val) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = size; }
            if (N == Dynamic) { this->n = size; }
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < Size(); ++i) {
                this->data[i] = val;
            }
        }
        Matrix(size_t size, Complex<T> val) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = size; }
            if (N == Dynamic) { this->n = size; }
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < Size(); ++i) {
                this->data[i] = val;
            }
        }
        Matrix(size_t nrows, size_t ncols, T val) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = nrows; }
            if (N == Dynamic) { this->n = ncols; }
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < this->m*this->n; ++i) {
                this->data[i] = val;
            }
        }
        Matrix(size_t nrows, size_t ncols, Complex<T> val) {
            this->m = nrows;
            this->n = ncols;
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < this->m*this->n; ++i) {
                this->data[i] = val;
            }
        }
        Matrix(std::initializer_list<T> list) {
            this->m = M;
            this->n = N;
            if (M == Dynamic && N == Dynamic) { this->m = list.size(); this->n = 1; }
            if (M != Dynamic && N == Dynamic) { this->n = (list.size() >= this->m ? list.size() / this->m : 1); }
            if (M == Dynamic && N != Dynamic) { this->m = (list.size() >= this->n ? list.size() / this->n : 1); }
            this->data = new Complex<T>[this->m*this->n];
            size_t i = 0;
            for (auto it = std::begin(list); it != std::end(list); ++it) {
                this->data[i] = *it;
                i += 1;
                if (i == Size())
                    break;
            }

            // Fill the remaining entries with zero.
            for (; i < Size(); ++i)
                this->data[i] = T(0);
        }
        Matrix(std::initializer_list<Complex<T>> list) {
            this->m = M;
            this->n = N;
            if (M == Dynamic && N == Dynamic) { this->m = list.size(); this->n = 1; }
            if (M != Dynamic && N == Dynamic) { this->n = (list.size() >= this->m ? list.size() / this->m : 1); }
            if (M == Dynamic && N != Dynamic) { this->m = (list.size() >= this->n ? list.size() / this->n : 1); }
            this->data = new Complex<T>[this->m*this->n];
            size_t i = 0;
            for (auto it = std::begin(list); it != std::end(list); ++it) {
                this->data[i] = *it;
                i += 1;
                if (i == Size())
                    break;
            }

            // Fill the remaining entries with zero.
            for (; i < Size(); ++i)
                this->data[i] = T(0);
        }
        Matrix(std::initializer_list<std::initializer_list<T>> list) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = list.size(); }
            if (N == Dynamic) {
                for (auto it = std::begin(list); it != std::end(list); ++it) {
                    if ((*it).size() < this->n || this->n == Dynamic)
                        this->n = (*it).size();
                }
            }
            this->data = new Complex<T>[this->m*this->n];
            size_t r = 0;
            for (auto it = std::begin(list); it != std::end(list); ++it) {
                size_t c = 0;
                for (auto it2 = std::begin(*it); it2 != std::end(*it); ++it2) {
                    (*this)(r,c) = *it2;
                    c += 1;
                }
                // Fill the remaining entries with zero.
                for (; c < NumColumns(); ++c)
                    (*this)(r,c) = T(0);

                r += 1;
                if (r == NumRows())
                    break;
            }

            // Fill the remaining entries with zero.
            for (; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) = T(0);
                }
            }
        }
        Matrix(std::initializer_list<std::initializer_list<Complex<T>>> list) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = list.size(); }
            if (N == Dynamic) {
                for (auto it = std::begin(list); it != std::end(list); ++it) {
                    if ((*it).size() < this->n || this->n == Dynamic)
                        this->n = (*it).size();
                }
            }
            this->data = new Complex<T>[this->m*this->n];
            size_t r = 0;
            for (auto it = std::begin(list); it != std::end(list); ++it) {
                size_t c = 0;
                for (auto it2 = std::begin(*it); it2 != std::end(*it); ++it2) {
                    (*this)(r,c) = *it2;
                    c += 1;
                }
                // Fill the remaining entries with zero.
                for (; c < NumColumns(); ++c)
                    (*this)(r,c) = T(0);

                r += 1;
                if (r == NumRows())
                    break;
            }

            // Fill the remaining entries with zero.
            for (; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) = T(0);
                }
            }
        }
        Matrix(const Matrix<T,M,N,Flags>& copy) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = copy.NumRows(); }
            if (N == Dynamic) { this->n = copy.NumColumns(); }
            this->data = new Complex<T>[this->m*this->n];

            for (size_t r = 0; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) = copy(r,c);
                }
            }
        }
        template<typename U,size_t P, size_t Q, unsigned int Flags2>
        Matrix(const Matrix<U,P,Q,Flags2>& other) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = other.NumRows(); }
            if (N == Dynamic) { this->n = other.NumColumns(); }
            this->data = new Complex<T>[this->m*this->n];
            if (NumRows() != other.NumRows() || NumColumns() != other.NumColumns())
                throw "Cannot assign matrix to a matrix of different size.";

            for (size_t r = 0; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) = (Complex<T>)other(r,c);
                }
            }
        }

        size_t NumRows() const { return this->m; }
        size_t NumColumns() const { return this->n; }
        size_t Size() const { return NumRows()*NumColumns(); }

        Matrix<T,1,N> GetRow(size_t r) const {
            Matrix<T,1,N> ret(NumColumns(), T(0));
            for (size_t c = 0; c < NumColumns(); ++c)
                ret(0,c) = (*this)(r,c);
            return ret;
        }
        Matrix<T,M,1> GetColumn(size_t c) const {
            Matrix<T,M,1> ret(NumRows(),T(0));
            for (size_t r = 0; r < NumRows(); ++r)
                ret(r,0) = (*this)(r,c);
            return ret;
        }
        template <size_t P, size_t Q>
        typename std::enable_if<(P==1||P==Dynamic)&&(Q==N||Q==Dynamic||N==Dynamic),void>::type SetRow(size_t r, const Matrix<T,P,Q>& row) {
            if (row.NumRows() != 1 || row.NumColumns() != NumColumns())
                throw "Expected a row vector in Matrix::SetRow()";
            for (size_t c = 0; c < NumColumns(); ++c)
                (*this)(r,c) = row(0,c);
        }
        template <size_t P, size_t Q>
        typename std::enable_if<(P==M||P==Dynamic||M==Dynamic)&&(Q==1||Q==Dynamic),void>::type SetColumn(size_t c, const Matrix<T,P,Q>& column) {
            if (column.NumRows() != NumRows() || column.NumColumns() != 1)
                throw "Expected a column vector in Matrix::SetColumn()";
            for (size_t r = 0; r < NumRows(); ++r)
                (*this)(r,c) = column(r,0);
        }
        // R_{r1} <-> R_{r2}
        void SwapRows(size_t r1, size_t r2) {
            for (size_t c = 0; c < NumColumns(); ++c) {
                Complex<T> temp = (*this)(r1,c);
                (*this)(r1,c) = (*this)(r2,c);
                (*this)(r2,c) = temp;
            }
        }
        // s*R_r -> R_r
        void ScaleRow(size_t r, Complex<T> s) {
            for (size_t c = 0; c < NumColumns(); ++c)
                (*this)(r,c) *= s;
        }
        // R_{r1}+s*R_{r2} -> R_{r1}
        void AddRows(size_t r1, size_t r2, Complex<T> s) {
            for (size_t c = 0; c < NumColumns(); ++c)
                (*this)(r1,c) += s*(*this)(r2,c);
        }

        /// Operators.
        // Access operators
        Complex<T> operator[] (size_t i) const { return this->data[i]; }
        Complex<T> & operator[] (size_t i) { return this->data[i]; }
        Complex<T> operator() (size_t r, size_t c) const {
            if (Flags & ColumnMajor)
                return this->data[c*this->m+r];
            return this->data[r*this->n+c];
        }
        Complex<T> & operator() (size_t r, size_t c) {
            if (Flags & ColumnMajor)
                return this->data[c*this->m+r];
            return this->data[r*this->n+c];
        }
        // Type conversion.
        template<typename U, typename std::enable_if<std::is_convertible<T,U>::value>::type* = nullptr>
        operator Matrix<U,M,N,Flags>() {
            return Matrix<U,M,N,Flags>(*this);
        }
        // Assignments
        template <size_t P, size_t Q, unsigned int Flags2>
        Matrix<T,M,N,Flags> & operator=(const Matrix<T,P,Q,Flags2>& other) {
            if (M == Dynamic) { this->m = other.NumRows(); }
            if (N == Dynamic) { this->n = other.NumColumns(); }
            if (NumRows() != other.NumRows() || NumColumns() != other.NumColumns())
                throw "Cannot assign matrix to a matrix of different size.";
            for (size_t r = 0; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) = other(r,c);
                }
            }
            return *this;
        }
        template <size_t P, size_t Q, unsigned int Flags2>
        Matrix<T,M,N,Flags> & operator+=(const Matrix<T,P,Q,Flags2>& other) {
            if (NumRows() != other.NumRows() || NumColumns() != other.NumColumns())
                throw "Cannot add two matrices of differing sizes.";
            for (size_t r = 0; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) += other(r,c);
                }
            }
            return *this;
        }
        template <size_t P, size_t Q, unsigned int Flags2>
        Matrix<T,M,N,Flags> & operator-=(const Matrix<T,P,Q,Flags2>& other) {
            if (NumRows() != other.NumRows() || NumColumns() != other.NumColumns())
                throw "Cannot subtract two matrices of differing sizes.";
            for (size_t r = 0; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) -= other(r,c);
                }
            }
            return *this;
        }
        template <size_t P, size_t Q, unsigned int Flags2>
        Matrix<T,M,N,Flags> & operator*=(const Matrix<T,P,Q,Flags2>& other) {
            if (NumColumns() != other.NumRows() || NumColumns() != other.NumColumns())
                throw "Cannot muliply two matrices due to size mismatch.";
            *this = (*this) * other;
            return *this;
        }
        Matrix<T,M,N,Flags> & operator*=(Complex<T> other) {
            for (size_t r = 0; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) *= other;
                }
            }
            return *this;
        }
        Matrix<T,M,N,Flags> & operator/=(Complex<T> other) {
            for (size_t r = 0; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) /= other;
                }
            }
            return *this;
        }
        // Binary operators.
        template <size_t P, size_t Q, unsigned int Flags2>
        friend Matrix<T,M,N,Flags> operator+(Matrix<T,M,N,Flags> a, const Matrix<T,P,Q,Flags2>& b) { return a += b; }
        template <size_t P, size_t Q, unsigned int Flags2>
        friend Matrix<T,M,N,Flags> operator-(Matrix<T,M,N,Flags> a, const Matrix<T,P,Q,Flags2>& b) { return a -= b; }
        friend Matrix<T,M,N,Flags> operator*(Matrix<T,M,N,Flags> a, const Complex<T>& b) { return a *= b; }
        friend Matrix<T,M,N,Flags> operator*(const Complex<T>& b, Matrix<T,M,N,Flags> a) { return a *= b; }
        friend Matrix<T,M,N,Flags> operator/(Matrix<T,M,N,Flags> a, const Complex<T>& b) { return a /= b; }
        template <size_t P, size_t Q, unsigned int Flags2>
        friend Matrix<T,M,Q,Flags> operator*(Matrix<T,M,N,Flags> a, const Matrix<T,P,Q,Flags2>& b) {
            if (a.NumColumns() != b.NumRows())
                throw "Cannot muliply two matrices due to size mismatch.";
            Matrix<T,M,Q,Flags> ret(a.NumRows(), b.NumColumns(), T(0));
            for (size_t r = 0; r < a.NumRows(); ++r) {
                for (size_t c = 0; c < b.NumColumns(); ++c) {
                    for (size_t i = 0; i < a.NumColumns(); ++i) {
                        ret(r,c) += a(r,i)*b(i,c);
                    }
                }
            }
            return ret;
        }
        // Unary operators.
        Matrix<T,M,N,Flags> operator-() {
            Matrix<T,M,N,Flags> ret(NumRows(),NumColumns(),T(0));
            for (size_t r = 0; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    ret(r,c) = -(*this)(r,c);
                }
            }
            return ret;
        }
        // Output operator
        friend std::ostream& operator<<(std::ostream& out, const Matrix<T,M,N,Flags>& m) {
            out << m.NumRows() << "x" << m.NumColumns() << std::endl;
            unsigned int longest = 5;
            for (size_t r = 0; r < m.NumRows(); ++r) {
                for (size_t c = 0; c < m.NumColumns(); ++c) {
                    std::stringstream stream;
                    stream << m(r,c);
                    std::string asstring = stream.str();
                    if (asstring.length() > longest) {
                        longest = asstring.length();
                    }
                }
            }
            longest += 1;

            std::string padding = "";
            for (unsigned int i = 0; i <= longest; ++i)
                padding += " ";

            for (size_t r = 0; r < m.NumRows(); ++r) {
                for (size_t c = 0; c < m.NumColumns(); ++c) {
                    std::stringstream stream;
                    stream << m(r,c);
                    std::string asstring = stream.str();
                    out << padding.substr(0,longest-asstring.length()) << asstring.c_str();
                }
                out << std::endl;
            }
            return out;
        }
        // Comparison
        template <size_t P, size_t Q, unsigned int Flags2>
        friend bool operator==(const Matrix<T,M,N,Flags>& a, const Matrix<T,P,Q,Flags2>& b) {
            if (a.NumRows() != b.NumRows() || a.NumColumns() != b.NumColumns())
                return false;
            for (size_t r = 0; r < a.NumRows(); ++r) {
                for (size_t c = 0; c < a.NumColumns(); ++c) {
                    if (a(r,c) != b(r,c))
                        return false;
                }
            }
            return true;
        }
        template <size_t P, size_t Q, unsigned int Flags2>
        friend bool operator!=(const Matrix<T,M,N,Flags>& a, const Matrix<T,P,Q,Flags2>& b) { return !(a==b); }
    private:
        Complex<T> * data;
        size_t m, n;
    };

    template <typename T,size_t N, unsigned int Flags = 0>
    using SquareMatrix = Matrix<T,N,N,Flags>;

    using Matrix2i = Matrix<int,2,2>;
    using Matrix3i = Matrix<int,3,3>;
    using Matrix4i = Matrix<int,4,4>;
    using MatrixXi = Matrix<int,Dynamic,Dynamic>;
    using Matrix2f = Matrix<float,2,2>;
    using Matrix3f = Matrix<float,3,3>;
    using Matrix4f = Matrix<float,4,4>;
    using MatrixXf = Matrix<float,Dynamic,Dynamic>;
    using Matrix2d = Matrix<double,2,2>;
    using Matrix3d = Matrix<double,3,3>;
    using Matrix4d = Matrix<double,4,4>;
    using MatrixXd = Matrix<double,Dynamic,Dynamic>;
}
