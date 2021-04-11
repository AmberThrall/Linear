#pragma once
#include <sstream>
#include "Complex.h"

namespace Linear {
    const unsigned int Dynamic = 0;
    const unsigned int RowMajor = 0x0000;
    const unsigned int ColumnMajor = 0x0001;

    typedef unsigned int size_t;
    template<typename T, size_t M, size_t N, unsigned int Flags = 0>
    class Matrix {
    public:
        Matrix() : Matrix(0) {}
        Matrix(T val) {
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
            this->m = nrows;
            this->n = ncols;
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < this->m*this->n; ++i) {
                this->data[i] = Complex<T>(val, 0);
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
        template <typename=std::enable_if<(M!=Dynamic&&N!=Dynamic)>>
        Matrix(std::initializer_list<Complex<T>> list) {
            this->m = M;
            this->n = N;
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
        template <typename U, size_t P, size_t Q, unsigned int Flags2>
        Matrix(const Matrix<U,P,Q,Flags2>& copy) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = copy.NumRows(); }
            if (N == Dynamic) { this->n = copy.NumColumns(); }
            this->data = new Complex<T>[this->m*this->n];

            size_t r = 0;
            for (; r < NumRows() && r < copy.NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns() && c < copy.NumColumns(); ++c) {
                    (*this)(r,c) = (Complex<T>)copy(r,c);
                }
            }

            // Fill the remaining entries with zero.
            for (; r < NumRows(); ++r) {
                for (size_t c = 0; c < NumColumns(); ++c) {
                    (*this)(r,c) = T(0);
                }
            }
        }

        size_t NumRows() const { return this->m; }
        size_t NumColumns() const { return this->n; }
        size_t Size() const { return NumRows()*NumColumns(); }

        /// Operators.
        // Access operators
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
    private:
        Complex<T> * data;
        size_t m, n;
    };

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

    template <typename T,unsigned int N>
    using Vector = Matrix<T,N,1>;
    using Vector2i = Matrix<int,2,1>;
    using Vector3i = Matrix<int,3,1>;
    using Vector4i = Matrix<int,4,1>;
    using VectorXi = Matrix<int,Dynamic,1>;
    using Vector2f = Matrix<float,2,1>;
    using Vector3f = Matrix<float,3,1>;
    using Vector4f = Matrix<float,4,1>;
    using VectorXf = Matrix<float,Dynamic,1>;
    using Vector2d = Matrix<double,2,1>;
    using Vector3d = Matrix<double,3,1>;
    using Vector4d = Matrix<double,4,1>;
    using VectorXd = Matrix<double,Dynamic,1>;

    template <typename T,unsigned int N>
    using RowVector = Matrix<T,1,N>;
    using RowVector2i = Matrix<int,1,2>;
    using RowVector3i = Matrix<int,1,3>;
    using RowVector4i = Matrix<int,1,4>;
    using RowVectorXi = Matrix<int,1,Dynamic>;
    using RowVector2f = Matrix<float,1,2>;
    using RowVector3f = Matrix<float,1,3>;
    using RowVector4f = Matrix<float,1,4>;
    using RowVectorXf = Matrix<float,1,Dynamic>;
    using RowVector2d = Matrix<double,1,2>;
    using RowVector3d = Matrix<double,1,3>;
    using RowVector4d = Matrix<double,1,4>;
    using RowVectorXd = Matrix<double,1,Dynamic>;
}
