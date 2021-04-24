#pragma once
#include <sstream>
#include "Complex.h"
#include "Global.h"

namespace Linear {
    const unsigned int Dynamic = 0;
    const unsigned int RowMajor = 0x0000;
    const unsigned int ColumnMajor = 0x0001;

    /**
     * Class for matrices.
     * @param T Type to store the real and imaginary part of each entry as.
     * @param M Number of rows. Use Dynamic to allow this value to change over time.
     * @param N Number of column. Use Dynamic to allow this value to change over time.
     * @param Flags Whether to use row major storage or column major storage. (default = row major).
     */
    template<typename T, size_t M, size_t N, unsigned int Flags = 0>
    class Matrix {
    public:
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ where \f$a_{ij}=x\f$ for all \f$i,j\f$. If M or N are set to Dynamic, a 0x0 matrix will be created.
         * @param x Real number (default = 0)
         */
        Matrix(T x = T(0)) {
            if (M == Dynamic || N == Dynamic) {
                this->m = 0;
                this->n = 0;
                this->data = NULL;
                return;
            }
            this->m = M;
            this->n = N;
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < Size(); ++i) {
                this->data[i] = x;
            }
        }
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ where \f$a_{ij}=z\f$ for all \f$i,j\f$. If M or N are set to Dynamic, a 0x0 matrix will be created.
         * @param val Complex number
         */
        Matrix(Complex<T> z) {
            if (M == Dynamic || N == Dynamic) {
                this->m = 0;
                this->n = 0;
                this->data = NULL;
                return;
            }
            this->m = M;
            this->n = N;
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < Size(); ++i) {
                this->data[i] = z;
            }
        }
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ where \f$a_{ij}=x\f$ for all \f$i,j\f$. If M or N are set to Dynamic, they will be set to size. Otherwise,
         * the size parameter is ignored.
         * @param size Size to set M or N if they are set to Dynamic.
         * @param x Real number (default = 0)
         */
        Matrix(size_t size, T x) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = size; }
            if (N == Dynamic) { this->n = size; }
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < Size(); ++i) {
                this->data[i] = x;
            }
        }
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ where \f$a_{ij}=z\f$ for all \f$i,j\f$. If M or N are set to Dynamic, they will be set to size. Otherwise,
         * the size parameter is ignored.
         * @param size Size to set M or N if they are set to Dynamic.
         * @param z Complex number
         */
        Matrix(size_t size, Complex<T> z) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = size; }
            if (N == Dynamic) { this->n = size; }
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < Size(); ++i) {
                this->data[i] = z;
            }
        }
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ where \f$a_{ij}=x\f$ for all \f$i,j\f$. If M is set to Dynamic, it will be set to nrows.
         * If N is set to Dynamic, it will be set to ncols. Otherwise, the nrows and ncols parameters are ignored.
         * @param nrows Size to set M if it is set to Dynamic.
         * @param ncols Size to set N if it is set to Dynamic.
         * @param x Real number (default = 0)
         */
        Matrix(size_t nrows, size_t ncols, T x) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = nrows; }
            if (N == Dynamic) { this->n = ncols; }
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < this->m*this->n; ++i) {
                this->data[i] = x;
            }
        }
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ where \f$a_{ij}=z\f$ for all \f$i,j\f$. If M is set to Dynamic, it will be set to nrows.
         * If N is set to Dynamic, it will be set to ncols. Otherwise, the nrows and ncols parameters are ignored.
         * @param nrows Size to set M if it is set to Dynamic.
         * @param ncols Size to set N if it is set to Dynamic.
         * @param z Complex number
         */
        Matrix(size_t nrows, size_t ncols, Complex<T> z) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = nrows; }
            if (N == Dynamic) { this->n = ncols; }
            this->data = new Complex<T>[this->m*this->n];
            for (size_t i = 0; i < this->m*this->n; ++i) {
                this->data[i] = z;
            }
        }
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ by filling in entries from list. The order they are filled in depend on Flags. If the list
         * is longer than M*N, then the remaining entries are ignored. If the list is shorter than M*N, then the remaining entries are set
         * to zero.
         *
         * If both M and N is set to Dynamic, it will set M to the list's size and N to 1.
         * If just M is Dynamic or just N is Dynamic, it will set the Dynamic size to the list's size divided the non-Dynamic size.
         * However, if the list's size is less than the non-Dynamic size, it will set the Dynamic size to 1.
         * @param list Initializer list of real numbers
         */
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
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ by filling in entries from list. The order they are filled in depend on Flags. If the list
         * is longer than M*N, then the remaining entries are ignored. If the list is shorter than M*N, then the remaining entries are set
         * to zero.
         *
         * If both M and N is set to Dynamic, it will set M to the list's size and N to 1.
         * If just M is Dynamic or just N is Dynamic, it will set the Dynamic size to the list's size divided the non-Dynamic size.
         * However, if the list's size is less than the non-Dynamic size, it will set the Dynamic size to 1.
         * @param list Initializer list of Complex numbers
         */
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
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ by filling in entries from list. Each initializer list in list represents a row of A.
         * If that row is longer than N, than the remaining elements of that row are ignored. If that row is shorter than N, the remaining
         * elements are zero to zero.
         *
         * If M is set to Dynamic, it will set M to the number of rows in the list.
         * If N is set to Dynamic, it will set N to the longest row on the list.
         * @param list Initializer list of initializer lists of real numbers.
         */
        Matrix(std::initializer_list<std::initializer_list<T>> list) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = list.size(); }
            if (N == Dynamic) {
                for (auto it = std::begin(list); it != std::end(list); ++it) {
                    if ((*it).size() > this->n)
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
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ by filling in entries from list. Each initializer list in list represents a row of A.
         * If that row is longer than N, than the remaining elements of that row are ignored. If that row is shorter than N, the remaining
         * elements are zero to zero.
         *
         * If M is set to Dynamic, it will set M to the number of rows in the list.
         * If N is set to Dynamic, it will set N to the longest row on the list.
         * @param list Initializer list of initializer lists of complex numbers.
         */
        Matrix(std::initializer_list<std::initializer_list<Complex<T>>> list) {
            this->m = M;
            this->n = N;
            if (M == Dynamic) { this->m = list.size(); }
            if (N == Dynamic) {
                for (auto it = std::begin(list); it != std::end(list); ++it) {
                    if ((*it).size() > this->n)
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
        /**
         * Constructor.
         * Copies the MxN matrix other.
         *
         * @param other MxN Matrix
         */
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
        /**
         * Constructor.
         * Creates the M-by-N matrix \f$A\f$ by copying the matrix given by other. If Q < N, the remaining N-Q elements of each row will be
         * set to zero. If Q > N, the remaining Q-N elements of other is ignored. If P < M, the remaining M-P rows will be set to zero.
         * If P > M, the remaining P-M rows of other is ignored.
         *
         * If M is set to Dynamic, it will set M to P.
         * If N is set to Dynamic, it will set N to Q.
         * @param other PxQ Matrix
         */
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

        /**
         * Resizes a dynamic matrix.
         * If M or N is dynamic, this function becomes available. It resizes the matrix to newSizexN or MxnewSize depending on if
         * N or M is dynamic. If both M and N are dynamic, it becomes an newSizexnewSize matrix. If neither M nor N are dynamic, this
         * function throws an exception.
         * @param newSize New number of rows/columns
         */
        void Resize(size_t newSize) {
            if (M != Dynamic && N != Dynamic)
                throw "Cannot resize statically sized matrices.";

            size_t oldM = this->m, oldN = this->n;
            if (M == Dynamic) { this->m = newSize; }
            else { this->m = M; }
            if (N == Dynamic) { this->n = newSize; }
            else { this->n = N; }

            if (this->m == 0 || this->n == 0) {
                this->m = 0;
                this->n = 0;
                if (this->data != NULL)
                    delete this->data;
                this->data = NULL;
                return;
            }

            Complex<T> * oldData = this->data;
            this->data = new Complex<T>[this->m*this->n];
            for (size_t r = 0; r < this->m; ++r) {
                for (size_t c = 0; c < this->n; ++c) {
                    if (r >= oldM || c >= oldN) {
                        if (Flags & ColumnMajor)
                            this->data[c*this->m+r] = T(0);
                        else
                            this->data[r*this->n+c] = T(0);
                    }
                    else if (oldData != NULL) {
                        if (Flags & ColumnMajor)
                            this->data[c*this->m+r] = oldData[c*oldM+r];
                        else
                            this->data[r*this->n+c] = oldData[r*oldN+c];
                    }
                }
            }
            if (oldData != NULL)
                delete oldData;
        }
        /**
         * Resizes a dynamic matrix.
         * If M or N are dynamic, this function becomes available. If both M and N are dynamic, it resizes the matrix to newMxnewN.
         * Otherwise it resizes the matrix to newMxN or MxnewN depending if M or N is dynamic. If neither M nor N are dynamic, this
         * function throws an exception.
         * @param newM New number of rows. Ignored if M is not dynamic
         * @param newN New number of columns. Ignored if N is not dynamic
         */
        void Resize(size_t newM, size_t newN) {
            if (M != Dynamic && N != Dynamic)
                throw "Cannot resize statically sized matrices.";

            size_t oldM = this->m, oldN = this->n;
            if (M == Dynamic) { this->m = newM; }
            else { this->m = M; }
            if (N == Dynamic) { this->n = newN; }
            else { this->n = N; }

            if (this->m == 0 || this->n == 0) {
                this->m = 0;
                this->n = 0;
                if (this->data != NULL)
                    delete this->data;
                this->data = NULL;
                return;
            }

            Complex<T> * oldData = this->data;
            this->data = new Complex<T>[this->m*this->n];
            for (size_t r = 0; r < this->m; ++r) {
                for (size_t c = 0; c < this->n; ++c) {
                    if (r >= oldM || c >= oldN) {
                        if (Flags & ColumnMajor)
                            this->data[c*this->m+r] = T(0);
                        else
                            this->data[r*this->n+c] = T(0);
                    }
                    else if (oldData != NULL) {
                        if (Flags & ColumnMajor)
                            this->data[c*this->m+r] = oldData[c*oldM+r];
                        else
                            this->data[r*this->n+c] = oldData[r*oldN+c];
                    }
                }
            }
            if (oldData != NULL)
                delete oldData;
        }

        /**
         * @return Number of rows
         */
        size_t NumRows() const { return this->m; }
        /**
         * @return Number of columns
         */
        size_t NumColumns() const { return this->n; }
        /**
         * @return Number of cells (Equivalent to NumRows()*NumColumns())
         */
        size_t Size() const { return NumRows()*NumColumns(); }

        /**
         * Takes row r of the matrix and returns it. If \f$r\ge M\f$, an exception is thrown.
         * @param r Row index
         * @return 1xN Matrix
         */
        Matrix<T,1,N> GetRow(size_t r) const {
            if (r >= NumRows())
                throw "Cannot get row from matrix. Index out of bounds.";
            Matrix<T,1,N> ret(NumColumns(), T(0));
            for (size_t c = 0; c < NumColumns(); ++c)
                ret(0,c) = (*this)(r,c);
            return ret;
        }
        /**
         * Takes column c of the matrix and returns it. If \f$c\ge N\f$, an exception is thrown.
         * @param c Column index
         * @return Mx1 Matrix
         */
        Matrix<T,M,1> GetColumn(size_t c) const {
            if (c >= NumColumns())
                throw "Cannot get row from matrix. Index out of bounds.";
            Matrix<T,M,1> ret(NumRows(),T(0));
            for (size_t r = 0; r < NumRows(); ++r)
                ret(r,0) = (*this)(r,c);
            return ret;
        }
        /**
         * Sets row r of the matrix to row. If \f$r\ge M\f$, an exception is thrown. If \f$P\ne 1\f$ or \f$Q\ne N\f$, an exception is thrown.
         * @param r Row index
         * @param row PxQ Matrix
         */
        template <size_t P, size_t Q>
        typename std::enable_if<(P==1||P==Dynamic)&&(Q==N||Q==Dynamic||N==Dynamic),void>::type SetRow(size_t r, const Matrix<T,P,Q>& row) {
            if (r >= NumRows())
                throw "Cannot set row in matrix. Index out of bounds.";
            if (row.NumRows() != 1)
                throw "Expected a row vector in Matrix::SetRow()";
            if (row.NumColumns() != NumColumns())
                throw "Cannot set row in matrix. Size mismatch.";
            for (size_t c = 0; c < NumColumns(); ++c)
                (*this)(r,c) = row(0,c);
        }
        /**
         * Sets column c of the matrix to row. If \f$c\ge N\f$, an exception is thrown. If \f$P\ne M\f$ or \f$Q\ne 1\f$, an exception is thrown.
         * @param c Column index
         * @param column PxQ Matrix
         */
        template <size_t P, size_t Q>
        typename std::enable_if<(P==M||P==Dynamic||M==Dynamic)&&(Q==1||Q==Dynamic),void>::type SetColumn(size_t c, const Matrix<T,P,Q>& column) {
            if (c >= NumColumns())
                throw "Cannot set column in matrix. Index out of bounds.";
            if (column.NumColumns() != 1)
                throw "Expected a vector in Matrix::SetColumn()";
            if (column.NumRows() != NumRows())
                throw "Cannot set row in matrix. Size mismatch.";
            for (size_t r = 0; r < NumRows(); ++r)
                (*this)(r,c) = column(r,0);
        }
        /**
         * Swaps rows r1 and r2. If \f$r1\ge M\f$ or \f$r2\ge M\f$ an exception is thrown.
         * @param r1 Row 1 index
         * @param r2 Row 2 index
         */
        void SwapRows(size_t r1, size_t r2) {
            if (r1 >= NumRows() || r2 >= NumRows())
                throw "Cannot swap rows. Index out of bounds.";
            for (size_t c = 0; c < NumColumns(); ++c) {
                Complex<T> temp = (*this)(r1,c);
                (*this)(r1,c) = (*this)(r2,c);
                (*this)(r2,c) = temp;
            }
        }
        /**
         * Scales row r by s, i.e., \f$a_{ri}=sa_{ri}\f$ for \f$0\le i<N\f$. If \f$r\ge M\f$ an exception is thrown.
         * @param r Row index
         * @param s Complex number
         */
        void ScaleRow(size_t r, Complex<T> s) {
            if (r >= NumRows())
                throw "Cannot scale row. Index out of bounds.";
            for (size_t c = 0; c < NumColumns(); ++c)
                (*this)(r,c) *= s;
        }
        /**
         * Scales row r2 by s then adds it to r1, i.e., \f$a_{r1,i}=a_{r1,i}+sa_{r2,i}\f$ for \f$0\le i<N\f$. If \f$r1\ge M\f$ or \f$r2\ge M\f$ an exception is thrown.
         * @param r1 Row 1 index
         * @param r2 Row 2 index
         * @param s Complex number
         */
        void AddRows(size_t r1, size_t r2, Complex<T> s) {
            if (r1 >= NumRows() || r2 >= NumRows())
                throw "Cannot add rows. Index out of bounds.";
            for (size_t c = 0; c < NumColumns(); ++c)
                (*this)(r1,c) += s*(*this)(r2,c);
        }

        /// Operators.
        // Access operators
        /**
         * Returns the ith entry of the data. This accesses the data directly, the result will depend on if the matrix is row major or column major.
         * @param i Index
         * @return Complex number
         */
        Complex<T> operator[] (size_t i) const {
            if (i >= Size())
                throw "Cannot access matrix. Index out of bounds.";
            return this->data[i];
        }
        Complex<T> & operator[] (size_t i) {
            if (i >= Size())
                throw "Cannot access matrix. Index out of bounds.";
            return this->data[i];
        }
        /**
         * Returns the the data located in row r and column c. If the matrix is column major, r still represents the row and c
         * still represents the column.
         * @param r Row index
         * @param c Column index
         * @return Complex number
         */
        Complex<T> operator() (size_t r, size_t c) const {
            if (r >= NumRows() || c >= NumColumns())
                throw "Cannot access matrix. Row and column out of range.";
            if (Flags & ColumnMajor)
                return this->data[c*this->m+r];
            return this->data[r*this->n+c];
        }
        Complex<T> & operator() (size_t r, size_t c) {
            if (r >= NumRows() || c >= NumColumns())
                throw "Cannot access matrix. Row and column out of range.";
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
        /**
         * Set the matrix to other. If M or N are Dynamic, they are set to P and Q respectively. Otherwise, if
         * \f$M\ne P\f$ or \f$N\ne Q\f$ an exception is raised.
         * @param other PxQ Matrix
         * @return PxQ Matrix
         */
        template <size_t P, size_t Q, unsigned int Flags2>
        Matrix<T,M,N,Flags> & operator=(const Matrix<T,P,Q,Flags2>& other) {
            if (M == Dynamic && N == Dynamic) { Resize(other.NumRows(), other.NumColumns()); }
            else if (M == Dynamic) { Resize(other.NumRows()); }
            else if (N == Dynamic) { Resize(other.NumColumns()); }
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
        /**
         * Computes the MxN matrix \f$C=A+B\f$ defined by \f$c_{ij}=a_{ij}+b_{ij}\f$.
         * If \f$M\ne P\f$ or \f$N\ne Q\f$ an exception is raised.
         * @param A MxN Matrix
         * @param B PxQ Matrix
         * @return MxN Matrix
         */
        template <size_t P, size_t Q, unsigned int Flags2>
        friend Matrix<T,M,N,Flags> operator+(Matrix<T,M,N,Flags> A, const Matrix<T,P,Q,Flags2>& B) { return A += B; }
        /**
         * Computes the MxN matrix \f$C=A-B\f$ defined by \f$c_{ij}=a_{ij}-b_{ij}\f$.
         * If \f$M\ne P\f$ or \f$N\ne Q\f$ an exception is raised.
         * @param A MxN Matrix
         * @param B PxQ Matrix
         * @return MxN Matrix
         */
        template <size_t P, size_t Q, unsigned int Flags2>
        friend Matrix<T,M,N,Flags> operator-(Matrix<T,M,N,Flags> A, const Matrix<T,P,Q,Flags2>& B) { return A -= B; }
        /**
         * Computes the MxN matrix \f$B=sA\f$ defined by \f$b_{ij}=sa_{ij}\f$.
         * @param A MxN Matrix
         * @param s Complex number
         * @return MxN Matrix
         */
        friend Matrix<T,M,N,Flags> operator*(Matrix<T,M,N,Flags> A, const Complex<T>& s) { return A *= s; }
        /**
         * Computes the MxN matrix \f$B=sA\f$ defined by \f$b_{ij}=sa_{ij}\f$.
         * @param s Complex number
         * @param A MxN Matrix
         * @return MxN Matrix
         */
        friend Matrix<T,M,N,Flags> operator*(const Complex<T>& s, Matrix<T,M,N,Flags> A) { return A *= s; }
        /**
         * Computes the MxN matrix \f$B=A/s\f$ defined by \f$b_{ij}=a_{ij}/s\f$.
         * @param A MxN Matrix
         * @param s Complex number
         * @return MxN Matrix
         */
        friend Matrix<T,M,N,Flags> operator/(Matrix<T,M,N,Flags> A, const Complex<T>& s) { return A /= s; }
        /**
         * Computes the MxQ matrix \f$C=AB\f$ defined by \f$c_{ij}=\sum_{k=0}^{N-1}a_{ik}b_{kj}\f$.
         * If \f$N\ne P\f$ an exception is raised.
         * @param A MxN Matrix
         * @param B PxQ Matrix
         * @return MxQ Matrix
         */
        template <size_t P, size_t Q, unsigned int Flags2>
        friend Matrix<T,M,Q,Flags> operator*(const Matrix<T,M,N,Flags>& A, const Matrix<T,P,Q,Flags2>& B) {
            if (A.NumColumns() != B.NumRows())
                throw "Cannot muliply two matrices due to size mismatch.";
            Matrix<T,M,Q,Flags> ret(A.NumRows(), B.NumColumns(), T(0));
            for (size_t r = 0; r < A.NumRows(); ++r) {
                for (size_t c = 0; c < B.NumColumns(); ++c) {
                    for (size_t i = 0; i < A.NumColumns(); ++i) {
                        ret(r,c) += A(r,i)*B(i,c);
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
        /**
         * Takes a matrix and writes to out in a formatted manner.
         * @param out Reference to std::ostream.
         * @param m MxN Matrix.
         */
        friend std::ostream& operator<<(std::ostream& out, const Matrix<T,M,N,Flags>& m) {
            if (m.NumRows() == 1) {
                out << "(1x"<<m.NumColumns()<<")[";
                for (size_t i = 0; i < m.Size(); ++i)
                    out << (i > 0 ? ", " : "") << m[i];
                out << "]";
                return out;
            }

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
        /**
         * Compares two matrices entrywise.
         * @param a MxN Matrix
         * @param b PxQ Matrix
         * @return If \f$M\ne P\f$ or \f$N\ne Q\f$ it will return false. It returns true if and only if each entry is equal.
         */
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
        /**
         * Compares two matrices entrywise.
         * @param a MxN Matrix
         * @param b PxQ Matrix
         * @return If \f$M\ne P\f$ or \f$N\ne Q\f$ it will return true. It returns false if and only if each entry is equal.
         */
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
