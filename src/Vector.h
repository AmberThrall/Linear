#pragma once
#include <vector>
#include <iostream>
#include "Matrix.h"
#include "Types.h"
#include "Global.h"

namespace Linear {
    template<typename T, size_t N>
    using Vector = Matrix<T,N,1>;
    using Vector2i = Vector<int,2>;
    using Vector3i = Vector<int,3>;
    using Vector4i = Vector<int,4>;
    using VectorXi = Vector<int,Dynamic>;
    using Vector2f = Vector<float,2>;
    using Vector3f = Vector<float,3>;
    using Vector4f = Vector<float,4>;
    using VectorXf = Vector<float,Dynamic>;
    using Vector2d = Vector<double,2>;
    using Vector3d = Vector<double,3>;
    using Vector4d = Vector<double,4>;
    using VectorXd = Vector<double,Dynamic>;

    template<typename T, size_t N>
    using RowVector = Matrix<T,1,N>;
    using RowVector2i = RowVector<int,2>;
    using RowVector3i = RowVector<int,3>;
    using RowVector4i = RowVector<int,4>;
    using RowVectorXi = RowVector<int,Dynamic>;
    using RowVector2f = RowVector<float,2>;
    using RowVector3f = RowVector<float,3>;
    using RowVector4f = RowVector<float,4>;
    using RowVectorXf = RowVector<float,Dynamic>;
    using RowVector2d = RowVector<double,2>;
    using RowVector3d = RowVector<double,3>;
    using RowVector4d = RowVector<double,4>;
    using RowVectorXd = RowVector<double,Dynamic>;

    /**
     * Constructs a subvector of a column vector. If P<1 or off+P>N, an exception is thrown.
     *
     * Example: w is column vector {2,3}
     *
     *     Vector4d v = { 1, 2, 3, 4 };
     *     Vector2d w = SubVector<2>(v, 1);
     *
     * @param P Length for subvector
     * @param v Column vector of size N
     * @param off Offset to start the subvector (default = 0)
     * @return Column vector of size P
     */
    template <size_t P, typename T, size_t N>
    Vector<T,P> SubVector(const Vector<T,N>& v, size_t off = 0) {
        if (off+P > v.Length())
            throw "Cannot create subvector, indices out of bounds.";
        Vector<T,P> ret(T(0));
        for (size_t i = 0; i < ret.Length(); ++i) {
            ret[i] = v[off+i];
        }
        return ret;
    }
    /**
     * Constructs a subvector of a column vector. If size<1 or off+size>N, an exception is thrown.
     *
     * Example: w is dynamic column vector {2,3}
     *
     *     Vector4d v = { 1, 2, 3, 4 };
     *     VectorXd w = SubVector(v, 2, 1);
     *
     * @param v Column vector of size N
     * @param size Length for subvector
     * @param off Offset to start the subvector (default 0)
     * @return Column vector of size size
     */
    template <typename T, size_t N>
    Vector<T,Dynamic> SubVector(const Vector<T,N>& v, size_t size, size_t off = 0) {
        if (off+size > v.Length())
            throw "Cannot create subvector, indices out of bounds.";
        Vector<T,Dynamic> ret(size,T(0));
        for (size_t i = 0; i < ret.Length(); ++i) {
            ret[i] = v[off+i];
        }
        return ret;
    }
    /**
     * Constructs a subvector of a row vector. If P<1 or off+P>N, an exception is thrown.
     *
     * Example: w is row vector {2,3}
     *
     *     RowVector4d v = { 1, 2, 3, 4 };
     *     RowVector2d w = SubVector<2>(v, 1);
     *
     * @param P Length for subvector
     * @param v Row vector of size N
     * @param off Offset to start the subvector (default = 0)
     * @return Row vector of size P
     */
    template <size_t P, typename T, size_t N>
    RowVector<T,P> SubVector(const RowVector<T,N>& v, size_t off = 0) {
        return Transpose(SubVector<T,P>(Transpose(v), off));
    }
    /**
     * Constructs a subvector of a row vector. If size<1 or off+size>N, an exception is thrown.
     *
     * Example: w is dynamic row vector {2,3}
     *
     *     RowVector4d v = { 1, 2, 3, 4 };
     *     RowVectorXd w = SubVector<2>(v, 2, 1);
     *
     * @param v Row vector of size N
     * @param size Length for subvector
     * @param off Offset to start the subvector (default = 0)
     * @return Row vector of size size
     */
    template <typename T, size_t N>
    RowVector<T,Dynamic> SubVector(const RowVector<T,N>& v, size_t size, size_t off = 0) {
        return Transpose(SubVector<T>(Transpose(v), size, off));
    }
    /**
     * Handles subvector in the size 1 case. If P is not one or off is not zero, an exception is raised.
     *
     * Example: w is vector {1}
     *
     *     Vector<double,1> v = { 1 };
     *     Vector<double,1> w = SubVector<1>(v);
     *
     * @param P Length for subvector
     * @param v Vector of length 1
     * @param off Offset to start the subvector (default = 0)
     * @return v
     */
    template <size_t P, typename T>
    Vector<T,P> SubVector(const Vector<T,1>& v, size_t off = 0) {
        if (off != 0)
            throw "Cannot create subvector, indices out of bounds.";
        Vector<T,P> ret(v[0]);
        return ret;
    }
    /**
     * Handles subvector in the size 1 case. If size is not one or off is not zero, an exception is raised.
     *
     * Example: w is dynamic vector {1}
     *
     *     Vector<double,1> v = { 1 };
     *     VectorXd w = SubVector(v, 1);
     *
     * @param v Row vector of size 1
     * @param size Size for subvector
     * @param off Offset to start the subvector (default = 0)
     * @return v
     */
    template <typename T>
    Vector<T,Dynamic> SubVector(const Vector<T,1>& v, size_t size, size_t off = 0) {
        if (off != 0)
            throw "Cannot create subvector, indices out of bounds.";
        Vector<T,Dynamic> ret(size, v[0]);
        return ret;
    }

    /**
     * Computes the dot product of two vectors a and b both of size N, \f$\sum{i=0}^{N-1}\overline{a_i}b_i\f$. If either a or b is not a vector,
     * or they have differing sizes, an exception is thrown.
     *
     * Example: z=1*3+2*4=11
     *
     *     Vector2d a = {1, 2};
     *     Vector2d b = {3, 4};
     *     Complexd z = Dot(a,b);
     *
     * @param a Row/column vector
     * @param b Row/column vector
     * @return Complex number
     */
    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    Complex<T> Dot(const Matrix<T,M,N,Flags>& a, const Matrix<T,P,Q,Flags2>& b) {
        if (!IsVector(a) || !IsVector(b))
            throw "Dot product is only defined for vectors.";
        if (a.Length() != b.Length())
            throw "Cannot take the dot product of two different sized vectors.";
        Complex<T> ret;
        for (size_t i = 0; i < a.Length(); ++i) {
            ret += Conjugate(a[i])*b[i];
        }
        return ret;
    }

    /**
     * Normalize a vector, \f$v=\frac{v}{\|v\|}\f$.
     *
     * Example: w is the column vector {1/5, 2/5}
     *
     *     Vector2d v = {1, 2};
     *     Vector2d w = Normalize(v);
     *
     * @param v Row/column vector
     * @return Column vector if v is a column vector, otherwise row vector.
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> Normalize(Matrix<T,M,N,Flags> v) {
        if (!IsVector(v))
            throw "Normalize is only defined for vectors.";
        if (v.Length() == 0)
            return v;

        T norm = Norm(v);
        if (norm < T(Tol)) {
            for (size_t i = 0; i < v.Length(); ++i)
                v[i] = 0;
            return v;
        }
        return v / norm;
    }

    /**
     * Computes the cross product of two vectors a and b both of length 3. The cross product is defined by
     * \f\[a\times b = (a_0i+a_1j+a_2k)\times (b_0i+b_1j+b_2k)=(a_1b_2-a_2b_1)i + (a_2b_0-a_0b_2)j + (a_0b_1-a_1b_0)k.\f\]
     * If either a or b is not a vector, or they aren't length 3, an exception is thrown.
     *
     * Example: c is the column vector {-3,6,-3}
     *
     *     Vector3d a = {1, 2, 3};
     *     Vector3d b = {4, 5, 6};
     *     Vector3d c = Cross(a,b);
     *
     * @param a Row/column vector
     * @param b Row/column vector
     * @return Column vector if a is a column vector, otherwise row vector.
     */
    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    Matrix<T,M,N,Flags> Cross(const Matrix<T,M,N,Flags>& a, const Matrix<T,P,Q,Flags2>& b) {
        if (!IsVector(a) || !IsVector(b) || a.Length() != 3 || b.Length() != 3)
            throw "The cross product is only defined for two 3-dimensional vectors.";
        Matrix<T,M,N,Flags> ret(a.NumRows(), a.NumColumns(), T(0));
        ret[0] = a[1]*b[2] - a[2]*b[1];
        ret[1] = a[2]*b[0] - a[0]*b[2];
        ret[2] = a[0]*b[1] - a[1]*b[0];
        return ret;
    }

    /**
     * Computes the vector projection \f$\proj_{onto}v = \frac{Dot(onto,v)}{Dot(onto,onto)}onto\f$. If onto is the zero vector,
     * by definition this function will return the zero vector.
     * If either v or onto is not a vector, or they have different lengths, an exception is thrown.
     *
     * Example: w is the column vector {2,0}
     *
     *     Vector2d v = {2, 2};
     *     Vector2d onto = {1, 0};
     *     Vector2d w = Proj(v,onto);
     *
     * @param v Row/column vector
     * @param onto Row/column vector
     * @return Column vector if onto is a column vector, otherwise row vector.
     */
    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    Matrix<T,P,Q,Flags2> Proj(const Matrix<T,M,N,Flags>& v, const Matrix<T,P,Q,Flags2>& onto) {
        if (!IsVector(v) || !IsVector(onto))
            throw "Vector projection doesn't work on matrices.";
        if (v.Length() != onto.Length())
            throw "Cannot project vector onto a vector of differing dimension.";
        if (v.Length() == 0)
            throw "Cannot project empty vector.";
        if (Norm(onto) < T(Tol)) { // Special case: Proj_0(v) = 0
            for (size_t i = 0; i < onto.Length(); ++i)
                onto[i] = T(0);
            return onto;
        }
        return (Dot(onto,v)/Dot(onto,onto))*onto;
    }

    /**
     * Performs Gram-Schmidt on a list of column vectors.
     * If all vectors don't have the same length, an exception is thrown.
     *
     * Example: vectors is the list of column vectors {0.477,0.894} and {0.894,-0.477}
     *
     *     Vector2d v1 = { 1, 2 };
     *     Vector2d v2 = { 3, 4 };
     *     std::vector<Vector2d> vectors = { v1, v2 };
     *     vectors = GramSchmidt(vectors);
     *
     * @param v List of column vectors
     * @return List of column vectors
     */
    template <typename T, size_t N>
    std::vector<Vector<T,N>> GramSchmidt(const std::vector<Vector<T,N>>& v) {
        if (v.size() == 0)
            return v;
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i].Length() != v[0].Length())
                throw "All vectors must be the same length in order to perform Gram-Schmidt.";
        }
        if (v[0].Length() == 0)
            return v;

        std::vector<Vector<T,N>> u;
        for (size_t k = 0; k < v.size(); ++k) {
            u.push_back(v[k]);
            for (size_t i = 0; i < k; ++i) {
                u[k] -= Proj(u[k], u[i]);
            }
            u[k] = Normalize(u[k]);
        }
        return u;
    }
    /**
     * Performs Gram-Schmidt on a list of row vectors.
     * If all vectors don't have the same length, an exception is thrown.
     *
     * Example: vectors is the list of row vectors {0.477,0.894} and {0.894,-0.477}
     *
     *     RowVector2d v1 = { 1, 2 };
     *     RowVector2d v2 = { 3, 4 };
     *     std::vector<RowVector2d> vectors = { v1, v2 };
     *     vectors = GramSchmidt(vectors);
     *
     * @param v List of row vectors
     * @return List of row vectors
     */
    template <typename T, size_t N>
    std::vector<RowVector<T,N>> GramSchmidt(const std::vector<RowVector<T,N>>& v) {
        if (v.size() == 0)
            return v;

        std::vector<Vector<T,N>> vT;
        for (size_t i = 0; i < v.size(); ++i)
            vT.push_back(Transpose(v[i]));
        vT = GramSchmidt(vT);
        std::vector<RowVector<T,N>> res;
        for (size_t i = 0; i < v.size(); ++i)
            res.push_back(Transpose(vT[i]));
        return res;
    }
    /**
     * Performs Gram-Schmidt on the columns of A.
     *
     * Example: B is the 2x2 matrix {{0.477,0.894}, {0.894,-0.477}}
     *
     *     Matrix2d A = { {1, 3}, {2, 4} };
     *     Matrix2d B = GramSchmidt(A);
     *
     * @param A MxN Matrix
     * @return MxN matrix
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> GramSchmidt(Matrix<T,M,N,Flags> A) {
        if (A.NumEntries() == 0)
            return A;

        std::vector<Vector<T,M>> columns;
        for (size_t i = 0; i < A.NumColumns(); ++i)
            columns.push_back(A.GetColumn(i));
        columns = GramSchmidt(columns);
        for (size_t i = 0; i < A.NumColumns(); ++i)
            A.SetColumn(i, columns[i]);
        return A;
    }
    /**
     * Performs Gram-Schmidt on a list of length 1 vectors.
     *
     * Example: vectors is the list of length 1 vectors {1} and {1}
     *
     *     Vector<double, 1> v1 = { 1 };
     *     Vector<double, 1> v2 = { 2 };
     *     std::vector<Vector<double,1>> vectors = { v1, v2 };
     *     vectors = GramSchmidt(vectors);
     *
     * @param v List of vectors all of length 1
     * @return List of vectors all of length 1
     */
    template <typename T>
    std::vector<Vector<T,1>> GramSchmidt(const std::vector<Vector<T,1>>& v) {
        if (v.size() == 0)
            return v;
        
        std::vector<Vector<T,1>> res;
        for (size_t i = 0; i < v.size(); ++i) {
            res.push_back(Normalize(v[i]));
        }
        return res;
    }
}
