#pragma once
#include <vector>
#include <iostream>
#include "Matrix.h"
#include "Basics.h"
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

    template<typename T, size_t M, size_t N, unsigned int Flags>
    bool IsVector(const Matrix<T,M,N,Flags>& m) {
        return (m.NumColumns() == 1 || m.NumRows() == 1);
    }

    template <size_t P, typename T, size_t N>
    typename std::enable_if<(P>0),Vector<T,P>>::type SubVector(const Vector<T,N>& v, size_t off = 0) {
        if (off+P > v.Size())
            throw "Cannot create subvector, indices out of bounds.";
        Vector<T,P> ret(T(0));
        for (size_t i = 0; i < ret.Size(); ++i) {
            ret[i] = v[off+i];
        }
        return ret;
    }
    template <typename T, size_t N>
    Vector<T,Dynamic> SubVector(const Vector<T,N>& v, size_t size, size_t off = 0) {
        if (off+size > v.Size())
            throw "Cannot create subvector, indices out of bounds.";
        Vector<T,Dynamic> ret(size,T(0));
        for (size_t i = 0; i < ret.Size(); ++i) {
            ret[i] = v[off+i];
        }
        return ret;
    }
    template <size_t P, typename T, size_t N>
    typename std::enable_if<(P>0),RowVector<T,P>>::type SubVector(const RowVector<T,N>& v, size_t off = 0) {
        return Transpose(SubVector<T,P>(Transpose(v), off));
    }
    template <typename T, size_t N>
    RowVector<T,Dynamic> SubVector(const RowVector<T,N>& v, size_t size, size_t off = 0) {
        return Transpose(SubVector<T>(Transpose(v), size, off));
    }
    template <size_t P, typename T>
    typename std::enable_if<(P==1),Vector<T,1>>::type SubVector(const RowVector<T,1>& v, size_t off = 0) {
        if (off != 0)
            throw "Cannot create subvector, indices out of bounds.";
        return v;
    }
    template <typename T>
    Vector<T,Dynamic> SubVector(const RowVector<T,1>& v, size_t size, size_t off = 0) {
        if (off != 0 || size != 1)
            throw "Cannot create subvector, indices out of bounds.";
        return v;
    }

    template<typename T,size_t N>
    typename std::enable_if<(N>0), Vector<T,N>>::type Basis(size_t i) {
        Vector<T,N> ret(T(0));
        ret[i] = T(1);
        return ret;
    }
    template<typename T>
    Vector<T,Dynamic> Basis(size_t dim, size_t i) {
        Vector<T,Dynamic> ret(dim, T(0));
        ret[i] = T(1);
        return ret;
    }

    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    typename std::enable_if<((M==1||N==1||M==Dynamic||N==Dynamic)&&(P==1||Q==1||P==Dynamic||Q==Dynamic)),Complex<T>>::type
    Dot(const Matrix<T,M,N,Flags>& a, const Matrix<T,P,Q,Flags2>& b) {
        if (!IsVector(a) || !IsVector(b))
            throw "Dot product is only defined for vectors.";
        if (a.Size() != b.Size())
            throw "Cannot take the dot product of two different sized vectors.";
        Complex<T> ret;
        for (size_t i = 0; i < a.Size(); ++i) {
            ret += Conjugate(a[i])*b[i];
        }
        return ret;
    }

    template<typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==1||N==1||M==Dynamic||N==Dynamic), T>::type SquaredLength(const Matrix<T,M,N,Flags>& v) {
        if (!IsVector(v))
            throw "Cannot calculate the length of a matrix.";
        return Dot(v,v).Re;
    }
    template<typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==1||N==1||M==Dynamic||N==Dynamic), T>::type Length(const Matrix<T,M,N,Flags>& v) {
        return Sqrt(SquaredLength(v)).Re;
    }

    template <typename T, size_t N>
    Vector<T,N> Normalize(Vector<T,N> v) {
        T len = Length(v);
        if (len < T(Tol)) {
            for (size_t i = 0; i < v.Size(); ++i)
                v[i] = 0;
            return v;
        }
        return v / len;
    }
    template <typename T, size_t N>
    RowVector<T,N> Normalize(RowVector<T,N> v) {
        return Transpose(Normalize(Transpose(v)));
    }
    template <typename T>
    Vector<T,1> Normalize(Vector<T,1> v) {
        v[0] = T(1);
        return v;
    }

    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    typename std::enable_if<((M==1||N==1||M==Dynamic||N==Dynamic)&&(P==1||Q==1||P==Dynamic||Q==Dynamic)),Matrix<T,M,N,Flags>>::type
    Cross(const Matrix<T,M,N,Flags>& a, const Matrix<T,P,Q,Flags2>& b) {
        if (!IsVector(a) || !IsVector(b))
            throw "Cross product is only defined for vectors.";
        if (a.Size() != 3 || b.Size() != 3)
            throw "The cross product is only defined for the 3-dimensional vectors.";
        Matrix<T,M,N,Flags> ret(a.NumRows(), a.NumColumns(),T(0));
        ret[0] = a[1]*b[2] - a[2]*b[1];
        ret[1] = a[2]*b[0] - a[0]*b[2];
        ret[2] = a[0]*b[1] - a[1]*b[0];
        return ret;
    }

    template<typename T, size_t M, size_t N, unsigned int Flags, size_t P, size_t Q, unsigned int Flags2>
    typename std::enable_if<((M==1||N==1||M==Dynamic||N==Dynamic)&&(P==1||Q==1||P==Dynamic||Q==Dynamic)),Matrix<T,P,Q,Flags2>>::type
    Proj(const Matrix<T,M,N,Flags>& v, const Matrix<T,P,Q,Flags2>& onto) {
        if (!IsVector(v) || !IsVector(onto))
            throw "Vector projection doesn't work on matrices.";
        if (v.Size() != onto.Size())
            throw "Cannot project vector onto a vector of differing dimension.";
        if (Length(onto) < T(Tol)) { // Special case: Proj_0(v) = 0
            for (size_t i = 0; i < onto.Size(); ++i)
                onto[i] = T(0);
            return onto;
        }
        return (Dot(onto,v)/Dot(onto,onto))*onto;
    }

    template <typename T, size_t N>
    std::vector<Vector<T,N>> GramSchmidt(const std::vector<Vector<T,N>>& v) {
        if (v.size() == 0)
            return v;
        std::vector<Vector<T,N>> u;
        for (size_t k = 0; k < v.size(); ++k) {
            u.push_back(v[k]);
            for (size_t i = 0; i < k; ++i) {
                u[k] -= Proj(u[k], u[i]);
            }
            T len = Length(u[k]);
            if (len >= T(Tol))
                u[k] = u[k]/len;
        }
        return u;
    }
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
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Matrix<T,M,N,Flags> GramSchmidt(Matrix<T,M,N,Flags> m) {
        std::vector<Vector<T,M>> columns;
        for (size_t i = 0; i < m.NumColumns(); ++i)
            columns.push_back(m.GetColumn(i));
        columns = GramSchmidt(columns);
        for (size_t i = 0; i < m.NumColumns(); ++i)
            m.SetColumn(i, columns[i]);
        return m;
    }
    template <typename T>
    std::vector<Vector<T,1>> GramSchmidt(const std::vector<Vector<T,1>>& v) {
        if (v.size() == 0)
            return v;
        std::vector<Vector<T,1>> res;
        res.push_back(Normalize(v[0]));
        return res;
    }
}
