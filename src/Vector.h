#pragma once
#include "Matrix.h"

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

    template<typename T, size_t P, size_t Q>
    typename std::enable_if<(P==Q||P==Dynamic||Q==Dynamic),Complex<T>>::type Dot(const Vector<T,P>& a, const Vector<T,Q>& b) {
        if (a.Size() != b.Size())
            throw "Cannot take the dot product of two different sized vectors.";
        Complex<T> ret;
        for (size_t i = 0; i < a.Size(); ++i) {
            ret += a[i].Conjugate()*b[i];
        }
        return ret;
    }
    template<typename T, size_t P, size_t Q>
    typename std::enable_if<(P==Q||P==Dynamic||Q==Dynamic),Complex<T>>::type Dot(const Vector<T,P>& a, const RowVector<T,Q>& b) {
        return Dot(a, Transpose(b));
    }
    template<typename T, size_t P, size_t Q>
    typename std::enable_if<(P==Q||P==Dynamic||Q==Dynamic),Complex<T>>::type Dot(const RowVector<T,P>& a, const Vector<T,Q>& b) {
        return Dot(Transpose(a), b);
    }
    template<typename T, size_t P, size_t Q>
    typename std::enable_if<(P==Q||P==Dynamic||Q==Dynamic),Complex<T>>::type Dot(const RowVector<T,P>& a, const RowVector<T,Q>& b) {
        return Dot(Transpose(a), Transpose(b));
    }

    template <typename T, size_t N>
    T SquaredLength(Vector<T,N> v) {
        return Dot(v,v).Re;
    }
    template <typename T, size_t N>
    T SquaredLength(RowVector<T,N> v) {
        return Dot(v,v).Re;
    }
    template <typename T, size_t N>
    T Length(Vector<T,N> v) {
        return std::sqrt(SquaredLength(v));
    }
    template <typename T, size_t N>
    T Length(RowVector<T,N> v) {
        return std::sqrt(SquaredLength(v));
    }

    template <typename T, size_t N>
    Vector<T,N> Normalize(Vector<T,N> v) {
        return v / Length(v);
    }
    template <typename T, size_t N>
    RowVector<T,N> Normalize(RowVector<T,N> v) {
        return v / Length(v);
    }

    template<typename T, size_t P, size_t Q>
    typename std::enable_if<(P==3||P==Dynamic)&&(Q==3||Q==Dynamic),Vector<T,3>>::type Cross(Vector<T,P> a, Vector<T,Q> b) {
        if (a.Size() != 3 || b.Size() != 3)
            throw "The cross product is only defined for the 3-dimensional vectors.";
        Vector<T,3> ret(T(0));
        ret[0] = a[1]*b[2] - a[2]*b[1];
        ret[1] = a[2]*b[0] - a[0]*b[2];
        ret[2] = a[0]*b[1] - a[1]*b[0];
        return ret;
    }
    template<typename T, size_t P, size_t Q>
    typename std::enable_if<(P==3||P==Dynamic)&&(Q==3||Q==Dynamic),Vector<T,3>>::type Cross(Vector<T,P> a, RowVector<T,Q> b) {
        return Cross(a, Transpose(b));
    }
    template<typename T, size_t P, size_t Q>
    typename std::enable_if<(P==3||P==Dynamic)&&(Q==3||Q==Dynamic),RowVector<T,3>>::type Cross(RowVector<T,P> a, Vector<T,Q> b) {
        return Transpose(Cross(Transpose(a), b));
    }
    template<typename T, size_t P, size_t Q>
    typename std::enable_if<(P==3||P==Dynamic)&&(Q==3||Q==Dynamic),RowVector<T,3>>::type Cross(RowVector<T,P> a, RowVector<T,Q> b) {
        return Transpose(Cross(Transpose(a), Transpose(b)));
    }
}
