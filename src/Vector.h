#pragma once
#include "Matrix.h"

namespace Linear {
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
