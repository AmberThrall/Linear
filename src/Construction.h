#pragma once
#include <vector>
#include <array>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"

namespace Linear {
    template <typename T, size_t N>
    typename std::enable_if<(N>0), SquareMatrix<T,N>>::type Identity() {
        SquareMatrix<T,N> ret(T(0));
        for (size_t i = 0; i < N; ++i)
            ret(i,i) = T(1);
        return ret;
    }
    template<typename T>
    SquareMatrix<T,Dynamic> Identity(size_t n) {
        SquareMatrix<T,Dynamic> ret(n, n, T(0));
        for (size_t i = 0; i < n; ++i)
            ret(i,i) = T(1);
        return ret;
    }

    template<typename T, typename... Ts>
    SquareMatrix<T,sizeof...(Ts)+1> Diag(Complex<T> first, Ts... ts) {
        SquareMatrix<T,sizeof...(Ts)+1> ret(sizeof...(Ts)+1, sizeof...(Ts)+1, T(0));
        std::array<Complex<T>,sizeof...(Ts)> arr = {ts...};
        ret(0,0) = first;
        for (size_t i = 1; i < sizeof...(Ts); ++i)
            ret(i,i) = arr[i];
        return ret;
    }
    template<typename T>
    SquareMatrix<T,Dynamic> Diag(std::vector<Complex<T>> v) {
        SquareMatrix<T,Dynamic> ret(v.size(), v.size(), T(0));
        for (size_t i = 0; i < v.size(); ++i)
            ret(i,i) = v[i];
        return ret;
    }
    template<typename T, size_t N>
    SquareMatrix<T,N> Diag(Vector<T,N> v) {
        SquareMatrix<T,N> ret(v.NumRows(), v.NumRows(), T(0));
        for (size_t i = 0; i < v.NumRows(); ++i)
            ret(i,i) = v(i,0);
        return ret;
    }
    template<typename T, size_t N>
    SquareMatrix<T,N> Diag(RowVector<T,N> v) {
        return Diag(Transpose(v));
    }
}
