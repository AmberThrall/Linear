#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    Matrix3f a = {
        {1,0,0}, {2,3,0}, {4,5,6}
    };
    std::cout << "a = " << a << std::endl;
    std::cout << "a^T = " << Transpose(a) << std::endl;
    std::cout << "tr(a) = " << Trace(a) << std::endl;
    std::cout << "det(a) = " << Determinant(a) << std::endl;
    std::cout << "a^{-1} = " << Inverse(a) << std::endl;
    std::cout << "a*a^{-1} = " << a*Inverse(a) << std::endl;

    Matrix2f b = SubMatrix<2,2>(a);
    std::cout << "b = " << b << std::endl;
    std::cout << "b^T = " << Transpose(b) << std::endl;
    std::cout << "tr(b) = " << Trace(b) << std::endl;
    std::cout << "det(b) = " << Determinant(b) << std::endl;
    std::cout << "b^{-1} = " << Inverse(b) << std::endl;
    std::cout << "b*b^{-1} = " << b*Inverse(b) << std::endl;

    Matrix<float,3,6> c = {
        {1,0,0,1,0,0}, {2,3,0,0,1,0}, {4,5,6,0,0,1}
    };
    std::cout << "c = " << c << std::endl;
    std::cout << "rref(c) = " << RREF(c) << std::endl;
}
