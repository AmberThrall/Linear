#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    Matrix3f A = {
        {1,0,0}, {2,3,0}, {4,5,6}
    };
    Vector3f a = { 1, 1, 1 };
    Vector3f b = { 1, 2, 3 };
    RowVector3f c = Transpose(a);
    RowVector3f d = Transpose(b);

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "a+b = " << a+b << std::endl;
    std::cout << "Ab = " << A*b << std::endl;

    std::cout << "c = " << c << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "c+d = " << c+d << std::endl;
    std::cout << "dA = " << d*A << std::endl;

    std::cout << "len(a) = " << Length(a) << std::endl;
    std::cout << "len(b) = " << Length(b) << std::endl;
    std::cout << "dot(a,b) = " << Dot(a,b) << std::endl;
    std::cout << "a x b = " << Cross(a,b) << std::endl;

    std::cout << "Normalize(c) = " << Normalize(c) << std::endl;

    Matrix3f gramschmidt = { {1,0,0}, {2,3,0}, {4,5,6} };
    Matrix3f orthogonal = GramSchmidt(gramschmidt);
    std::cout << "gramschmidt before = " << gramschmidt << std::endl;
    std::cout << "gramschmidt after = " << orthogonal << std::endl;
    std::cout << "gramschmidt test = " << Transpose(orthogonal)*orthogonal << std::endl;
}
