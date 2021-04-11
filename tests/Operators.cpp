#include "../src/Matrix.h"
#include <iostream>
using namespace Linear;

void test(Matrix3f a) {
    a *= -1;
}

int main() {
    Matrix3f a(1);
    Matrix3f b = {
        {1,2,3}, {4,5,6}, {7,8,9}
    };

    test(a);
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "a+b = " << a+b << std::endl;
    std::cout << "a-b = " << a-b << std::endl;
    std::cout << "a*2 = " << a*2 << std::endl;
    std::cout << "a/2 = " << a/2 << std::endl;
    std::cout << "-a = " << -a << std::endl;
    std::cout << "a*b = " << a*b << std::endl;
    std::cout << "b*a = " << b*a << std::endl;
    std::cout << "a==b = " << (a==b ? "true" : "false") << std::endl;
    std::cout << "a!=b = " << (a!=b ? "true" : "false") << std::endl;
}
