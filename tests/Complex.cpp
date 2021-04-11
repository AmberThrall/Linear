#include "../src/Complex.h"
#include <iostream>
using namespace Linear;

int main() {
    Complexf z(1,1);
    Complexf w(2,-3);

    std::cout << "z = " << z << std::endl;
    std::cout << "w = " << w << std::endl;
    std::cout << "z+w = " << z+w << std::endl;
    std::cout << "z+2 = " << z+2 << std::endl;
    std::cout << "z-w = " << z-w << std::endl;
    std::cout << "z-2 = " << z-2 << std::endl;
    std::cout << "z*w = " << z*w << std::endl;
    std::cout << "z*2 = " << z*2 << std::endl;
    std::cout << "z/w = " << z/w << std::endl;
    std::cout << "w/z = " << w/z << std::endl;
    std::cout << "z/2 = " << z/2 << std::endl;
    std::cout << "2/z = " << 2/z << std::endl;
    std::cout << "\\bar{z} = " << z.Conjugate() << std::endl;
    std::cout << "sqrt(z) = " << Complexf::Sqrt(z) << std::endl;
    std::cout << "exp(z) = " << Complexf::Exp(z) << std::endl;
    std::cout << "log(z) = " << Complexf::Log(z) << std::endl;
    std::cout << "pow(z,w) = " << Complexf::Pow(z,w) << std::endl;
    std::cout << "pow(w,z) = " << Complexf::Pow(w,z) << std::endl;
    std::cout << "z/0 = " << z/0 << std::endl;
}
