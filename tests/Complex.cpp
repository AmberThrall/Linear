#include "../src/Complex.h"
#include <iostream>
using namespace Linear;

int main() {
    Complexf z(1,1);
    Complexf w(2,-3);

    std::cout << "z = " << z << std::endl;
    std::cout << "w = " << w << std::endl;
    std::cout << "-z = " << -z << std::endl;
    std::cout << "z+w = " << z+w << std::endl;
    std::cout << "z+2 = " << z+2 << std::endl;
    std::cout << "z-w = " << z-w << std::endl;
    std::cout << "z-2 = " << z-2 << std::endl;
    std::cout << "w*i = " << w*Complexf::i() << std::endl;
    std::cout << "z*w = " << z*w << std::endl;
    std::cout << "z*2 = " << z*2 << std::endl;
    std::cout << "z/w = " << z/w << std::endl;
    std::cout << "w/z = " << w/z << std::endl;
    std::cout << "z/2 = " << z/2 << std::endl;
    std::cout << "2/z = " << 2/z << std::endl;
    std::cout << "\\bar{z} = " << Conjugate(z) << std::endl;
    std::cout << "sqrt(z) = " << Sqrt(z) << std::endl;
    std::cout << "exp(z) = " << Exp(z) << std::endl;
    std::cout << "log(z) = " << Log(z) << std::endl;
    std::cout << "log(e,i) = " << Log(Exp(Complexf(1,0)), Complexf::i()) << std::endl;
    std::cout << "pow(z,w) = " << Pow(z,w) << std::endl;
    std::cout << "pow(w,z) = " << Pow(w,z) << std::endl;

    std::cout << std::endl << "Trig. Functions: " << std::endl;
    std::cout << "asin(sin(z)) = " << ASin(Sin(z)) << std::endl;
    std::cout << "acos(cos(z)) = " << ACos(Cos(z)) << std::endl;
    std::cout << "atan(tan(z)) = " << ATan(Tan(z)) << std::endl;
    std::cout << "acsc(csc(z)) = " << ACsc(Csc(z)) << std::endl;
    std::cout << "asec(sec(z)) = " << ASec(Sec(z)) << std::endl;
    std::cout << "acot(cot(z)) = " << ACot(Cot(z)) << std::endl;

    std::cout << std::endl << "Hyperbolic Trig. Functions: " << std::endl;
    std::cout << "asinh(sinh(z)) = " << ASinh(Sinh(z)) << std::endl;
    std::cout << "acosh(cosh(z)) = " << ACosh(Cosh(z)) << std::endl;
    std::cout << "atanh(tanh(z)) = " << ATanh(Tanh(z)) << std::endl;
    std::cout << "acsch(csch(z)) = " << ACsch(Csch(z)) << std::endl;
    std::cout << "asech(sech(z)) = " << ASech(Sech(z)) << std::endl;
    std::cout << "acoth(coth(z)) = " << ACoth(Coth(z)) << std::endl;

    Complexi x = z+Complexf(0.5,0.5);
    std::cout << "x = " << x << std::endl;
}
