#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
        Matrix3f a = {
            {1,0,0}, {2,3,0}, {4,5,6}
        };
        Vector3f v = {1,2,3};

        std::cout << "a = " << a << std::endl;
        std::cout << "||a|| = " << Norm(a) << std::endl;
        std::cout << "||a||_{2,2} = " << EntrywiseNorm(a) << std::endl;
        std::cout << "||a||_{max} = " << MaxNorm(a) << std::endl;
        std::cout << "||a||_{inf} = " << InfinityNorm(a) << std::endl;
        std::cout << "v = " << v << std::endl;
        std::cout << "||v|| = " << Norm(v) << std::endl;
        std::cout << "a.*a = " << EntrywiseProduct(a,a) << std::endl;
        std::cout << "a./a = " << EntrywiseDivision(a,a) << std::endl;
        std::cout << "Abs(a) = " << Abs(a) << std::endl;
        std::cout << "Sign(a) = " << Sign(a) << std::endl;
        std::cout << "exp(a) = " << Exp(a) << std::endl;
        std::cout << "pow(a,13) = " << Pow(a,13) << std::endl;

        Matrix3d one = One<double, 3, 3>();
        std::cout << std::endl << "Trig. Functions: " << std::endl;
        std::cout << "||asin(sin(one))||_{inf} = " << InfinityNorm(ASin(Sin(one))) << std::endl;
        std::cout << "||acos(cos(one))||_{inf} = " << InfinityNorm(ACos(Cos(one))) << std::endl;
        std::cout << "||atan(tan(one))||_{inf} = " << InfinityNorm(ATan(Tan(one))) << std::endl;
        std::cout << "||acsc(csc(one))||_{inf} = " << InfinityNorm(ACsc(Csc(one))) << std::endl;
        std::cout << "||asec(sec(one))||_{inf} = " << InfinityNorm(ASec(Sec(one))) << std::endl;
        std::cout << "||acot(cot(one))||_{inf} = " << InfinityNorm(ACot(Cot(one))) << std::endl;

        std::cout << std::endl << "Hyperbolic Trig. Functions: " << std::endl;
        std::cout << "||asinh(sinh(one))||_{inf} = " << InfinityNorm(ASinh(Sinh(one))) << std::endl;
        std::cout << "||acosh(cosh(one))||_{inf} = " << InfinityNorm(ACosh(Cosh(one))) << std::endl;
        std::cout << "||atanh(tanh(one))||_{inf} = " << InfinityNorm(ATanh(Tanh(one))) << std::endl;
        std::cout << "||acsch(csch(one))||_{inf} = " << InfinityNorm(ACsch(Csch(one))) << std::endl;
        std::cout << "||asech(sech(one))||_{inf} = " << InfinityNorm(ASech(Sech(one))) << std::endl;
        std::cout << "||acoth(coth(one))||_{inf} = " << InfinityNorm(ACoth(Coth(one))) << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
