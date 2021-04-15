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
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
