#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
        Matrix3f a = {
            {12,-51,4}, {6,167,-68}, {-4,24,-41}
        };

        std::pair<Matrix3f,Matrix3f> qr = QR(a);
        std::cout << "a = " << a << std::endl;
        std::cout << "q = " << qr.first << std::endl;
        std::cout << "r = " << qr.second << std::endl;
        std::cout << "qr = " << qr.first*qr.second << std::endl;
        std::cout << "q^Tq = " << Transpose(qr.first)*qr.first << std::endl;

        Matrix3f b = {
            {1,2,3},{4,5,6},{7,8,9}
        };
        std::tuple<Matrix3f,Matrix3f,Matrix3f> lup = LUP(b);
        std::cout << "b = " << b << std::endl;
        std::cout << "l = " << std::get<0>(lup) << std::endl;
        std::cout << "u = " << std::get<1>(lup) << std::endl;
        std::cout << "p = " << std::get<2>(lup) << std::endl;
        std::cout << "pb = " << std::get<2>(lup)*b << std::endl;
        std::cout << "lu = " << std::get<0>(lup)*std::get<1>(lup) << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
