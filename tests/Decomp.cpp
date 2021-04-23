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

        Matrix3f B = {
            {1,2,3},{4,5,6},{7,8,9}
        };
        std::tuple<Matrix3f,Matrix3f,Matrix3f> lup = LUP(B);
        std::cout << "B = " << B << std::endl;
        std::cout << "L = " << std::get<0>(lup) << std::endl;
        std::cout << "U = " << std::get<1>(lup) << std::endl;
        std::cout << "P = " << std::get<2>(lup) << std::endl;
        std::cout << "PB = " << std::get<2>(lup)*B << std::endl;
        std::cout << "LU = " << std::get<0>(lup)*std::get<1>(lup) << std::endl;

        Matrix3f c = {
            {4,12,-16}, {12,37,-43}, {-16,-43,98}
        };

        Matrix3f l = LL(c);
        std::cout << "c = " << c << std::endl;
        std::cout << "l = " << l << std::endl;
        std::cout << "ll* = " << l*ConjugateTranspose(l) << std::endl;

        std::pair<Matrix3f,Matrix3f> ld = LDL(c);
        std::cout << "c = " << c << std::endl;
        std::cout << "l = " << ld.first << std::endl;
        std::cout << "d = " << ld.second << std::endl;
        std::cout << "ldl* = " << ld.first*ld.second*ConjugateTranspose(ld.first) << std::endl;

        Matrix2f d = {
            {1,0}, {1,3}
        };
        std::tuple<Matrix2f,Matrix2f,Matrix2f> qd = Eigendecomposition(d);
        std::cout << "d = " << d << std::endl;
        std::cout << "q = " << std::get<0>(qd) << std::endl;
        std::cout << "d = " << std::get<1>(qd) << std::endl;
        std::cout << "q^{-1} = " << std::get<2>(qd) << std::endl;
        std::cout << "qdq^{-1} = " << std::get<0>(qd)*std::get<1>(qd)*std::get<2>(qd) << std::endl;

        Matrix<float,4,5> e = {
            {1,0,0,0,2}, {0,0,3,0,0}, {0,0,0,0,0}, {0,2,0,0,0}
        };
        std::tuple<Matrix4f,Matrix<float,4,5>,Matrix<float,5,5>> svd = SVD(e);
        std::cout << "e = " << e << std::endl;
        std::cout << "u = " << std::get<0>(svd) << std::endl;
        std::cout << "s = " << std::get<1>(svd) << std::endl;
        std::cout << "v* = " << ConjugateTranspose(std::get<2>(svd)) << std::endl;
        std::cout << "usv* = " << std::get<0>(svd)*std::get<1>(svd)*ConjugateTranspose(std::get<2>(svd)) << std::endl;
        std::cout << "uu* = " << std::get<0>(svd)*ConjugateTranspose(std::get<0>(svd)) << std::endl;
        std::cout << "vv* = " << std::get<2>(svd)*ConjugateTranspose(std::get<2>(svd)) << std::endl;

        Matrix4d F  = {
            {1,2,3,4}, {5,6,7,8}, {9,10,11,12}, {13,14,15,16}
        };
        std::pair<Matrix4d,Matrix4d> qh = Hessenberg(F);
        std::cout << "F = " << F << std::endl;
        std::cout << "Q = " << qh.first << std::endl;
        std::cout << "H = " << qh.second << std::endl;
        std::cout << "QHQ* = " << qh.first*qh.second*ConjugateTranspose(qh.first) << std::endl;
        std::cout << "QQ* = " << qh.first*ConjugateTranspose(qh.first);

        std::pair<Matrix4d,Matrix4d> schur = Schur(F);
        std::cout << "F = " << F << std::endl;
        std::cout << "Q = " << schur.first << std::endl;
        std::cout << "U = " << schur.second << std::endl;
        std::cout << "QUQ* = " << schur.second*schur.first*ConjugateTranspose(schur.second) << std::endl;
        std::cout << "UU* = " << schur.second*ConjugateTranspose(schur.second);
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
