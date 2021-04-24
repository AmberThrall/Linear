#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
        Matrix3f A = {
            {12,-51,4}, {6,167,-68}, {-4,24,-41}
        };

        QR<float,3,3> qr(A);
        std::cout << "A = " << A << std::endl;
        std::cout << "Q = " << qr.Q << std::endl;
        std::cout << "R = " << qr.R << std::endl;
        std::cout << "QR = " << qr.Q*qr.R << std::endl;
        std::cout << "Q*Q = " << ConjugateTranspose(qr.Q)*qr.Q << std::endl;
        std::cout << "=======================" << std::endl;

        Matrix3f B = {
            {1,2,3},{4,5,6},{7,8,9}
        };
        LUP<float,3> lup(B);
        std::cout << "B = " << B << std::endl;
        std::cout << "L = " << lup.L << std::endl;
        std::cout << "U = " << lup.U << std::endl;
        std::cout << "P = " << lup.P << std::endl;
        std::cout << "PB = " << lup.P*B << std::endl;
        std::cout << "LU = " << lup.L*lup.U << std::endl;
        std::cout << "=======================" << std::endl;

        Matrix3d C = {
            {4,12,-16}, {12,37,-43}, {-16,-43,98}
        };

        Cholesky<double,3> cholesky(C);
        std::cout << "C = " << C << std::endl;
        std::cout << "L = " << cholesky.L << std::endl;
        std::cout << "D = " << cholesky.D << std::endl;
        std::cout << "LDL* = " << cholesky.L*cholesky.D*cholesky.Lh << std::endl;
        std::cout << "=======================" << std::endl;

        Matrix2d E = {
            {1,0}, {1,3}
        };
        Eigendecomposition<double,2> eigen(E);
        std::cout << "E = " << E << std::endl;
        std::cout << "Q = " << eigen.Q << std::endl;
        std::cout << "D = " << eigen.D << std::endl;
        std::cout << "Q^{-1} = " << eigen.Qinv << std::endl;
        std::cout << "QDQ^{-1} = " << eigen.Q*eigen.D*eigen.Qinv << std::endl;
        std::cout << "=======================" << std::endl;

        Matrix4d F = {
            {1,2,3,4}, {5,6,7,8}, {9,10,11,12}, {13,14,15,16}
        };
        Hessenberg<double,4> hess(F);
        std::cout << "F = " << F << std::endl;
        std::cout << "Q = " << hess.Q << std::endl;
        std::cout << "H = " << hess.H << std::endl;
        std::cout << "QHQ* = " << hess.Q*hess.H*hess.Qh << std::endl;
        std::cout << "QQ* = " << hess.Q*hess.Qh << std::endl;
        std::cout << "=======================" << std::endl;

        Schur<double,4> schur(F);
        std::cout << "F = " << F << std::endl;
        std::cout << "Q = " << schur.Q << std::endl;
        std::cout << "U = " << schur.U << std::endl;
        std::cout << "QUQ* = " << schur.Q*schur.U*schur.Qh << std::endl;
        std::cout << "QQ* = " << schur.Q*schur.Qh;
        std::cout << "=======================" << std::endl;

        Matrix<double,4,5> G = {
            {1,0,0,0,2}, {0,0,3,0,0}, {0,0,0,0,0}, {0,2,0,0,0}
        };
        SVD<double,4,5> svd(G);
        std::cout << "G = " << G << std::endl;
        std::cout << "U = " << svd.U << std::endl;
        std::cout << "S = " << svd.S << std::endl;
        std::cout << "V* = " << svd.Vh << std::endl;
        std::cout << "USV* = " << svd.U*svd.S*svd.Vh << std::endl;
        std::cout << "UU* = " << svd.U*ConjugateTranspose(svd.U) << std::endl;
        std::cout << "V*V = " << svd.Vh*ConjugateTranspose(svd.Vh) << std::endl;
        std::cout << "=======================" << std::endl;

        SVD<double,4,5,0,THIN_SVD> svdThin(G);
        std::cout << "G = " << G << std::endl;
        std::cout << "U = " << svdThin.U << std::endl;
        std::cout << "S = " << svdThin.S << std::endl;
        std::cout << "V* = " << svdThin.Vh << std::endl;
        std::cout << "USV = " << svdThin.U*svdThin.S*svdThin.Vh << std::endl;
        std::cout << "=======================" << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
