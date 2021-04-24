#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
        Matrix3d A = {
            {1,0,0}, {2,3,0}, {4,5,6}
        };

        std::cout << "A = " << A << std::endl;
        std::vector<Eigenpair<double,3>> eigenpairs_a = Eigen(A);
        for (size_t i = 0; i < eigenpairs_a.size(); ++i) {
            std::cout << "(" << eigenpairs_a[i].value << "," << Transpose(eigenpairs_a[i].vector) << ")" << std::endl;
        }
        for (size_t i = 0; i < eigenpairs_a.size(); ++i) {
            std::cout << "||Av_i/lambda_i-v_i|| = " << Norm((A*eigenpairs_a[i].vector)/eigenpairs_a[i].value-eigenpairs_a[i].vector) << std::endl;
        }
        std::cout << "==================" << std::endl;

        Matrix3d B = { {-4,14,0}, {-5,13,0}, {-1,0,2} };

        std::cout << "B = " << B << std::endl;
        std::vector<Eigenpair<double,3>> eigenpairs_b = Eigen(B);
        for (size_t i = 0; i < eigenpairs_b.size(); ++i) {
            std::cout << "(" << eigenpairs_b[i].value << "," << Transpose(eigenpairs_b[i].vector) << ")" << std::endl;
        }
        for (size_t i = 0; i < eigenpairs_b.size(); ++i) {
            std::cout << "||Bv_i/lambda_i-v_i|| = " << Norm((B*eigenpairs_b[i].vector)/eigenpairs_b[i].value-eigenpairs_b[i].vector) << std::endl;
        }
        std::cout << "==================" << std::endl;

        Matrix3d C = { {1,0,0}, {1,2,0}, {1,2,2} };

        std::cout << "C = " << C << std::endl;
        std::vector<Eigenpair<double,3>> eigenpairs_c = Eigen(C);
        for (size_t i = 0; i < eigenpairs_c.size(); ++i) {
            std::cout << "(" << eigenpairs_c[i].value << "," << Transpose(eigenpairs_c[i].vector) << ")" << std::endl;
        }
        for (size_t i = 0; i < eigenpairs_c.size(); ++i) {
            std::cout << "||Cv_i/lambda_i-v_i|| = " << Norm((C*eigenpairs_c[i].vector)/eigenpairs_c[i].value-eigenpairs_c[i].vector) << std::endl;
        }
        std::cout << "==================" << std::endl;

        Vector<double,5> coeff = {-1,2,-3,4,1};
        Matrix4d D = Companion(coeff);
        Vector4d eigenvalues_d = Eigenvalues(D);
        std::cout << "coeff = " << Transpose(coeff) << std::endl;
        std::cout << "D = " << D << std::endl;
        std::cout << "eig(D) = ";
        for (size_t i = 0; i < eigenvalues_d.Length(); ++i)
            std::cout << (i > 0 ? ", " : "") << eigenvalues_d[i];
        std::cout << std::endl;
        std::cout << "charpoly(D) = " << Transpose(CharPoly(D)) << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
