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
        for (unsigned int i = 0; i < eigenpairs_a.size(); ++i) {
            std::cout << "(" << eigenpairs_a[i].value << "," << Transpose(eigenpairs_a[i].vector) << ")" << std::endl;
        }

        Matrix3d B = { {-4,14,0}, {-5,13,0}, {-1,0,2} };

        std::cout << "B = " << B << std::endl;
        std::vector<Eigenpair<double,3>> eigenpairs_b = Eigen(B);
        for (unsigned int i = 0; i < eigenpairs_b.size(); ++i) {
            std::cout << "(" << eigenpairs_b[i].value << "," << Transpose(eigenpairs_b[i].vector) << ")" << std::endl;
        }

        Vector<double,5> coeff = {-1,2,-3,4,1};
        Matrix4d C = Companion(coeff);
        std::vector<Eigenpair<double,4>> eigenpairs_c = Eigen(C);
        std::cout << "coeff = " << Transpose(coeff) << std::endl;
        std::cout << "C = " << C << std::endl;
        std::cout << "eig(C) = ";
        for (unsigned int i = 0; i < eigenpairs_c.size(); ++i)
            std::cout << (i > 0 ? ", " : "") << eigenpairs_c[i].value;
        std::cout << std::endl;
        std::cout << "charpoly(C) = " << CharPoly(C) << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
