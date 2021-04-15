#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
        Matrix3d a = {
            {1,0,0}, {2,3,0}, {4,5,6}
        };

        std::cout << "a = " << a << std::endl;
        std::vector<std::pair<Complexd,Vector3d>> eigenpairs_a = Eigen(a);
        for (unsigned int i = 0; i < eigenpairs_a.size(); ++i) {
            std::cout << "(" << eigenpairs_a[i].first << "," << Transpose(eigenpairs_a[i].second) << ")" << std::endl;
        }

        Matrix3d b = { {-4,14,0}, {-5,13,0}, {-1,0,2} };

        std::cout << "b = " << b << std::endl;
        std::vector<std::pair<Complexd,Vector3d>> eigenpairs_b = Eigen(b);
        for (unsigned int i = 0; i < eigenpairs_b.size(); ++i) {
            std::cout << "(" << eigenpairs_b[i].first << "," << Transpose(eigenpairs_b[i].second) << ")" << std::endl;
        }

        Vector<double,5> coeff = {-1,2,-3,4,1};
        Matrix4d c = Companion(coeff);
        std::vector<std::pair<Complexd,Vector4d>> eigenpairs_c = Eigen(c);
        std::cout << "coeff = " << Transpose(coeff) << std::endl;
        std::cout << "c = " << c << std::endl;
        std::cout << "eig(c) = ";
        for (unsigned int i = 0; i < eigenpairs_c.size(); ++i)
            std::cout << (i > 0 ? ", " : "") << eigenpairs_c[i].first;
        std::cout << std::endl;
        std::cout << "charpoly(c) = " << CharPoly(c) << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
