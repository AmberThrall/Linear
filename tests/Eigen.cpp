#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
        Matrix3f a = {
            {1,0,0}, {2,3,0}, {4,5,6}
        };

        std::cout << "a = " << a << std::endl;
        std::vector<std::pair<Complexf,Vector3f>> eigenpairs_a = Eigen(a);
        for (unsigned int i = 0; i < eigenpairs_a.size(); ++i) {
            std::cout << "(" << eigenpairs_a[i].first << "," << Transpose(eigenpairs_a[i].second) << ")" << std::endl;
        }

        Matrix3f b = { {-4,14,0}, {-5,13,0}, {-1,0,2} };

        std::cout << "b = " << b << std::endl;
        std::vector<std::pair<Complexf,Vector3f>> eigenpairs_b = Eigen(b);
        for (unsigned int i = 0; i < eigenpairs_b.size(); ++i) {
            std::cout << "(" << eigenpairs_b[i].first << "," << Transpose(eigenpairs_b[i].second) << ")" << std::endl;
        }

        Matrix4f c = { {4,1,-2,2}, {1,2,0,1}, {-2,0,3,-2}, {2,1,-2,-1} };
        std::cout << "c = " << c << std::endl;
        std::cout << "Tridiagonalization(c) = " << Tridiagonalization(c) << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
