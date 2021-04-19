#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
        Matrix<double,4,6> A = {
            {1,0,-3,0,2,-8}, {0,1,5,0,-1,4}, {0,0,0,1,7,-9}, {0,0,0,0,0,0}
        };
        std::cout << "A = " << A << std::endl;
        std::vector<Vector<double,6>> null_basis = Nullspace(A);
        for (unsigned int i = 0; i < null_basis.size(); ++i) {
            std::cout << "null_basis["<<i<<"] = " << Transpose(null_basis[i]) << std::endl;
            std::cout << "A*null_basis["<<i<<"] = " << Transpose(A*null_basis[i]) << std::endl;
        }
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
