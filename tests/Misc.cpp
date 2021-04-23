#include "../src/Linear.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
        Matrix<double,4,6> A = {
            {1,0,-3,0,2,-8}, {0,1,5,0,-1,4}, {0,0,0,1,7,-9}, {0,0,0,0,0,0}
        };
        std::cout << "A = " << A << std::endl;
        std::vector<Vector<double,6>> null_basis = NullSpace(A);
        for (unsigned int i = 0; i < null_basis.size(); ++i) {
            std::cout << "null_basis["<<i<<"] = " << Transpose(null_basis[i]) << std::endl;
            std::cout << "A*null_basis["<<i<<"] = " << Transpose(A*null_basis[i]) << std::endl;
        }

        Matrix4d B = {
            {1,3,1,4}, {2,7,3,9}, {1,5,3,1}, {1,2,0,8}
        };
        std::cout << "B = " << B << std::endl;
        std::vector<Vector4d> cs_basis = ColumnSpace(B);
        for (unsigned int i = 0; i < null_basis.size(); ++i)
            std::cout << "cs_basis["<<i<<"] = " << Transpose(cs_basis[i]) << std::endl;

        unsigned int rank = Rank(B);
        unsigned int nullity = Nullity(B);
        std::cout << "rank(B)+nullity(B) = " << rank << " + " << nullity << " = " << rank+nullity << std::endl;

        Matrix<double,3,4> C = {
            {1,2,2,2}, {2,4,6,8}, {3,6,8,10}
        };
        Vector3d b = {1,5,6};
        Vector4d x = Solve(C, b);
        std::cout << "C = " << C << std::endl;
        std::cout << "b = " << Transpose(b) << std::endl;
        std::cout << "x = " << Transpose(x) << std::endl;
        std::cout << "Cx = " << Transpose(C*x) << std::endl;

        Matrix<double,3,2> D = {
            {1,2}, {2,4}, {3,6}
        };
        Vector2d x2 = {1, 2};
        b = D*x2;
        x2 = Solve(D, b);
        std::cout << "D = " << D << std::endl;
        std::cout << "b = " << Transpose(b) << std::endl;
        std::cout << "x2 = " << Transpose(x2) << std::endl;
        std::cout << "Dx = " << Transpose(D*x2) << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
