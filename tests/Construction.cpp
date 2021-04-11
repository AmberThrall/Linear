#include "../src/Construction.h"
#include <iostream>
using namespace Linear;

int main() {
    //MatrixXf eye0 = Identity<float,0>();
    Matrix3f eye1 = Identity<float,3>();
    std::cout << "eye1 = " << eye1 << std::endl;
    MatrixXf eye2 = Identity<float>(2);
    std::cout << "eye2 = " << eye2 << std::endl;

    Matrix3f diag1 = Diag<float>(1,2,3);
    std::cout << "diag1 = " << diag1 << std::endl;
    MatrixXf diag2 = Diag<float>(std::vector<Complex<float>>({1,2,3}));
    std::cout << "diag2 = " << diag2 << std::endl;
    Vector3f v = {1,2,3};
    Matrix3f diag3 = Diag(v);
    std::cout << "diag3 = " << diag3 << std::endl;
    RowVector3f row_v = {1,2,3};
    Matrix3f diag4 = Diag(row_v);
    std::cout << "diag4 = " << diag4 << std::endl;
    SquareMatrix<float,1> diag5 = Diag<float>(1);
    std::cout << "diag5 = " << diag5 << std::endl;
}
