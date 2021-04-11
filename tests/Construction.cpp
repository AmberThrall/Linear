#include "../src/Construction.h"
#include <iostream>
using namespace Linear;

int main() {
    try {
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
        MatrixXf diag5 = Diag(eye1, diag1, eye2);
        std::cout << "diag5 = " << diag5 << std::endl;

        Matrix2f zero1 = Zero<float,2,2,RowMajor>();
        std::cout << "zero1 = " << zero1 << std::endl;
        MatrixXf zero2 = Zero<float>(2,3);
        std::cout << "zero2 = " << zero2 << std::endl;

        Matrix2f one1 = One<float,2,2>();
        std::cout << "one1 = " << one1 << std::endl;
        MatrixXf one2 = One<float>(2,3);
        std::cout << "one2 = " << one2 << std::endl;

        Matrix2f constant1 = Constant<float,2,2>(Complex<float>::i());
        std::cout << "constant1 = " << constant1 << std::endl;
        MatrixXf constant2 = Constant<float>(2,3,Complex<float>::i());
        std::cout << "constant2 = " << constant2 << std::endl;

        Matrix2f kronecker1 = {1,2,3,4};
        Matrix2f kronecker2 = {0,5,6,7};
        std::cout << "kronecker(kronecker1,kronecker2) = " << Kronecker(kronecker1,kronecker2) << std::endl;
        std::cout << "KroneckerSum(kronecker1,kronecker2) = " << KroneckerSum(kronecker1,kronecker2) << std::endl;

        Matrix2f blocktl1 = Constant<float,2,2>(1);
        Matrix<float,2,3> blocktr1 = Constant<float,2,3>(2);
        Matrix<float,2,3> blockbl1 = Constant<float,2,3>(3);
        Matrix2f blockbr1 = Constant<float,2,2>(4);
        Matrix<float,4,5> block1 = Block(blocktl1, blocktr1, blockbl1, blockbr1);
        std::cout << "block1 = " << block1 << std::endl;
        Matrix2f blocktl2 = blocktl1;
        Vector2f blocktr2 = { 2, 2 };
        RowVector2f blockbl2 = { 3, 3 };
        Complexf blockbr2(4);
        Matrix3f block2 = Block(blocktl2, blocktr2, blockbl2, blockbr2);
        Matrix3f block3 = Block(blockbr2, blockbl2, blocktr2, blocktl2);
        std::cout << "block2 = " << block2 << std::endl;
        std::cout << "block3 = " << block3 << std::endl;

        Matrix3f companion = CompanionMatrix(std::vector<Complexf>({1,2,3}));
        std::cout << "companion = " << companion << std::endl;
    }
    catch (const char* what) {
        std::cerr << "Error: " << what << std::endl;
    }
}
