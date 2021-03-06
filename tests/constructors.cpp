#include "../src/Matrix.h"
#include <iostream>
using namespace Linear;

int main() {
    Matrix3f matrix_static_1(1);
    matrix_static_1(1,1) = 2;
    std::cout << "matrix_static_1 = " << matrix_static_1 << std::endl;
    std::cout << "size(matrix_static_1) = " << matrix_static_1.Size() << std::endl;

    Matrix2f matrix_static_2(Complexf(1,1));
    std::cout << "matrix_static_2 = " << matrix_static_2 << std::endl;
    std::cout << "size(matrix_static_2) = " << matrix_static_2.Size() << std::endl;

    Matrix3f matrix_static_3 = {
        1, 2, 3, 4, 5, 6, 7, 8, 9
    };
    std::cout << "matrix_static_3 = " << matrix_static_3 << std::endl;
    std::cout << "size(matrix_static_3) = " << matrix_static_3.Size() << std::endl;

    Matrix3f matrix_static_4 = {
        {1,2,3}, {4,5,6}, {7,8,9}
    };
    std::cout << "matrix_static_4 = " << matrix_static_4 << std::endl;
    std::cout << "size(matrix_static_4) = " << matrix_static_4.Size() << std::endl;

    Matrix<float,Dynamic,2> matrix_dynamic_1(2, 1);
    std::cout << "matrix_dynamic_1 = " << matrix_dynamic_1 << std::endl;
    std::cout << "size(matrix_dynamic_1) = " << matrix_dynamic_1.Size() << std::endl;

    Matrix<float,2,Dynamic> matrix_dynamic_2(2, 1);
    std::cout << "matrix_dynamic_2 = " << matrix_dynamic_2 << std::endl;
    std::cout << "size(matrix_dynamic_2) = " << matrix_dynamic_2.Size() << std::endl;

    MatrixXf matrix_dynamic_3(2, 3, 1);
    std::cout << "matrix_dynamic_3 = " << matrix_dynamic_3 << std::endl;
    std::cout << "size(matrix_dynamic_3) = " << matrix_dynamic_3.Size() << std::endl;

    Matrix<float,Dynamic,3> matrix_dynamic_4 = {
        1, 2, 3, 4, 5, 6, 7, 8, 9
    };
    std::cout << "matrix_dynamic_4 = " << matrix_dynamic_4 << std::endl;
    std::cout << "size(matrix_dynamic_4) = " << matrix_dynamic_4.Size() << std::endl;

    Matrix<float,3,Dynamic> matrix_dynamic_5 = {
        1, 2, 3, 4, 5, 6, 7, 8, 9
    };
    std::cout << "matrix_dynamic_5 = " << matrix_dynamic_5 << std::endl;
    std::cout << "size(matrix_dynamic_5) = " << matrix_dynamic_5.Size() << std::endl;

    MatrixXf matrix_dynamic_6 = {
        1, 2, 3, 4, 5, 6, 7, 8, 9
    };
    std::cout << "matrix_dynamic_6 = " << matrix_dynamic_6 << std::endl;
    std::cout << "size(matrix_dynamic_6) = " << matrix_dynamic_6.Size() << std::endl;

    MatrixXf matrix_dynamic_7 = {
        {1,2,3}, {4,5,6}, {7,8,9}
    };
    std::cout << "matrix_dynamic_7 = " << matrix_dynamic_7 << std::endl;
    std::cout << "size(matrix_dynamic_7) = " << matrix_dynamic_7.Size() << std::endl;

    Matrix3f matrix_copy_1(matrix_static_3);
    std::cout << "matrix_copy_1 = " << matrix_copy_1 << std::endl;
    std::cout << "size(matrix_copy_1) = " << matrix_copy_1.Size() << std::endl;

    Matrix3f matrix_copy_2(matrix_dynamic_7);
    std::cout << "matrix_copy_2 = " << matrix_copy_2 << std::endl;
    std::cout << "size(matrix_copy_2) = " << matrix_copy_2.Size() << std::endl;

    MatrixXf matrix_copy_3(matrix_static_3);
    std::cout << "matrix_copy_3 = " << matrix_copy_3 << std::endl;
    std::cout << "size(matrix_copy_3) = " << matrix_copy_3.Size() << std::endl;

    Matrix3d matrix_copy_4(matrix_static_3);
    std::cout << "matrix_copy_4 = " << matrix_copy_4 << std::endl;
    std::cout << "size(matrix_copy_4) = " << matrix_copy_4.Size() << std::endl;

    MatrixXd matrix_empty_1;
    std::cout << "matrix_empty_1 = " << matrix_empty_1 << std::endl;
    std::cout << "size(matrix_empty_1) = " << matrix_empty_1.Size() << std::endl;

    Matrix<double,Dynamic,2> matrix_empty_2;
    std::cout << "matrix_empty_2 = " << matrix_empty_2 << std::endl;
    std::cout << "size(matrix_empty_2) = " << matrix_empty_2.Size() << std::endl;

    Matrix<double,2,Dynamic> matrix_empty_3;
    std::cout << "matrix_empty_3 = " << matrix_empty_3 << std::endl;
    std::cout << "size(matrix_empty_3) = " << matrix_empty_3.Size() << std::endl;
}
