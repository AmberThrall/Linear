#pragma once
#include <vector>
#include <stack>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"
#include "Construction.h"
#include "Global.h"
#include "Types.h"

namespace Linear {
    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)&&(P==N||P==Dynamic||N==Dynamic)), std::pair<Complex<T>,Vector<T,P>>>::type
    PowerIteration(const Matrix<T,M,N,Flags>& A, Vector<T,P> b0, unsigned int max_iterations) {
        if (!IsSquare(A))
            throw "Eigenvalues are only defined for square matrices.";
        if (A.NumRows() != b0.NumRows())
            throw "Size mismatch in PowerIteration().";
        if (Length(b0) < T(Tol))
            throw "Cannot use the zero vector as the initial vector in PowerIteration().";

        Complex<T> lambda;
        Vector<T,N> q = b0/Length(b0);
        Vector<T,N> b = A*q;
        RowVector<T,N> q2 = ConjugateTranspose(q);
        RowVector<T,N> b2 = ConjugateTranspose(b);
        T res = T(Tol+1);
        for (unsigned int i = 0; res >= T(Tol) && i < max_iterations; ++i) {
            T len = Length(b);
            if (len > T(Tol))
                q = b/len;
            else
                return std::make_pair(lambda,b);
            b = A*q;
            lambda = Dot(q,b);

            b2 = q2*A;
            q2 = b2/Length(b2);
            T costheta = Abs(Dot(q2,q));
            if (costheta < T(Tol)) // Eigenvalues multiplicity >1 ?
                break;
            res = Length(b-lambda*q)/costheta;
        }
        return std::make_pair(lambda, b);
    }

    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)&&(P==N||P==Dynamic||N==Dynamic)), Vector<T,P>>::type
    InverseIteration(const Matrix<T,M,N,Flags>& A, Vector<T,P> b0, Complex<T> mu, unsigned int max_iterations) {
        if (!IsSquare(A))
            throw "Eigenvectors are only defined for square matrices.";
        if (A.NumRows() != b0.NumRows())
            throw "Size mismatch in PowerIteration().";
        if (Length(b0) < T(Tol))
            throw "Cannot use the zero vector as the initial vector in InverseIteration().";

        SquareMatrix<T,N,Flags> eye = Identity<T,Flags>(A.NumRows());
        SquareMatrix<T,N,Flags> B = Inverse(A-mu*eye);
        Vector<T,N> b = b0;
        for (unsigned int i = 0; i < max_iterations; ++i) {
            Vector<T,N> bnext = B*b;
            T len = Length(bnext);
            if (len > T(Tol))
                b = bnext/len;
            else
                return b;
        }
        return b;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    std::vector<Vector<T,N>> Nullspace(const Matrix<T,M,N,Flags>& m) {
        std::vector<Vector<T,N>> basis;

        SquareMatrix<T,N,Flags> eye = Identity<T>(m.NumColumns());
        Matrix<T,(M==Dynamic||N==Dynamic?Dynamic:M+N),N,Flags> aug = RowAugmented(m, eye);
        aug = Transpose(RREF(Transpose(aug)));
        SquareMatrix<T,N,Flags> c = SubMatrix(aug, m.NumColumns(), m.NumColumns(), m.NumRows(), 0);
        for (size_t col = 0; col < m.NumColumns(); ++col) {
            bool add = true;
            for (size_t row = 0; row < m.NumRows(); ++row) {
                if (aug(row, col) != T(0)) {
                    add = false;
                    break;
                }
            }
            if (add)
                basis.push_back(c.GetColumn(col));
        }
        return basis;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)), std::vector<std::pair<Complex<T>,Vector<T,N>>>>::type
    WielandtDeflationAlgorithm(const Matrix<T,M,N,Flags>& a, unsigned int num_simulations = 15) {
        if (!IsSquare(a))
            throw "Eigenvalues are only defined for square matrices.";

        // Final step: 1x1 matrix.
        if (a.NumRows() == 1) {
            std::vector<std::pair<Complex<T>,Vector<T,N>>> eigenpairs;
            eigenpairs.push_back(std::make_pair(a(0,0), Vector<T,1>(T(1))));
            return eigenpairs;
        }

        // Step 1: Computer the dominant eigenvalue and vector.
        std::vector<std::pair<Complex<T>,Vector<T,N>>> eigenpairs;
        Vector<T,N> b = Random<T>(a.NumRows(),1,Complex<T>(T(-1),T(-1)), Complex<T>(T(1),T(1)));
        std::pair<Complex<T>,Vector<T,N>> pair = PowerIteration(a, b, num_simulations);
        // Sanity check:
        Complex<T> lambda1 = pair.first;
        Vector<T,N> x1 = pair.second;

        // Step 2: Find the p such x1[p] is maximal.
        size_t p = 0;
        T max = T(0);
        for (size_t i = 0; i < a.NumRows(); ++i) {
            if (max < Abs(x1[i])) {
                p = i;
                max = Abs(x1[i]);
            }
        }
        if (max > T(Tol))
            x1 = x1 / x1[p];

        // Step 3/4: Compute Ap and remove row p and column p
        RowVector<T,N> ap = a.GetRow(p);
        SquareMatrix<T,Dynamic,Flags> Ap = RemoveRowAndColumn(a - x1*ap, p, p);

        // Step 5: Repeat.
        eigenpairs.push_back(std::make_pair(lambda1, x1));
        std::vector<std::pair<Complex<T>,Vector<T,Dynamic>>> res = WielandtDeflationAlgorithm(Ap, num_simulations);
        for (size_t i = 0; i < res.size(); ++i) {
            Vector<T,N> yi(a.NumRows(), T(0));
            for (size_t j = 0; j < yi.Size(); ++j) {
                if (j == p)
                    yi[j] = T(0);
                else if (j > p)
                    yi[j] = res[i].second[j-1];
                else
                    yi[j] = res[i].second[j];
            }
            yi = yi + (Dot(ap, yi)/(res[i].first-lambda1))*x1;
            eigenpairs.push_back(std::make_pair(res[i].first, yi));
        }
        return eigenpairs;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)), std::vector<std::pair<Complex<T>,Vector<T,N>>>>::type
    EigenvectorHelper(const Matrix<T,M,N,Flags>& A, std::vector<Complex<T>> eigenvalues) {
        if (!IsSquare(A))
            throw "Eigenvectors are only defined for square matrices.";

        bool isReal = IsReal(A);
        SquareMatrix<T,Dynamic,Flags> eye = Identity<T>(A.NumColumns());
        std::vector<std::pair<Complex<T>,Vector<T,N>>> eigenpairs;
        for (size_t i = 0; i < eigenvalues.size(); ++i) {
            size_t multiplicity = 1;
            for (size_t j = i+1; j < eigenvalues.size(); ++j) {
                if (eigenvalues[i] == eigenvalues[j]) {
                    multiplicity += 1;
                    eigenvalues.erase(eigenvalues.begin()+j);
                    j -= 1;
                }
            }

            std::vector<Vector<T,N>> basis = Nullspace(A - eigenvalues[i]*eye);
            if (basis.size() != multiplicity) {
                Vector<T,N> b0;
                if (!isReal || !IsReal(eigenvalues[i]))
                    b0 = Random<T>(A.NumRows(), 1, Complex<T>(T(-1),T(-1)), Complex<T>(T(1),T(1)));
                else
                    b0 = Random<T>(A.NumRows(), 1, T(-1), T(1));
                Vector<T,N> eigenvector = InverseIteration(A, b0, eigenvalues[i], 25);
                for (size_t j = 0; j < multiplicity; ++j)
                    eigenpairs.push_back(std::make_pair(eigenvalues[i], eigenvector));
            }
            else {
                for (size_t j = 0; j < multiplicity; ++j)
                    eigenpairs.push_back(std::make_pair(eigenvalues[i], basis[j]));
            }
        }

        return eigenpairs;
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<((M==N||M==Dynamic||N==Dynamic)), std::vector<std::pair<Complex<T>,Vector<T,N>>>>::type
    Eigen(const Matrix<T,M,N,Flags>& A, unsigned int num_simulations = 15) {
        if (!IsSquare(A))
            throw "Eigenvalues are only defined for square matrices.";

        std::vector<std::pair<Complex<T>,Vector<T,N>>> eigenpairs;
        if (IsTriangular(A)) {
            std::vector<Complex<T>> eigenvalues;
            for (size_t i = 0; i < A.NumRows(); ++i) {
                eigenvalues.push_back(A(i,i));
            }
            return EigenvectorHelper(A, eigenvalues);
        }

        if (A.NumRows() == 2) {
            Complex<T> trA = Trace(A), detA = Determinant(A);
            Complex<T> lambda1 = (trA+Sqrt(trA*trA-4*detA))/2;
            Complex<T> lambda2 = (trA-Sqrt(trA*trA-4*detA))/2;
            return EigenvectorHelper(A, std::vector<Complex<T>>({lambda1, lambda2}));
        }

        std::pair<SquareMatrix<T,N,Flags>,SquareMatrix<T,N,Flags>> schur = Schur(A);
        bool schurFailed = false;
        std::vector<Complex<T>> eigenvalues;
        for (size_t i = 0; i < A.NumRows(); ++i) {
            if (i == A.NumRows()-1 || Abs(schur.first(i+1,i)) < T(Tol)) { // 1x1 block.
                // Check it's all zeros below
                for (size_t r=i+2; r < A.NumRows(); ++r) {
                    if (Abs(schur.first(r,i)) >= T(Tol)) {
                        schurFailed = true;
                        break;
                    }
                }
                if (schurFailed)
                    break;
                eigenvalues.push_back(schur.first(i,i));
            }
            else if (i < A.NumRows()-1 && Abs(schur.first(i+1,i)) >= T(Tol)) { // 2x2 block.
                // Check it's all zeros below
                for (size_t c=i; c <= i+1; ++c) {
                    for (size_t r=i+2; r < A.NumRows(); ++r) {
                        if (Abs(schur.first(r,c)) >= T(Tol)) {
                            schurFailed = true;
                            break;
                        }
                    }
                }
                if (schurFailed)
                    break;

                SquareMatrix<T,2,Flags> block = SubMatrix<2,2>(schur.first, i, i);
                std::vector<std::pair<Complex<T>,Vector<T,2>>> eig = Eigen(block);
                eigenvalues.push_back(eig[0].first);
                eigenvalues.push_back(eig[1].first);
                i += 1;
            }
            else { // Error?
                schurFailed = true;
                break;
            }
        }
        if (schurFailed == false && eigenvalues.size() == A.NumRows()) {
            return EigenvectorHelper(A, eigenvalues);
        }

        // Generic case via Deflation and power iteration.
        return WielandtDeflationAlgorithm(A, num_simulations);
    }

    template <typename T>
    void VietaFormulaHelper(size_t data[], size_t start, size_t end, size_t index, size_t k, Complex<T> & sum, const std::vector<Complex<T>>& roots) {
        if (index == k) {
            Complex<T> prod(T(1));
            for (size_t i = 0; i < k; ++i)
                prod *= roots[data[i]];
            sum += prod;
            return;
        }

        for (size_t i = start; i <= end && end - i + 1 >= k - index; ++i) {
            data[index] = i;
            VietaFormulaHelper(data, i+1, end, index+1, k, sum, roots);
        }
    }

    template <typename T, size_t M, size_t N, unsigned int Flags>
    typename std::enable_if<(M==N||M==Dynamic||N==Dynamic), Vector<T,(N==Dynamic?Dynamic:N+1)>>::type CharPoly(const Matrix<T,M,N,Flags>& a) {
        if (!IsSquare(a))
            throw "Characteristic polynomial is only defined for square matrices.";

        std::vector<std::pair<Complex<T>,Vector<T,N>>> eigenpairs = Eigen(a);
        std::vector<Complex<T>> eigenvalues;
        for (size_t i = 0; i < eigenpairs.size(); ++i)
            eigenvalues.push_back(eigenpairs[i].first);
        Vector<T,(N==Dynamic?Dynamic:N+1)> c(a.NumRows()+1, T(0));
        for (size_t k = 1; k <= a.NumRows(); ++k) {
            size_t data[k];
            VietaFormulaHelper(data, 0, a.NumRows()-1, 0, k, c[k-1], eigenvalues);
            if (k % 2 == 1)
                c[k-1] *= -1;
        }
        c[a.NumRows()] = T(1);

        return c;
    }
}
