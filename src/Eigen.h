#pragma once
#include <vector>
#include <stack>
#include "Matrix.h"
#include "Vector.h"
#include "Basics.h"
#include "Construction.h"
#include "Decomp.h"
#include "Global.h"
#include "Types.h"
#include "Misc.h"

namespace Linear {
    /**
     * Struct for eigenpairs.
     * This struct contains an eigenvalue and it's corresponding eigenvector of a matrix A.
     * @param T Type to store vector entries as.
     * @param N Size of eigenvector.
    */
    template <typename T, size_t N>
    struct Eigenpair {
        Complex<T> value; /*!< Eigenvalue */
        Vector<T,N> vector; /*!< Eigenvector corresponding to value */
    };

    /**
     * Performs power iteration on A.
     * Power iteration, defined by the sequence \f$b_{k+1}=\frac{Ab_k}{\|Ab_k\|}\f$ is a simple algorithm that computes the greatest
     * (in magnitude) eigenvalue of A along with it's corresponding eigenvector. Convergence is often slow.
     *
     * More information: https://en.wikipedia.org/wiki/Power_iteration
     * @param A MxN matrix
     * @param b0 Initial vector of length P to start the sequence.
     * @param max_iterations Maximum number of iterations to perform
     * @return Pair \f$(\lambda,v)\f$ such that \f$Av\approx \lambda v\f$
     */
    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P>
    typename std::enable_if<(P==N||P==Dynamic||N==Dynamic), Eigenpair<T,P>>::type
    PowerIteration(const Matrix<T,M,N,Flags>& A, Vector<T,P> b0, unsigned int max_iterations) {
        if (!IsSquare(A))
            throw "Eigenvalues are only defined for square matrices.";
        if (A.NumRows() != b0.Length())
            throw "Size mismatch in PowerIteration().";
        if (Norm(b0) < T(Tol))
            throw "Cannot use the zero vector as the initial vector in PowerIteration().";
        if (A.NumEntries() == 0)
            throw "Cannot perform power iteration on a Mx0 or 0xN matrix.";

        Complex<T> lambda;
        Vector<T,N> q = b0/Norm(b0);
        Vector<T,N> b = A*q;
        RowVector<T,N> q2 = ConjugateTranspose(q);
        RowVector<T,N> b2 = ConjugateTranspose(b);
        T res = T(Tol+1);
        for (unsigned int i = 0; res >= T(Tol) && i < max_iterations; ++i) {
            q = Normalize(b);
            b = A*q;
            lambda = Dot(q,b);

            b2 = q2*A;
            q2 = Normalize(b2);
            T costheta = Abs(Dot(q2,q));
            if (costheta < T(Tol)) // Eigenvalues multiplicity >1 ?
                break;
            res = Norm(b-lambda*q)/costheta;
        }

        Eigenpair<T,N> pair;
        pair.value = lambda;
        pair.vector = b;
        return pair;
    }
    /**
     * Performs inverse iteration on A with guess mu.
     * Inverse iteration, defined by the sequence \f$b_{k+1}=\frac{(A-\mu I)^{-1}b_k}{\|(A-\mu I)^{-1}b_k\|}\f$ is a simple algorithm that computes the
     * eigenvalue of A closest to mu along with it's corresponding eigenvector. Convergence rate depends on how close mu is.
     *
     * More information: https://en.wikipedia.org/wiki/Inverse_iteration
     * @param A MxN matrix
     * @param b0 Initial vector of length P to start the sequence
     * @param mu Original guess at an eigenvalue
     * @param max_iterations Maximum number of iterations to perform
     * @return Vector v such that \f$Av\approx\lambda v\f$ where \f$\lambda\f$ is the eigenvalue of A closest to mu
     */
    template <typename T, size_t M, size_t N, unsigned int Flags, size_t P>
    typename std::enable_if<(P==N||P==Dynamic||N==Dynamic), Vector<T,P>>::type
    InverseIteration(const Matrix<T,M,N,Flags>& A, Vector<T,P> b0, Complex<T> mu, unsigned int max_iterations) {
        if (!IsSquare(A))
            throw "Eigenvectors are only defined for square matrices.";
        if (A.NumRows() != b0.Length())
            throw "Size mismatch in PowerIteration().";
        if (Norm(b0) < T(Tol))
            throw "Cannot use the zero vector as the initial vector in InverseIteration().";
        if (A.NumEntries() == 0)
            throw "Cannot perform inverse iteration on a Mx0 or 0xN matrix.";

        SquareMatrix<T,N,Flags> eye = Identity<T,Flags>(A.NumRows());
        SquareMatrix<T,N,Flags> B = Inverse(A-mu*eye);
        Vector<T,N> b = b0;
        for (unsigned int i = 0; i < max_iterations; ++i) {
            b = Normalize(B*b);
        }
        return b;
    }

    /**
     * Performs a combination of deflation and power iteration to get every eigenpair of A
     * Since power iteration returns the largest eigenvalue, if we want to compute the eigenvalues we need to use deflation.
     * The general idea is that once an eigenvalue is calculated, we reduce A to a slightly smaller matrix which has the same eigenvalues
     * except removing the one we just calculated. Repeating this process eventually computes all eigenvalues.
     * @param A MxN matrix
     * @param max_iterations Maximum number of iterations to perform for each call of PowerIteration
     * @return Column vector of eigenvalues.
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Vector<T,N> WielandtDeflationAlgorithm(const Matrix<T,M,N,Flags>& A, unsigned int max_iterations = 100) {
        if (!IsSquare(A))
            throw "Eigenvalues are only defined for square matrices.";
        if (A.NumEntries() == 0) {
            Vector<T,N> eigenvalues(T(0));
            return eigenvalues;
        }

        // Final step: 1x1 matrix.
        if (A.NumRows() == 1) {
            return A;
        }

        // Step 1: Computer the dominant eigenvalue.
        Vector<T,N> eigenvalues(A.NumRows(),1);
        std::vector<Eigenpair<T,N>> eigenpairs;
        Vector<T,N> b0 = Random<T>(A.NumRows(),1,Complex<T>(T(-1),T(-1)), Complex<T>(T(1),T(1)));
        Eigenpair<T,N> pair = PowerIteration(A, b0, max_iterations);
        eigenvalues[A.NumRows()-1] = pair.value;
        Vector<T,N> x = pair.vector;

        // Step 2: Find the p such x1[p] is maximal.
        size_t p = 0;
        T max = T(0);
        for (size_t i = 0; i < A.NumRows(); ++i) {
            if (max < Abs(x[i])) {
                p = i;
                max = Abs(x[i]);
            }
        }
        if (max > T(Tol))
            x = x / x[p];

        // Step 3/4: Compute Ap and remove row p and column p
        RowVector<T,N> ap = A.GetRow(p);
        SquareMatrix<T,Dynamic,Flags> Ap = RemoveRowAndColumn(A - x*ap, p, p);

        // Step 5: Repeat.
        Vector<T,(N==Dynamic?Dynamic:N-1)> res = WielandtDeflationAlgorithm(Ap, max_iterations);
        for (size_t i = 0; i < res.Length(); ++i)
            eigenvalues[i] = res[i];
        return eigenvalues;
    }

    /**
     * Calculates all eigenvalues of A
     *
     * This function handles various cases. First, if the matrix is either upper or lower triangular, then the eigenvalues are simply across the
     * diagonal. Allowing us to pull them directly out of the matrix and use either the Nullspace or InverseIteration to calculate their corresponding
     * eigenvectors. Second, if the matrix is 2x2, then a easy direct formula exists to calculate the eigenvalues. Namely
     * \f$(tr(A)\pm\sqrt{tr(A)^2-4det(A)})/2\f$. Third, if neither of the first two cases handles the matrix, we attempt to calculate the
     * Schur decomposition of A. If that fails, we resort to calling WielandtDeflationAlgorithm.
     * @param A MxN matrix
     * @return Column vector of eigenvalues.
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Vector<T,N> Eigenvalues(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            throw "Eigenvalues are only defined for square matrices.";
        if (A.NumEntries() == 0) {
            Vector<T,N> eigenvalues(T(0));
            return eigenvalues;
        }

        Vector<T,N> eigenvalues(A.NumRows(), T(0));
        if (IsTriangular(A)) {
            for (size_t i = 0; i < A.NumRows(); ++i)
                eigenvalues[i] = A(i,i);
            return eigenvalues;
        }

        if (A.NumRows() == 2) {
            Complex<T> trA = Trace(A), detA = Determinant(A);
            eigenvalues[0] = (trA+Sqrt(trA*trA-4*detA))/2;
            eigenvalues[1] = (trA-Sqrt(trA*trA-4*detA))/2;
            return eigenvalues;
        }

        Schur<T,N,Flags> schur(A);
        bool schurSucceeded = true;
        for (size_t i = 0; i < A.NumRows(); ++i) {
            if (i == A.NumRows()-1 || Abs(schur.U(i+1,i)) < T(Tol)) { // 1x1 block.
                // Check it's all zeros below
                for (size_t r=i+2; r < A.NumRows(); ++r) {
                    if (Abs(schur.U(r,i)) >= T(Tol)) {
                        schurSucceeded = false;
                        break;
                    }
                }
                if (!schurSucceeded)
                    break;
                eigenvalues[i] = schur.U(i,i);
            }
            else if (i < A.NumRows()-1 && Abs(schur.U(i+1,i)) >= T(Tol)) { // 2x2 block.
                // Check it's all zeros below
                for (size_t c=i; c <= i+1; ++c) {
                    for (size_t r=i+2; r < A.NumRows(); ++r) {
                        if (Abs(schur.U(r,c)) >= T(Tol)) {
                            schurSucceeded = false;
                            break;
                        }
                    }
                }
                if (!schurSucceeded)
                    break;

                SquareMatrix<T,2,Flags> block = SubMatrix<2,2>(schur.U, i, i);
                Vector<T,2> eig = Eigenvalues(block);
                eigenvalues[i] = eig[0];
                eigenvalues[i+1] = eig[1];
                i += 1;
            }
            else { // Error?
                schurSucceeded = false;
                break;
            }
        }
        if (schurSucceeded) {
            return eigenvalues;
        }

        // Generic case via Deflation and power iteration.
        return WielandtDeflationAlgorithm(A, 100);
    }

    /**
     * Calculates all eigenpairs of A
     *
     * This function handles various cases. First, if the matrix is either upper or lower triangular, then the eigenvalues are simply across the
     * diagonal. Allowing us to pull them directly out of the matrix and use either the Nullspace or InverseIteration to calculate their corresponding
     * eigenvectors. Second, if the matrix is 2x2, then a easy direct formula exists to calculate the eigenvalues. Namely
     * \f$(tr(A)\pm\sqrt{tr(A)^2-4det(A)})/2\f$. Third, if neither of the first two cases handles the matrix, we attempt to calculate the
     * Schur decomposition of A. If that fails, we resort to calling WielandtDeflationAlgorithm.
     * @param A MxN matrix
     * @return List of Eigenpairs
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    std::vector<Eigenpair<T,N>> Eigen(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            throw "Eigenvalues are only defined for square matrices.";
        if (A.NumEntries() == 0) {
            std::vector<Eigenpair<T,N>> eigenpairs;
            return eigenpairs;
        }

        Vector<T,N> eigenvalues = Eigenvalues(A);

        bool isAReal = IsReal(A);

        SquareMatrix<T,Dynamic,Flags> eye = Identity<T>(A.NumColumns());
        std::vector<Eigenpair<T,N>> eigenpairs;
        std::vector<size_t> repeats;
        for (size_t i = 0; i < eigenvalues.Length(); ++i) {
            // Check if we need to skip this one.
            bool skip = false;
            for (size_t j = 0; j < repeats.size(); ++j) {
                if (repeats[j] == i) {
                    skip = true;
                    break;
                }
            }
            if (skip)
                continue;

            // Find the multiplicity of lambda.
            size_t multiplicity = 1;
            for (size_t j = i+1; j < eigenvalues.Length(); ++j) {
                if (eigenvalues[j] == eigenvalues[i]) {
                    multiplicity += 1;
                    repeats.push_back(j);
                }
            }

            // Get the eigenvector.
            Eigenpair<T,N> pair;
            pair.value = eigenvalues[i];
            std::vector<Vector<T,N>> basis = NullSpace(A - eigenvalues[i]*eye);
            if (basis.size() >= 1) {
                for (size_t j = 0; j < multiplicity; ++j) {
                    if (j < basis.size())
                        pair.vector = Normalize(basis[j]);
                    eigenpairs.push_back(pair);
                }
            }
            else {
                if (Determinant(A-eigenvalues[i]*eye) != T(0)) {
                    Vector<T,N> b0;
                    if (!isAReal || !IsReal(eigenvalues[i]))
                        b0 = Random<T>(A.NumRows(),1,Complex<T>(-1,-1),Complex<T>(1,1));
                    else
                        b0 = Random<T>(A.NumRows(),1,T(-1),T(1));
                    pair.vector = Normalize(InverseIteration(A, b0, eigenvalues[i], 25));
                }
                else
                    pair.vector = Zero<T>(A.NumColumns(),1);

                for (size_t j = 0; j < multiplicity; ++j)
                    eigenpairs.push_back(pair);
            }
        }

        return eigenpairs;
    }

    template <typename T, size_t N>
    void VietaFormulaHelper(size_t data[], size_t start, size_t end, size_t index, size_t k, Complex<T> & sum, const Vector<T,N>& roots) {
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

    /**
     * Computes the characteristic polynomial.
     *
     * This returns the coefficients as a vector c of length N+1 such that \f$\det(A-tI)=c_0+c_1t+\dots+c_{N-1}t^{N-1}+c_Nt^N\f$.
     * @param A MxN matrix
     * @return Coefficients of characteristic polynomial
     */
    template <typename T, size_t M, size_t N, unsigned int Flags>
    Vector<T,(N==Dynamic?Dynamic:N+1)> CharPoly(const Matrix<T,M,N,Flags>& A) {
        if (!IsSquare(A))
            throw "Characteristic polynomial is only defined for square matrices.";
        if (A.NumEntries() == 0)
            throw "Cannot compute the characteristic polynomial of a 0x0, Mx0 or 0xN matrix.";

        // Check if it's the companion matrix.
        if (IsCompanion(A)) {
            Vector<T,(N==Dynamic?Dynamic:N+1)> c(A.NumRows()+1, T(0));
            for (size_t i = 0; i < A.NumRows(); ++i) {
                c[i] = -A(i,A.NumColumns()-1);
            }
            c[A.NumRows()] = 1;
            return c;
        }

        Vector<T,N> eigenvalues = Eigenvalues(A);
        Vector<T,(N==Dynamic?Dynamic:N+1)> c(A.NumRows()+1, T(0));
        for (size_t k = 1; k <= A.NumRows(); ++k) {
            size_t data[k];
            VietaFormulaHelper(data, 0, A.NumRows()-1, 0, k, c[k-1], eigenvalues);
            if (k % 2 == 1)
                c[k-1] *= -1;
        }
        c[A.NumRows()] = T(1);

        return c;
    }
}
