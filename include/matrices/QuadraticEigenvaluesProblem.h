//
// Created by panagiotis on 2/22/20.
//

#ifndef VOLESTI_QUADRATICEIGENVALUESPROBLEM_H
#define VOLESTI_QUADRATICEIGENVALUESPROBLEM_H


/// Uncomment to use Eigen solver for generalized eigenvalue problem; otherwise Spectra
#define EIGEN_EIGENVALUES_SOLVER

#include "DenseProductMatrix.h"
#include "../../external/Spectra/include/Spectra/GenEigsSolver.h"

/// Solves a quadratic eigenvalue problem
/// \tparam NT Numeric Type
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
template<typename NT, typename MT, typename VT>
class QuadraticEigenvaluesProblem {

};

/// Solves a quadratic eigenvalue problem, A template specialization for dense Eigen matrix and vector
/// \tparam NT
template<typename NT>
class QuadraticEigenvaluesProblem<NT, Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Matrix<NT,Eigen::Dynamic,1> > {
public:
    /// The type for Eigen Matrix
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// Transform the quadratic eigenvalue problem \[At^2 + Bt + c\] to
    /// the generalized eigenvalue problem X+lY.
    /// If the updateOnly flag is false, compute matrices X,Y from scratch;
    /// otherwise update them.
    /// \param[in] A
    /// \param[in] B
    /// \param[in] C
    /// \param[in, out] X
    /// \param[in, out] Y
    /// \param[in, out] updateOnly True if X,Y were previously computed and only B,C changed
    void linearization(const MT &A, const MT &B, const MT &C, MT &X, MT &Y, bool &updateOnly) {
        unsigned int matrixDim = A.rows();

        // check if the matrices X,Y are computed.
        //if yes, update them; otherwise compute them from scratch
        if (!updateOnly) {
            X.resize(2 * matrixDim, 2 * matrixDim);
            Y.resize(2 * matrixDim, 2 * matrixDim);

            Y.block(matrixDim, matrixDim, matrixDim, matrixDim) = -1 * C;
            Y.block(0, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
            Y.block(matrixDim, 0, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
            Y.block(0, 0, matrixDim, matrixDim) = A;

            X.block(0, matrixDim, matrixDim, matrixDim) = C;
            X.block(0, 0, matrixDim, matrixDim) = B;
            X.block(matrixDim, 0, matrixDim, matrixDim) = C;
            X.block(matrixDim, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
        } else {
            Y.block(matrixDim, matrixDim, matrixDim, matrixDim) = -1 * C;

            X.block(0, matrixDim, matrixDim, matrixDim) = C;
            X.block(0, 0, matrixDim, matrixDim) = B;
            X.block(matrixDim, 0, matrixDim, matrixDim) = C;
        }
    }

    /// Find the minimum positive eigenvalue of the quadratic eigenvalue problem \[At^2 + Bt + c\].
    /// First transform it to the generalized eigenvalue problem X+lY.
    /// If the updateOnly flag is false, compute matrices X,Y from scratch;
    /// otherwise only update them.
    /// \param[in] A Input matrix
    /// \param[in] B Input matrix
    /// \param[in] C Input matrix
    /// \param[in, out] X
    /// \param[in, out] Y
    /// \param[out] eigenvector The eigenvector corresponding to the minimum positive eigenvalue
    /// \param[in, out] updateOnly True if X,Y were previously computed and only B,C changed
    /// \return Minimum positive eigenvalue
    NT
    minimumPositiveEigenvalue(const MT &A, const MT &B, const MT &C, MT &X, MT &Y, VT &eigenvector, bool &updateOnly) {
        // perform linearization and create generalized eigenvalue problem X+lY
        linearization(A, B, C, X, Y, updateOnly);

        NT lambdaMinPositive = std::numeric_limits<NT>::max();
        int matrixDim = A.rows();

#if defined(EIGEN_EIGENVALUES_SOLVER)
        // use the Generalized eigenvalue solver of Eigen

        // compute generalized eigenvalues with Eigen solver
        Eigen::GeneralizedEigenSolver<MT> ges(X, -Y);

        // retrieve minimum positive eigenvalue
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();
        int index = 0;

        for (int i = 0; i < alphas.rows(); i++) {

            if (betas(i) == 0 || alphas(i).imag() != 0)
                continue;

            double lambda = alphas(i).real() / betas(i);
            if (lambda > 0 && lambda < lambdaMinPositive) {
                lambdaMinPositive = lambda;
                index = i;
            }
        }

        // retrieve corresponding eigenvector
        typename Eigen::GeneralizedEigenSolver<MT>::EigenvectorsType eivecs = ges.eigenvectors();
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivec = eivecs.col(index);

        eigenvector.resize(matrixDim);
        for (int i = 0; i < matrixDim; i++)
            eigenvector(i) =  eivec(matrixDim + i).real();
//        std::cout << lambdaMinPositive << " " << eigenvector.transpose() << "\n";fflush(stdout);
#else
        // Transform the problem to a standard eigenvalue problem and use the general eigenvalue solver of Spectra

        // This makes the transformation to standard eigenvalue problem. See class for more info.
        // We have the generalized problem  X + lY, or Xv = -lYv
        // create matrix M = -Y * X^[-1]
        Y *= -1;
        DenseProductMatrix<NT> M(&Y, &X);

        // This parameter is for Spectra. It must be larger than #(requested eigenvalues) + 2
        // and smaller than the size of matrix;
        int ncv = 3;

        // Prepare to solve Cx = (1/l)x
        // we want the smallest positive eigenvalue in the original problem,
        // so in this the largest positive eigenvalue;
        Spectra::GenEigsSolver<NT, Spectra::LARGEST_REAL, DenseProductMatrix<NT> > eigs(&M, 1, ncv);

        // compute
        eigs.init();
        eigs.compute();

        // reset matrix Y
        Y *= -1; // TODO find better way

        //retrieve result and invert to get required eigenvalue of the original problem
        if (eigs.info() != Spectra::SUCCESSFUL)
            return 0;

        lambdaMinPositive = 1/((eigs.eigenvalues())(0).real());

        // retrieve the eigenvector; we are at dimension 2*matrixDim,
        // while the eigenvector in the original quadratic problem had dimension MatrixDim.
        // We only need to keep the last #matrixDim coordinates.
        eigenvector.resize(matrixDim);
        for (int i = 0; i < matrixDim; i++)
            eigenvector(i) =  (eigs.eigenvectors()).col(0)(matrixDim + i).real();
#endif
//        std::cout << lambdaMinPositive << " " << eigenvector.transpose() << "\n";fflush(stdout);
        return lambdaMinPositive;
    }
};

#endif //VOLESTI_QUADRATICEIGENVALUESPROBLEM_H
