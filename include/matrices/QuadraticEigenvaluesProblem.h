//
// Created by panagiotis on 2/22/20.
//

#ifndef VOLESTI_QUADRATICEIGENVALUESPROBLEM_H
#define VOLESTI_QUADRATICEIGENVALUESPROBLEM_H

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
class QuadraticEigenvaluesProblem<NT, Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Matrix<NT,Eigen::Dynamic,1> >{
public:
    /// The type for Eigen Matrix
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

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
    void linearization(const MT& A, const MT& B, const MT& C, MT& X, MT& Y, bool& updateOnly) {
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
        }
        else {
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
    NT minimumPositiveEigenvalue(const MT& A, const MT& B, const MT& C, MT& X, MT& Y, VT& eigenvector, bool& updateOnly) {
        // perform linearization and create generalized eigenvalue problem X+lY
        linearization(A, B, C, X, Y, updateOnly);

        NT lambdaMinPositive = std::numeric_limits<NT>::max();
        int matrixDim = A.rows();

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

        return lambdaMinPositive;
    }

};

#endif //VOLESTI_QUADRATICEIGENVALUESPROBLEM_H
