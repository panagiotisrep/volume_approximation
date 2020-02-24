//
// Created by panagiotis on 2/24/20.
//

#ifndef VOLESTI_DENSEPRODUCTMATRIX_H
#define VOLESTI_DENSEPRODUCTMATRIX_H

/// A wrapper class for dense Eigen matrices in Spectra
/// This class will be the wrapper to use the Spectra nonsymemmetric standard eigenvalue Cx = lx solver to
/// solve a generalized eigenvalue Ax = lBx.
/// In particular, this class represents the product @f[ C = B^-1 A @f]
///
/// \tparam NT Numeric Type
template<typename NT>
class DenseProductMatrix {
public:
    /// Eigen matrix type
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    /// Eigen vector type
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// The number of rows
    int _rows;
    /// The number of cols
    int _cols;

    /// Pointer to matrix A
    MT *A;
    /// Pointer to matrix B
    MT *B;

    /// The decomposition we will use
    typedef Eigen::FullPivLU<MT> Decomposition;

    /// The LU decomposition of B
    Decomposition Blu;

    /// Constructs an object of this class and computes the LU decomposition of B.
    ///
    /// \param[in] A The matrix A
    /// \param[in] B The matrix B
    DenseProductMatrix(MT *A, MT *B) : A(A), B(B) {
        Blu = Decomposition(*B);
        _rows = A->rows();
        _cols = B->cols();
    }

    ///Required by Spectra
    /// \return The number of rows
    int rows() {
        return _rows;
    }

    ///Required by Spectra
    /// \return The number of columns
    int cols() {
        return _cols;
    }

    /// Required by Spectra.
    /// Computed the product Cx = y, i.e. @f[ (B^-1 A)v = y@$]. But B = LU, so Ax = LUy.
    /// Let Ax = v, then LUy = v. Then Lw = v and finally Uy = w to get y;
    /// \param[in] x_in
    /// \param[out] y_out
    void perform_op(NT const * x_in, NT* y_out) {

        // Declaring the vector like this, we don't copy the values of v and after to w
        Eigen::Map<VT> const x(const_cast<double*>(x_in), _rows);
        VT const v = *A * x;

        Eigen::Map<VT> y(y_out, _rows);
        y = Blu.solve(v);
    }
};
#endif //VOLESTI_DENSEPRODUCTMATRIX_H
