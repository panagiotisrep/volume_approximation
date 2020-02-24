//
// Created by panagiotis on 2/22/20.
//

#ifndef VOLESTI_LMI_H
#define VOLESTI_LMI_H

/// This class handles a linear matrix inequality of the form \[A_0 +  \sum x_i A_i\]
/// @tparam NT Numeric Type
/// @tparam MT Matrix Type
/// @tparam VT Vector Type
template<typename NT, typename MT, typename VT>
class LMI {
    /// The matrices A_0, A_i
    std::vector<MT> matrices;

    /// The dimension of the vector x
    unsigned int d;

    /// The size of the matrices A_i
    unsigned int m;


    /// Creates A LMI object
    /// \param[in] matrices The matrices A_0, A_i
    LMI(std::vector<MT>& matrices) {
        typename std::vector<MT>::iterator it = matrices.begin();

        while (it!=matrices.end()) {
            this->matrices.push_back(*it);
        }

        d = matrices.size() - 1;
        m = matrices[0].rows();
    }

    /// \returns The dimension of vector x
    unsigned int dimension() const {
        return d;
    }

    /// \returns The size of the matrices
    unsigned int sizeOfMatrices() const {
        return m;
    }

    /// Evaluate \[A_0 + \sum x_i A_i \]
    /// \param[in] x The input vector
    /// \param[out] ret The output matrix
    void evaluate(const VT& x, MT& ret) const {
    }

    /// Compute  \[x_1*A_1 + ... + x_n A_n]
    /// \param[in] x Input vector
    /// \param[out] res Output matrix
    void evaluateWithoutA0(const VT& x, MT& ret) const {

    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] p Input parameter
    /// \param[in] Input vector: lmi(p)*e = 0, e != 0
    /// \param[out] ret The normalized gradient of the determinant of the LMI at p
    void normalizedDeterminantGradient(const VT& p, const VT& e, VT& ret) {}
};


/// This class handles a linear matrix inequality of the form \[A_0 +  \sum x_i A_i\]
/// A template specialization for dense Eigen matrices and vectors
/// @tparam NT Numeric Type
template<typename NT>
class LMI<NT, Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Matrix<NT,Eigen::Dynamic,1> > {
    public:
    /// Eigen matrix type
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    /// Eigen vector type
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    /// The matrices A_0, A_i
    std::vector<MT> matrices;

    /// The dimension of the vector x
    unsigned int d;

    /// The size of the matrices A_i
    unsigned int m;

    LMI(){}

    /// Creates A LMI object
    /// \param[in] matrices The matrices A_0, A_i
    LMI(std::vector<MT>& matrices) {
        typename std::vector<MT>::iterator it = matrices.begin();

        while (it!=matrices.end()) {
            this->matrices.push_back(*it);
            it++;
        }

        d = matrices.size() - 1;
        m = matrices[0].rows();
    }

    /// \returns The dimension of vector x
    unsigned int dimension() const {
        return d;
    }

    /// \return The matrices A0, A1, ..., Ad
    std::vector<MT> getMatrices() const {
        return matrices;
    }

    /// \returns The size of the matrices
    unsigned int sizeOfMatrices() const {
        return m;
    }

    /// Evaluate A_0 + \[A_0 + \sum x_i A_i \]
    /// \param[in] x The input vector
    /// \param[out] ret The output matrix
    void evaluate(const VT& x, MT& ret) {
        typename std::vector<MT>::iterator it = matrices.begin();

        ret = *it;
        it++;
        int i=0;

        while (it != matrices.end()) {
            ret += *it * x(i);
            it++;
            i++;
        }
    }

    /// Compute  \[x_1*A_1 + ... + x_n A_n]
    /// \param[in] x Input vector
    /// \param[out] res Output matrix
    void evaluateWithoutA0(const VT& x, MT& res)  {
        res = MT::Zero(m, m);
        typename std::vector<MT>::iterator it;

        int i = 0;
        it = matrices.begin();
        ++it; // skip A0
        for (; it != matrices.end(); it++, i++)
            res.noalias() += x(i) * (*it);
    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] p Input parameter
    /// \param[in] Input vector: lmi(p)*e = 0, e != 0
    /// \param[out] ret The normalized gradient of the determinant of the LMI at p
    void normalizedDeterminantGradient(const VT& p, const VT& e, VT& ret) {
        ret.resize(d);

        // i-th coordinate of the determinant is e^T * A_i * e
        for (int i = 0; i < d; i++) {
            ret(i) = e.dot(matrices[i+1] * e);
        }

        ret.normalize();
    }
};

#endif //VOLESTI_LMI_H
