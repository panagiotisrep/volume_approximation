//
// Created by panagiotis on 2/23/20.
//

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "samplers.h"
#include "rounding.h"
#include "sample_only.h"
#include "EigenvaluesProblems.h"
#include "spectrahedron.h"
#include "LMI.h"
#include "sdp_generator.h"

template <typename NT, typename MT, typename VT>
void writeSDPAFormatFile(std::ostream &os, const LMI<NT, MT, VT> &lmi, const VT &objectiveFunction) {
    int dim = lmi.dimension();

    std::vector<MT> matrices = lmi.getMatrices();
    MT A0 = matrices[0];
    matrices.erase(matrices.begin());

    os << dim << "\n";
    os << 1 << "\n";
    os << A0.rows() << "\n";

    os << objectiveFunction.transpose() << "\n";

    for (int i = 0; i < A0.rows(); i++)
        os << A0.row(i) << "\n";

    for (MT matrix : matrices)
        for (int i = 0; i < matrix.rows(); i++)
            os << -1 * matrix.row(i) << "\n";
}

//' @export
// [[Rcpp::export]]
void generator_sdp(Rcpp::Nullable<int> nn = R_NilValue,
                   Rcpp::Nullable<int> mm = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI <NT, MT, VT> LMI;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
    unsigned int n = Rcpp::as<int>(nn);

    SPECTRAHEDRON spectrahedron;
    spectrahedron = generateSDP2<NT>(Rcpp::as<int>(nn), Rcpp::as<int>(mm));

    Point c = get_direction<RNGType, Point, NT>(n);

    std::filebuf fb;
    std::string bar = "_";
    std::string txt = ".txt";
    fb.open("sdp_"+bar+std::to_string(n)+bar+std::to_string(Rcpp::as<int>(mm))+txt, std::ios::out);
    std::ostream os(&fb);
    writeSDPAFormatFile<NT, MT, VT>(os, spectrahedron.getLMI(), c.getCoefficients());

    return;
}