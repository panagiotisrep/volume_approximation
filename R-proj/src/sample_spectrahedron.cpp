//
// Created by panagiotis on 2/23/20.
//

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "samplers.h"
#include "rounding.h"
#include "sample_only.h"
#include "QuadraticEigenvaluesProblem.h"
#include "LMI.h"
#include "spectrahedron.h"
#include "HMC_RandomWalk.h"
#include "sdp_generator.h"
#include "SDPA_FormatManager.h"

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix HMC(Rcpp::Nullable<Rcpp::Reference> SP = R_NilValue,
        Rcpp::Nullable<Rcpp::CharacterVector> file = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> objectiveFunction = R_NilValue,
                                  Rcpp::Nullable<unsigned int> N = R_NilValue,
                                  Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                                  Rcpp::Nullable<double> Temperature = R_NilValue,
                                    Rcpp::Nullable<double> Diameter = R_NilValue){


    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
    typedef HMC_RandomWalk<Point, MT, VT, RNGType > HMC_RandomWalk;

    LMI<NT, MT, VT> lmi;
    VT c;

    if (SP.isNotNull()) {
        if (!objectiveFunction.isNotNull())
            throw Rcpp::exception("Wrong input!");

        int dimension = Rcpp::as<Rcpp::Reference>(SP).field("dimension");

        c.setZero(dimension);
        std::vector<NT> obj = Rcpp::as<std::vector<NT> >(objectiveFunction);
        for (unsigned int j=0; j<dimension; j++){
            c(j) = obj[j];
        }

        std::vector<MT> matrices = Rcpp::as<Rcpp::Reference>(SP).field("matrices");
        lmi = LMI<NT, MT, VT>(matrices);
    } else if (file.isNotNull()) {
        std::ifstream inp;
        inp.open(Rcpp::as<std::string>(file), std::ifstream::in);

        loadSDPAFormatFile(inp, lmi, c);
    } else
        throw Rcpp::exception("Wrong input!");

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    SPECTRAHEDRON spectrahedron(lmi);
    unsigned int n = spectrahedron.dimension();

    RNGType rng(seed);

    int walkL, NN;
    NT T, diam;
    if(walk_length.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_length);
    if(N.isNotNull()) NN = Rcpp::as<unsigned int>(N);
    if(Temperature.isNotNull()) T = Rcpp::as<NT>(Temperature);
    if(Diameter.isNotNull()) diam = Rcpp::as<NT>(Diameter);

//    SP.ComputeInnerBall(diam, radius);

    std::list<Point> randPoints;
    
    Point _c(c);
    HMC_RandomWalk::Settings settings(walkL, rng, _c, T, diam);
    HMC_RandomWalk hmc(settings);
    Point p(n);
    hmc.sample(spectrahedron, p, NN, randPoints);


    MT RetMat(n, NN);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
        RetMat.col(jj) = (*rpit).getCoefficients();

    return Rcpp::wrap(RetMat);


}