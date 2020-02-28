//
// Created by panagiotis on 2/23/20.
//

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Spectra/GenEigsSolver.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "samplers.h"
#include "rounding.h"
#include "sample_only.h"
#include "EigenvaluesProblems.h"
#include "LMI.h"
#include "spectrahedron.h"
#include "HMC_RandomWalk.h"
#include "sdp_generator.h"
#include "SDPA_FormatManager.h"
#include "CoordinateDirectionsHitAndRun_RandomWalk.h"
#include "SimulatedAnnealing.h"

/// An enumeration of the available random walks
enum RandomWalk {
    coordHitAndRun, HMC_Boltzmann
};


/// Find out which random walk to use
/// \param[in] _randomWalk Input from user from R
/// \return A member of enum RandomWalk indicating which random walk to use
RandomWalk whichRandomWalk(Rcpp::Nullable<Rcpp::CharacterVector> const & _randomWalk) {
    if (Rcpp::as<std::string>(_randomWalk).compare(std::string("CDHR")) == 0)
        return coordHitAndRun;
    else if (Rcpp::as<std::string>(_randomWalk).compare(std::string("HMC")) == 0)
        return HMC_Boltzmann;
    else
        throw Rcpp::exception("Unknown walk type!");
}


/// Loads a spectrahedron from a file, in SDPA format
/// \tparam Spectrahedron Must be a spectrahedron
/// \tparam NT Numeric type
/// \tparam MT Matrix type
/// \tparam VT Vector type
/// \param[out] spectrahedron Return object
/// \param[out] objectiveFunction
/// \param[in] inputFile Name of the file containing spectrahedron data
template<class Spectrahedron, typename NT, typename MT,  typename VT>
void loadSpectrahedronFromFile(Spectrahedron & spectrahedron, VT& objectiveFunction, Rcpp::Nullable<Rcpp::CharacterVector> const & inputFile) {
    LMI<NT, MT, VT> lmi;

    if (!inputFile.isNotNull())
        throw Rcpp::exception("Wrong input!");

    std::ifstream inp;
    inp.open(Rcpp::as<std::string>(inputFile), std::ifstream::in);

    loadSDPAFormatFile(inp, lmi, objectiveFunction);
    spectrahedron = Spectrahedron(lmi);
}


/// Samples the spectrahedron uniformly using coordinate directions hit and run
/// \tparam Spectrahedron
/// \tparam Point
/// \tparam MT Matrix type
/// \tparam VT Vector type
/// \param[in] spectrahedron A spectrahedron
/// \param[in] interiorPoint A point in the interior of the spectrahedron
/// \param[in] numPoints The number of points to sample
/// \param[in] walkLength The length of the walk
/// \param[out] randPoints The list which will hold the samples
template<class Spectrahedron, class Point, typename MT, typename VT>
void sampleCDHR(Spectrahedron& spectrahedron, Point& interiorPoint, int numPoints, int walkLength, std::list<Point>& randPoints) {
    typedef CoordinateDirectionsHitAndRun_RandomWalk<Point, MT, VT, RNGType > CDHR;
    typedef boost::mt19937 RNGType;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    CDHR cdhr(walkLength, rng);
    cdhr.sample(spectrahedron, interiorPoint, numPoints, randPoints);
}


/// Prepare settings for the HMC Boltzmann random walk, i.e. e^{-cx/T}
/// \tparam Spectrahedron
/// \tparam HMCSettings
/// \tparam Point
/// \tparam NT Numeric Type
/// \tparam VT Vector Type
/// \param[in] A spectrahedron
/// \param[out] settings The settings for the HMC random walk
/// \param[in] objectiveFunction The direction c
/// \param[in] walkL The length of the random walk
/// \param[in] parameters Input from user from R. Must contain T and diameter
/// \param[in] interiorPoint A point in the interior of the spectrahedron
template<class Spectrahedron, typename HMCSettings, class Point, typename NT, typename VT>
void getHMCSettings(Spectrahedron& spectrahedron, HMCSettings &settings, VT const & objectiveFunction, const int walkL, Rcpp::Nullable<Rcpp::List> const & parameters, Point const & interiorPoint) {
    typedef boost::mt19937 RNGType;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    // retrieve T
    if (!Rcpp::as<Rcpp::List>(parameters).containsElementNamed("T"))  throw Rcpp::exception("Wrong input - Give T!");
    NT T = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["T"]);

    // retrieve diameter
    NT diameter;
    if (!Rcpp::as<Rcpp::List>(parameters).containsElementNamed("Diameter")) {
        // no diameter, estimate it
        diameter = spectrahedron.estimateDiameter(spectrahedron.dimension()*4, interiorPoint);
        std::cout << "No diameter was provided; estimated it at " << diameter << "\n";fflush(stdout);
    }
    else
        // get diameter from user
        diameter = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["Diameter"]);

    Point c(objectiveFunction);
    settings = HMCSettings(walkL, rng, c, T, diameter);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix HMC(Rcpp::Nullable<Rcpp::CharacterVector> spectrahedronInputFile = R_NilValue,
        Rcpp::Nullable<Rcpp::CharacterVector> random_Walk = R_NilValue,
        Rcpp::Nullable<unsigned int> N = R_NilValue,
        Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
        Rcpp::Nullable<Rcpp::List> parameters = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
    typedef HMC_RandomWalk<Point, MT, VT, RNGType > HMC_RandomWalk;
    typedef CoordinateDirectionsHitAndRun_RandomWalk<Point, MT, VT, RNGType > CDHR;

    // find out which random walk we will use
    RandomWalk randomWalk = whichRandomWalk(random_Walk);

    // load spectrahedron from file
    VT objectiveFunction;
    SPECTRAHEDRON spectrahedron;
    loadSpectrahedronFromFile<SPECTRAHEDRON, NT, MT, VT> (spectrahedron, objectiveFunction, spectrahedronInputFile);

    // a list to store the samples
    std::list<Point> randPoints;

    // get walk information
    int walkLength, numPoints;
    if(!walk_length.isNotNull() || !N.isNotNull()) throw Rcpp::exception("Wrong input!");
    walkLength = Rcpp::as<unsigned int>(walk_length);
    numPoints = Rcpp::as<unsigned int>(N);

    // get interiorPoint
    Point interiorPoint(spectrahedron.dimension());

    // sample points
    switch (randomWalk) {
        case coordHitAndRun:
            // Sample from uniform distribution with coordinate directions hit and run
            sampleCDHR<SPECTRAHEDRON, Point, MT, VT>(spectrahedron, interiorPoint, numPoints, walkLength, randPoints);
            break;
        case HMC_Boltzmann:
            // Sample from Boltzmann distribution with HMC
            HMC_RandomWalk::Settings settings;
            getHMCSettings<SPECTRAHEDRON, HMC_RandomWalk::Settings, Point, NT, VT>(spectrahedron, settings, objectiveFunction, walkLength, parameters, interiorPoint);
            HMC_RandomWalk hmc(settings);
            hmc.sample(spectrahedron, interiorPoint, numPoints, randPoints);
            break;
    }


    // return points to R
    unsigned int n = spectrahedron.dimension();
    MT RetMat(n, numPoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
        RetMat.col(jj) = (*rpit).getCoefficients();

    return Rcpp::wrap(RetMat);

}



//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix solveSDP(Rcpp::Nullable<Rcpp::CharacterVector> spectrahedronInputFile = R_NilValue,
                        Rcpp::Nullable<Rcpp::List> parameters = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
    typedef SimulatedAnnealing<Point, MT, VT> SA;

    // load spectrahedron from file
    VT objectiveFunction;
    SPECTRAHEDRON spectrahedron;
    loadSpectrahedronFromFile<SPECTRAHEDRON, NT, MT, VT> (spectrahedron, objectiveFunction, spectrahedronInputFile);

    NT error = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["error"]);
    int walkLength = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["walkLength"]);
    int maxNumSteps = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["maxNumSteps"]);
    NT k = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["k"]);

    SA::Settings settings(error, walkLength, maxNumSteps, k);
    SA simulatedAnnealing(&spectrahedron, Point(objectiveFunction), settings);

    Point x;
    simulatedAnnealing.solve(x);

    return Rcpp::wrap(x.getCoefficients());
}