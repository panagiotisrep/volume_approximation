//
// Created by panagiotis on 2/24/20.
//

#include "Eigen/Eigen"
#define VOLESTI_DEBUG
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "rotating.h"
#include "misc.h"
#include "linear_extensions.h"
#include "cooling_balls.h"
#include "cooling_hpoly.h"
#include "sample_only.h"
#include "exact_vols.h"
#include "EigenvaluesProblems.h"
#include "spectrahedron.h"
#include "EigenDenseMatrix.h"
#include "DenseProductMatrix.h"
#include "HMC_RandomWalk.h"
#include "SDPA_FormatManager.h"



int main(int argc, char* argv[]) {


    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
    typedef HMC_RandomWalk<Point, MT, VT, RNGType > HMC_RandomWalk;

    std::ifstream inp;
    inp.open(argv[1],std::ifstream::in);
    LMI<NT, MT, VT> lmi;
    VT c;
    loadSDPAFormatFile(inp, lmi, c);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    SPECTRAHEDRON spectrahedron(lmi);
    unsigned int n = spectrahedron.dimension();

    RNGType rng(seed);

    int walkL=1, NN=10, T=1, diam=1;

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

    std::cout << RetMat;

    return 0;
}