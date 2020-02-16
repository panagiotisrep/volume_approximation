//
// Created by panagiotis on 2/16/20.
//

#ifndef VOLESTI_RANDOMDIRECTIONSHITANDRUN_H
#define VOLESTI_RANDOMDIRECTIONSHITANDRUN_H

//#include <random/uniform_real_distribution.hpp>
//#include <random/uniform_int_distribution.hpp>
#include <ConvexBody.h>
#include <list>
#include "RandomWalk.h"
//#include "list"
//#include "RandomDirectionsHitAndRunSettings.h"
//#include "ConvexBody.h"
//#include "Eigen"
//#include "random_walk_misc.h"

template<class Point, typename RNGType>
class RandomDirectionsHitAndRun : public RandomWalk<Point,  RNGType> {
public:

    void sample(ConvexBody<Point, typename Point::FT> &convexBody, Point &point, const int numOfPoints, std::list <Point> &points,
                RandomWalkSettings<Point, RNGType>& settings) override {

        typedef typename Point::FT NT;
        unsigned int n = convexBody.dimension();
        RNGType &rng = settings.rng;

        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, n - 1);

        unsigned int rand_coord, rand_coord_prev;
        NT kapa, lambda;
        Point p_prev = point, p1(n), p2(n), v(n);

        v = get_direction<RNGType, Point, NT>(n);
        std::pair<NT, NT> bpair = convexBody.vline_intersect(point, v, settings.lamdas, settings.Av);
        lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        point += (lambda * v);


        for (unsigned int i = 1; i <= numOfPoints; ++i) {
            for (unsigned int j = 0; j < settings.walk_len; ++j) {
                v = get_direction<RNGType, Point, NT>(n);
                std::pair <NT, NT> bpair = convexBody.vline_intersect(point, v, settings.lamdas, settings.Av, lambda);
                lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
                point = (lambda * v) + point;
            }
            points.push_back(point);
        }
    }

};

#endif //VOLESTI_RANDOMDIRECTIONSHITANDRUN_H
