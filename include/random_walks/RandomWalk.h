//
// Created by panagiotis on 2/16/20.
//

#ifndef VOLESTI_RANDOMWALK_H
#define VOLESTI_RANDOMWALK_H

#include "RandomWalkSettings.h"

template<class Point, typename RNGType>
class RandomWalk {
public:
    RandomWalkSettings<Point, RNGType>& settings;

    virtual ~RandomWalk(){}
    RandomWalk(RandomWalkSettings<Point, RNGType>& settings) : settings(settings) {};
    virtual void sample(ConvexBody<Point, typename Point::FT>&, Point&, const int, std::list<Point>&) = 0;
};

#endif //VOLESTI_RANDOMWALK_H
