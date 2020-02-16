//
// Created by panagiotis on 2/16/20.
//

#ifndef VOLESTI_RANDOMWALK_H
#define VOLESTI_RANDOMWALK_H

//#include "RandomWalkSettings.h"

template<class Point, typename RNGType>
class RandomWalk {
public:
    virtual ~RandomWalk(){}
    virtual void sample(ConvexBody<Point, typename Point::FT>&, Point&, const int, std::list<Point>&, RandomWalkSettings<Point, RNGType>&) = 0;
};

#endif //VOLESTI_RANDOMWALK_H
