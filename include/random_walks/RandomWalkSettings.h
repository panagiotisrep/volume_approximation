//
// Created by panagiotis on 2/16/20.
//

#ifndef VOLESTI_RANDOMWALKSETTINGS_H
#define VOLESTI_RANDOMWALKSETTINGS_H

template<class Point, typename RNGType>
class RandomWalkSettings {
public:
    typedef typename Point::FT NT;
    typedef typename HPolytope<Point>::VT VT;

    int walk_len, dim;
    RNGType rng;
    VT lamdas, Av;

    RandomWalkSettings(const int walk_len, const int dim, RNGType& rng) {
        this->walk_len = walk_len;
        this->rng = rng;
        this->dim = dim;
    }
};

#endif //VOLESTI_RANDOMWALKSETTINGS_H
