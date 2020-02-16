//
// Created by panagiotis on 2/16/20.
//

#ifndef VOLESTI_RANDOM_WALK_MISC_H
#define VOLESTI_RANDOM_WALK_MISC_H



// Pick a random direction as a normilized vector
template <typename RNGType, typename Point, typename NT>
Point get_direction(const unsigned int dim) {

    boost::normal_distribution<> rdist(0,1);
    std::vector<NT> Xs(dim,0);
    NT normal = NT(0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType rng2 = var.rng;
    for (unsigned int i=0; i<dim; i++) {
        Xs[i] = rdist(rng);
        normal += Xs[i] * Xs[i];
    }
    normal=1.0/std::sqrt(normal);

    for (unsigned int i=0; i<dim; i++) {
        Xs[i] = Xs[i] * normal;
    }
    Point p(dim, Xs.begin(), Xs.end());
    return p;
}


// Pick a random point from a d-sphere
template <typename RNGType, typename Point, typename NT>
Point get_point_on_Dsphere(const unsigned int dim, const NT &radius){
    Point p = get_direction<RNGType, Point, NT>(dim);
    p = (radius == 0) ? p : radius * p;
    return p;
}


// Pick a random point from a d-ball
template <typename RNGType, typename Point, typename NT>
Point get_point_in_Dsphere(const unsigned int dim, const NT &radius){

    boost::random::uniform_real_distribution<> urdist(0,1);
    NT U;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng2(seed);
    Point p = get_direction<RNGType, Point, NT>(dim);
    U = urdist(rng2);
    U = std::pow(U, 1.0/(NT(dim)));
    p = (radius*U)*p;
    return p;
}

#endif //VOLESTI_RANDOM_WALK_MISC_H
