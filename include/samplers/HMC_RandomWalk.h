//
// Created by panagiotis on 2/19/20.
//

#ifndef HMCPOLYTOPES_HMC_RANDOMWALK_H
#define HMCPOLYTOPES_HMC_RANDOMWALK_H

template<class Point, typename RNGType>
class HMC_Settings {
public:
    typedef typename Point::FT NT;
    int walk_length;
    RNGType rng;
    Point c;
    NT temperature;
    NT diameter;

    HMC_Settings(int walkLength, RNGType rng, Point c, NT temperature, NT diameter) : walk_length(walkLength), rng(rng),
                                                                                      c(c), temperature(temperature),
                                                                                      diameter(diameter) {}
};

template<class Point, typename RNGType>
class HMC_Boltzmann {
public:
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic, Eigen::Dynamic> MT;

    HMC_Settings<Point, RNGType> settings;

    HMC_Boltzmann(HMC_Settings<Point, RNGType> & hmcSettings) : settings(hmcSettings) {}

    struct PrecomputedValues {
        VT ac;
        VT acPrime;
        Point minus_c_overT;
        Point minus_c_over2T;
        VT ac4;
        MT dotProducts;
        VT ap;

        PrecomputedValues(int dim, const Point& c, const NT temperature) {
            minus_c_over2T =  c/ (-2*temperature);
            minus_c_overT =  c/ (-temperature);
        }
    };

    void sample(HPolytope<Point>& polytope, Point& interiorPoint, const int pointsNum, std::list<Point>& points) {
        int n = polytope.dimension();
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, n - 1);
        RNGType& rng = settings.rng;
        PrecomputedValues precomputedValues(polytope.num_of_hyperplanes(), settings.c, settings.temperature);
        polytope.precomputeForHMC_Boltsmann(interiorPoint, settings.c, settings.temperature, precomputedValues.ac, precomputedValues.ac4, precomputedValues.acPrime, precomputedValues.ap, precomputedValues.dotProducts);

        for (unsigned int i = 1; i <= pointsNum; ++i) {
            for (unsigned int j = 0; j < settings.walk_length; ++j) {
                getNextPoint(polytope, interiorPoint, precomputedValues);
            }
            points.push_back(interiorPoint);
        }
    }


    void getNextPoint(HPolytope<Point> &polytope, Point &p, PrecomputedValues& precomputedValues) {

        /*************** INITIALIZE **************/

        RNGType &rng = settings.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        NT T = urdist(rng) * settings.diameter;
        const int n = polytope.dimension();
        const NT dl = 0.995;

        Point v = get_direction<RNGType, Point, NT>(n), p0 = p;
        std::pair<NT, int> pbpair;
        NT temperature = settings.temperature;
        int it = 0;
        VT av; // store the dot products of the facets of the polytope with direction v
        NT lambda_prev = NT(0);
        int facet;

        /************ THE FIRST BOUNDARY ORACLE CALL ************/

        pbpair = polytope.HMC_Boltzmann_intersect_first_call(p, av, v, precomputedValues.ac, precomputedValues.ap);

//        std::cout <<  p.getCoefficients()(0) << " " << p.getCoefficients()(1) <<"\n";
//        std::cout <<  v.getCoefficients()(0) << " " << v.getCoefficients()(1) <<"\n";
        std::cout << pbpair.first <<"\n";

        if (T <= pbpair.first) {
            precomputedValues.ap += precomputedValues.ac*(T*T) + av*T;
            p += T*T*precomputedValues.minus_c_over2T + T*v;
            return;
        }

        facet = pbpair.second;
        lambda_prev = dl * pbpair.first;

        p += (lambda_prev*lambda_prev)*precomputedValues.minus_c_over2T + lambda_prev*v;
        T -= lambda_prev;
        v = precomputedValues.minus_c_overT*lambda_prev + lambda_prev*v;
        polytope.compute_reflection(v, p, pbpair.second);

        /*************** THE REFLECTIONS ******************/
        // we need the first call and the reflections separately, to avoid unnecessary computations

        while (it<0.5*n) {
            pbpair = polytope.HMC_Boltzmann_intersect_not_first_call(p, v, precomputedValues.ac, precomputedValues.acPrime,
                    precomputedValues.ac4, av, precomputedValues.ap, lambda_prev, facet, precomputedValues.dotProducts);
            std::cout << pbpair.first <<"\n";
            if (T <= pbpair.first) {
                p += T*T*precomputedValues.minus_c_over2T + T*v;
                lambda_prev = T;
                precomputedValues.ap += precomputedValues.ac*(T*T) + av*T;
                break;
            }

            facet = pbpair.second;
            lambda_prev = dl * pbpair.first;
            p += (lambda_prev*lambda_prev)*precomputedValues.minus_c_over2T + lambda_prev*v;
            T -= lambda_prev;
            v = precomputedValues.minus_c_overT*lambda_prev + lambda_prev*v;
            polytope.compute_reflection(v, p, pbpair.second);

            it++;
        }

        precomputedValues.ap += precomputedValues.ac*(lambda_prev*lambda_prev) + av*lambda_prev;

//        if(it == 10*n) p = p0;
    }
};
#endif //HMCPOLYTOPES_HMC_RANDOMWALK_H
