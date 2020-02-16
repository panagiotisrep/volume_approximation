//
// Created by panagiotis on 2/16/20.
//

#ifndef VOLESTI_CONVEXBODY_H
#define VOLESTI_CONVEXBODY_H

//#include "Eigen"

template<class Point, typename NT>
class ConvexBody {
public:
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    virtual std::pair<NT, NT> line_intersect(Point&, Point&, VT&, VT&, bool) = 0;
    virtual std::pair<NT, NT> line_intersect(Point&, Point&, VT&, VT& ,const NT&, bool)  = 0;

 
    virtual unsigned int dimension() const = 0;
    virtual int num_of_hyperplanes() const = 0;
    virtual std::pair<NT, NT> vline_intersect(Point& p, Point& v , VT& lamdas, VT& Av, bool pos = false) {
        return line_intersect(p, v , lamdas, Av, pos);

    };
    virtual std::pair<NT, NT> vline_intersect(Point& p, Point& v , VT& lamdas, VT& Av, const NT& lamda, bool pos = false) {
        return line_intersect(p, v , lamdas, Av, lamda, pos);

    };
};

#endif //VOLESTI_CONVEXBODY_H
