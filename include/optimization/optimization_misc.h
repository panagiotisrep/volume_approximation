//
// Created by panagiotis on 2/12/20.
//

#ifndef VOLUME_APPROXIMATION_MISC_H
#define VOLUME_APPROXIMATION_MISC_H

#include "polytopes.h"
#include "Eigen"
#include <list>
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>
#include <queue>

namespace optimization {

    /**
     *  CONSTANTS
     */
    const double ZERO = 0.000000000001;

    /**
     * TYPEDEFS
     */
    /**
    * For the Eigen library
    */
    typedef double NT_MATRIX;
    typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, 1> VT;

    /**
     * VARIABLES
     */
    int STEPS;



    /**
    * Computes the relative error
    * @param approx
    * @param exact
    * @return
    */
    double relative_error(double approx, double exact) {
        return fabs((exact - approx) / exact);
    }


    /**
     * A class to save successive approximation values
     */
    class SlidingWindow {
    public:
        std::queue<double> evals;
        int evalsSize; // how many values to store
        int count;
        double min;

        SlidingWindow(int evalSize) {this->evalsSize = evalSize;count = 0;}

        void push(double eval) {
            if (count >= evalsSize) {
                evals.pop();
            }
            else
                count++;

            evals.push(eval);
            min = eval;
        }

        double getRelativeError() {
            if (count < evalsSize)
                return 1;

            return fabs((min - evals.front()) / min);
        }

        double getRelativeErrorLastRound(double min) {
            return fabs((min - evals.back()) / min);
        }


        void half() {
            for (int i=0 ; i<evals.size()/2 ; i++) {
                evals.pop();
                count--;
            }
        }
    };





    /**
     * Returns the index of the minimum element in the given vector
     * @param vec
     * @return
     */
    int minElementIndex(VT& vec) {
        int index = 0;
        double min = vec(0);

        for (int i=1 ; i<vec.rows() ; i++)
            if (vec(i) < min) {
                index = i;
                min = vec(i);
            }

        return index;
    }

    /**
     * Computes the euclidean distance between two points.
     *
     * @param v1
     * @param v2
     * @return
     */
    double euclideanDistance(VT& v1, VT& v2) {
        double sum = 0;

        for (int i=0 ; i<v1.rows() ; i++)
            sum += (v1(i) - v2(i))*(v1(i) - v2(i));

        return sqrt(sum);
    }


    template <class Point>
    double euclideanDistance(const Point& v1, const Point& v2) {
        double sum = 0;

        for (int i=0 ; i<v1.dimension() ; i++)
            sum += (v1[i] - v2[i])*(v1[i] - v2[i]);

        return sqrt(sum);
    }

    template<class Point>
    void getArithmeticMean(std::vector<Point> points, Point& mean) {
        int dim = points[0].dimension();
        mean = Point(dim);

        for (auto p : points)
            mean = mean + p;

        mean = mean / points.size();
    }

    template<class Point>
    void getArithmeticMean(std::list<Point> points, Point& mean) {
        int dim = points.front().dimension();
        mean = Point(dim);

        for (auto p : points)
            mean = mean + p;

        mean = mean / points.size();
    }

}

#endif //VOLUME_APPROXIMATION_MISC_H
