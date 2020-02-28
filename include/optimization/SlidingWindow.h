//
// Created by panagiotis on 2/28/20.
//

#ifndef VOLESTI_SLIDINGWINDOW_H
#define VOLESTI_SLIDINGWINDOW_H


/// Computes the relative error
/// \param[in] approx The approximated value
/// \param[in] exact The exact value
/// \return The relative error
double relativeError(double approx, double exact) {
    return std::fabs((exact - approx) / exact);
}


/// A sliding window, which allows to get the relative error of approximations
/// \tparam NT Numeric type
template<typename NT>
class SlidingWindow {
public:

    /// The stored approximations
    std::list<NT> approximations;
    /// The size of the window
    int windowSize;
    /// The number of entries in list approximations
    int numEntries;

    /// Constructor
    /// \param[in] windowSize The size of the window
    SlidingWindow(int windowSize) : windowSize(windowSize) {
        count = 0;
    }

    /// Adds an approximation in the window
    /// \param[in] approximation The new approximation
    void push(NT approximation) {
        // if window is full, remove the oldest value
        if (count >= windowSize) {
            approximations.pop_back();
        }
        else
            count++;

        approximation.push_front(approximation);
    }

    /// \return The relative error between the youngest and oldest approximations
    double getRelativeError() {
        if (count < windowSize)
            return 1;

        return relativeError(approximations.back(), approximations.front());
    }
};

#endif //VOLESTI_SLIDINGWINDOW_H
