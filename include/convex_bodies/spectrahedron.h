//
// Created by panagiotis on 2/22/20.
//

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "LMI.h"

/// This class manipulates a spectrahedron, described by a LMI
/// \tparam NT Numeric Type
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
template<typename NT, typename MT, typename VT>
class Spectrahedron {
public:

    /// Among successive calls of this class methods, we may need to pass data
    /// from one call to the next, to avoid repeating computations, or to efficiently update values
    struct PrecomputedValues {

        /// These flags indicate whether the corresponding matrices are computed
        /// if yes, we can use them and not compute them fro scratch
        bool computed_A = false;
        bool computed_C = false;
        bool computed_XY = false;

        /// The matrices the method positive_intersection receives from its previous call
        /// if the flag first_positive_intersection is true
        MT A, B, C, X, Y;

        /// In method positive_intersect, the distance we are computing corresponds
        /// to the minimum positive eigenvalue of a quadratic eigenvalue problem.
        /// This will hold the eigenvector for that eigenvalue
        VT eigenvector;

        /// Sets all flags to false
        void resetFlags() {
            computed_XY = computed_C = computed_A = false;
        }
    };

    /// The dimension of the spectrahedron
    unsigned int d;

    /// The linear matrix inequality that describes the spectrahedron
    LMI<NT, MT, VT> lmi;

    Spectrahedron() {}

    /// Creates a spectrahedron
    /// \param[in] lmi The linear matrix inequality that describes the spectrahedron
    Spectrahedron(const LMI<NT, MT, VT>& lmi) : lmi(lmi) {
        d = lmi.dimension();
    }


    /// Construct the quadratic eigenvalue problem \[At^2 + Bt + C \] for positive_intersect.
    /// A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// \param[in] a Input vector
    /// \param[in] b Input vector
    /// \param[in] c Input vector
    /// \param[in, out] precomputedValues Holds matrices A, C
    void createMatricesForPoistiveIntersection(const VT& a, const VT& b, const VT& c, PrecomputedValues& precomputedValues) {
        // check if matrices A, C are ready
        // if not compute them
        if (!precomputedValues.computed_A) {
            lmi.evaluateWithoutA0(a, precomputedValues.A);
        }

        if (!precomputedValues.computed_C) {
            lmi.evaluate(c, precomputedValues.C);
        }

        // compute Matrix B
        lmi.evaluateWithoutA0(b, precomputedValues.B);
    }

    /// Computes the distance d we must travel on the parametrized polynomial curve \[at^2 + bt + c \],
    /// assuming we start at t=0, and we start increasing t.
    /// We construct the quadratic eigenvalue problem \[At^2 + Bt + C \],
    /// where A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// Then we do a linearization and transform it to the generalized problem X+lY,
    /// which we pass to an external library.
    /// \param[in] a Input vector, the coefficient of t \[t^2\]
    /// \param[in] b Input vector, the coefficient of t
    /// \param[in] c Input Vector, the constant term
    /// \param[in, out] precomputedValues Data we move between successive method calls
    /// \returns The distance d
    NT positive_intersection(const VT& a, const VT& b, const VT& c, PrecomputedValues& precomputedValues) {
        unsigned int matrixDim = lmi.sizeOfMatrices();

        // create matrices A, B, C
        createMatricesForPoistiveIntersection(a, b, c, precomputedValues);

        // get the minimum positive eigenvalue of At^2 + Bt + C
        QuadraticEigenvaluesProblem<NT, MT, VT> quadraticEigenvaluesProblem;
        NT distance = quadraticEigenvaluesProblem.minimumPositiveEigenvalue(precomputedValues.A, precomputedValues.B, precomputedValues.C, precomputedValues.X,
                precomputedValues.Y, precomputedValues.eigenvector, precomputedValues.computed_XY);

        return distance;
    }

    /// Computed the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] point A point on the boundary of the spectrahedron
    /// \param[in] incomingDirection The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    /// \param[in] precomputedValues Must contain an eigenvalue needed to compute the reflection
    void compute_reflection(const VT& point, const VT& incomingDirection, VT& reflectedDirection, PrecomputedValues& precomputedValues) {

        // get the gradient of the determinant of the lmi at point
        VT grad;
        lmi.normalizedDeterminantGradient(point, precomputedValues.eigenvector, grad);

        // compute reflected direction
        // if v is original direction and s the surface normal,
        // reflected direction = v - 2 <v,s>*s

        NT dot = 2 * incomingDirection.dot(grad);
        reflectedDirection = incomingDirection - dot * grad;
    }

    /// \return The dimension of the spectrahedron
    unsigned int dimension() {
        return d;
    }

    /// \return The LMI describing this spectrahedron
    LMI<NT, MT, VT> getLMI() {
        return lmi;
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
