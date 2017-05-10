#ifndef _REFINER_H_
#define _REFINER_H_

#include "carchitect.h"

ILOSTLBEGIN


class Refiner : public CArchitect
{
public:
    Refiner(const InputInstance& inputInstance,
            const DoubleMatrix& M,
            const IntMatrix &e,
            const int Z,
            const unsigned int k,
            const bool rootNotFixed,
            const bool forceDiploid,
            const double bound_bar_f);

    double getDistance();

protected:
    /// Upper bound to the distance
    const double _bound_bar_f;

protected:
    /// Build constraints
    void buildConstraints();
    /// Build objective
    void buildObjective();

};

#endif // _REFINER_H_
