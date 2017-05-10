#ifndef _CARCHITECT_H_
#define _CARCHITECT_H_

#include "basecarchitect.h"

ILOSTLBEGIN

/// This class represent the architect that builds the ILP model for solving the
/// problem when the given input is composed of a collection of leaves and an
/// integer e representing the maximum copy-number.
class CArchitect : public BaseCArchitect
{
public:
    CArchitect(const InputInstance& inputInstance,
               const DoubleMatrix& M,
               const IntMatrix &e,
               const int Z,
               const unsigned int k,
               const bool rootNotFixed,
               const bool forceDiploid);

    /// Get solution cost
    int getDelta();

    static HotStart firstCompleteHotStart(const InputInstance& inputInstance, const IntMatrix& e, const unsigned int k);

protected:
    /// TODO
    IloBoolVar3Array _bar_y;
    /// TODO
    IntMatrix _num_z;
    /// TODO
    IloBoolVar4Array _z;
    /// TODO
    IloIntVar4Array _a;
    /// TODO
    IloIntVar4Array _d;
    /// TODO
    IloIntVar4Array _bar_a;
    /// TODO
    IloIntVar4Array _bar_d;

    /// Build variables
    void buildVariables();
    /// Build constraints
    virtual void buildConstraints();
    /// Construct tree
    void constructTree();
};


#endif // _CARCHITECT_H_
