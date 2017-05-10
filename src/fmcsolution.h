#ifndef _FMCSOLUTION_H
#define _FMCSOLUTION_H

#include <iostream>
#include <math.h>
#include <algorithm>

#include "basic_types.h"
#include "copynumbertree.h"
#include "inputinstance.h"

class FMCSolution
{
public:
    FMCSolution();
    
    FMCSolution(const CopyNumberTree& tree, const DoubleMatrix& M, const InputInstance& input);
    
    friend std::ostream& operator<<(std::ostream& out, const FMCSolution& instance);
    friend std::istream& operator>>(std::istream& in, FMCSolution& instance);

    const CopyNumberTree& getTree() const
    {
        return _tree;
    }

    const DoubleMatrix& getM() const
    {
        return _M;
    }

    const InputInstance& getInput() const
    {
        return _input;
    }

    void writeDOT(std::ostream& out) const;
    
private:
    /// Copy number tree
    CopyNumberTree _tree;
    /// Usage matrix
    DoubleMatrix _M;
    /// Input instance
    InputInstance _input;
};

std::ostream& operator<<(std::ostream& out, const FMCSolution& instance);
std::istream& operator>>(std::istream& in, FMCSolution& instance);

#endif // _FMCSOLUTION_H
