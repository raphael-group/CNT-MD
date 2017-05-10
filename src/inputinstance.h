#ifndef _INPUTINSTANCE_H
#define _INPUTINSTANCE_H

#include "basic_types.h"
#include <iostream>
#include <math.h>
#include <algorithm>

class InputInstance
{
public:
    InputInstance();
    
    friend std::ostream& operator<<(std::ostream& out, const InputInstance& instance);
    friend std::istream& operator>>(std::istream& in, InputInstance& instance);
    
    const Double3Array& F() const
    {
        return _F;
    }
    
    int m() const
    {
        return _m;
    }
    
    int numChr() const
    {
        return _num_chr;
    }
    
    const IntArray& n() const
    {
        return _n;
    }
    
    /// Extract maximum copy number from F
    int e() const;
    
    IntArray z() const;
    
private:
    Double3Array _F;
    int _m;
    int _num_chr;
    IntArray _n;
};

std::ostream& operator<<(std::ostream& out, const InputInstance& instance);
std::istream& operator>>(std::istream& in, InputInstance& instance);

#endif // _INPUTINSTANCE_H
