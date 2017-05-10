#ifndef _MARCHITECT_H_
#define _MARCHITECT_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <ilcplex/ilocplex.h>

#include "inputinstance.h"
#include "basic_types.h"

ILOSTLBEGIN

class MArchitect
{
public:
    MArchitect(const InputInstance& inputInstance,
               const Int3Array& C,
               const unsigned int k,
               const IntMatrix &e);

    virtual ~MArchitect()
    {
        _cplex.end();
        _model.end();
        _env.end();
    }

    void init();
    bool solve(const int timeLimit, const int memoryLimit);

    DoubleMatrix getM();

    double getObjValue()
    {
        return _cplex.getObjValue();
    }

    void exportModel(const std::string& filename) const;

    double getUB()
    {
        return _cplex.getObjValue();
    }

    double getLB()
    {
        return _cplex.getObjValue();
    }

    double getTime()
    {
        return _timer;
    }

protected:
    typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
    typedef IloArray<IloBoolVarMatrix> IloBoolVar3Array;
    typedef IloArray<IloBoolVar3Array> IloBoolVar4Array;
    typedef IloArray<IloBoolVar4Array> IloBoolVar5Array;

    typedef IloArray<IloNumVarArray> IloNumVarMatrix;
    typedef IloArray<IloNumVarMatrix> IloNumVar3Array;
    typedef IloArray<IloNumVar3Array> IloNumVar4Array;
    typedef IloArray<IloNumVar4Array> IloNumVar5Array;

    typedef IloArray<IloIntVarArray> IloIntVarMatrix;
    typedef IloArray<IloIntVarMatrix> IloIntVar3Array;
    typedef IloArray<IloIntVar3Array> IloIntVar4Array;
    typedef IloArray<IloIntVar4Array> IloIntVar5Array;

protected:
    /// Input instance
    const InputInstance& _inputInstance;
    /// Copy number matrix, just the leaves
    const Int3Array& _C;
    /// Number of leaves
    const unsigned int _k;
    /// Maximum copy number per chromosome, per position
    const IntMatrix& _e;
    /// Input frequencies from _inputInstance
    const Double3Array& _F;
    /// Number of chromosomes
    const unsigned int _numChr;
    /// Number of samples
    const unsigned int _m;
    /// Number of positions of each chromosome
    const IntArray& _n;
    /// CPLEX environment
    IloEnv _env;
    /// CPLEX model
    IloModel _model;
    /// CPLEX solver
    IloCplex _cplex;
    /// TODO
    IloNumVarMatrix _x;
    /// TODO
    IloNumVar3Array _bar_f;
    /// TODO
    IloExpr _obj;
    /// TODO
    double _timer;

    /// Build variables
    void buildVariables();
    /// Build constraints
    void buildConstraints();
    /// Build objective
    void buildObjective();
};

#endif // _MARCHITECT_H_
