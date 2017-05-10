#ifndef _BASECARCHITECT_H_
#define _BASECARCHITECT_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <ilcplex/ilocplex.h>
#include "copynumbertree.h"
#include "inputinstance.h"

#include "basic_types.h"

ILOSTLBEGIN

class BaseCArchitect
{
public:
    BaseCArchitect(const InputInstance& inputInstance,
                   const DoubleMatrix& M,
                   const IntMatrix &e,
                   const int Z,
                   const unsigned int k,
                   const bool rootNotFixed,
                   const bool forceDiploid);

    virtual ~BaseCArchitect()
    {
        _completeHotStart.end();
        _partialHotStart.end();
        _cplex.end();
        _model.end();
        _env.end();
    }

    void init();
    
    bool solve(const int timeLimit, const int memoryLimit, const int nrThreads);

    Int3Array getC();
    
    virtual int getDelta() = 0;

    double getObjValue()
    {
        return _cplex.getObjValue();
    }

    CopyNumberTree& getTree()
    {
        return _T;
    }

    void exportModel(const std::string& filename) const;

    int addCompleteHotStart(const HotStart& completeHotStart);
    int addPartialHotStart(const HotStart& partialHotStart);

    HotStart getCompleteHotStart();
    HotStart getPartialHotStart();

    double getUB()
    {
        return _cplex.getObjValue();
    }

    double getLB()
    {
        return _cplex.getBestObjValue();
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

    int remap_j(int i, int j) const
    {
        assert(0 <= i && i < _k - 1);
        return j - (i + 1);
    }

    /// Build variables
    virtual void buildVariables();
    /// Build constraints
    virtual void buildConstraints();
    /// Build objective
    virtual void buildObjective();
    /// Construct tree
    virtual void constructTree();
    /// Return a complete Hot Start


protected:
    /// Input instance
    const InputInstance& _inputInstance;
    /// Usage matrix
    const DoubleMatrix& _M;
    /// Maximum copy number per chromosome, per position
    const IntMatrix &_e;
    /// Maximum number of events per chromosome
    const int _Z;
    /// Number of leaves
    const unsigned int _k;
    /// Force the presence of the normal diploid clone
    const bool _forceDiploid;
    /// Do not fix the root to the normal diploid
    const bool _rootNotFixed;
    /// Input frequencies from _inputInstance
    const Double3Array& _F;
    /// Number of chromosomes
    const unsigned int _numChr;
    /// Number of samples
    const unsigned int _m;
    /// Number of positions of each chromosome
    const IntArray& _n;
    /// Number of leaves + inner nodes: 2k - 1
    const unsigned int _num_vertices;
    /// CPLEX environment
    IloEnv _env;
    /// CPLEX model
    IloModel _model;
    /// CPLEX solver
    IloCplex _cplex;
    /// Resulting copy number tree
    CopyNumberTree _T;
    /// TODO
    IloBoolVarMatrix _x;
    /// TODO
    IloIntVar3Array _y;
    /// TODO
    IloNumVar3Array _bar_f;
    /// TODO
    IloExpr _obj;
    /// TODO
    IloNumVarArray _completeHotStart;
    /// TODO
    IloNumVarArray _partialHotStart;
    /// TODO
    double _timer;
};

#endif // _BASECARCHITECT_H_
