#ifndef _TRIPLETARCHITECT_H_
#define _TRIPLETARCHITECT_H_

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
#include "basic_types.h"
#include "copynumbertree.h"

class TripletArchitect
{
public:
    TripletArchitect(const IntMatrix& y1, const IntMatrix& y2);
    
    void init();
    
    bool solve(const int timeLimit, const int memoryLimit);
    
    double getObjValue()
    {
        return _cplex.getObjValue();
    }
    
    const CopyNumberTree& getTree()
    {
        return _T;
    }
    
    static int cost(const IntMatrix& y1, const IntMatrix& y2);
    
    void exportModel(const std::string& filename) const;
    
private:
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
    
    void constructTree();
    
private:
    const IntMatrix& _y1;
    const IntMatrix& _y2;
    IntArray _e;
    const unsigned int _num_chr;
    const unsigned int _k;
    const unsigned int _num_vertices;
    IntArray _n;
    
    IloEnv _env;
    IloModel _model;
    IloCplex _cplex;
    CopyNumberTree _T;
    
    /// _x[i][j] == 1 if and only if there is an edge (i,j)
    IloBoolVarMatrix _x;
    IloIntVar3Array _y;
    IloBoolVar3Array _bar_y;
    IntMatrix _num_z;
    IloBoolVar4Array _z;
    
    IloExpr _obj;
    
    IloIntVar4Array _a;
    IloIntVar4Array _d;
    IloIntVar4Array _A;
    IloIntVar4Array _D;
    
    void buildVariables();
    void buildConstraints();
    void buildObjective();
};

#endif // _TRIPLETARCHITECT_H_
