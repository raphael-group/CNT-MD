#include "refiner.h"


Refiner::Refiner(const InputInstance& inputInstance,
                 const DoubleMatrix& M,
                 const IntMatrix &e,
                 const int Z,
                 const unsigned int k,
                 const bool rootNotFixed,
                 const bool forceDiploid,
                 const double bound_bar_f)
    : CArchitect(inputInstance, M, e, Z, k, rootNotFixed, forceDiploid)
    , _bound_bar_f(bound_bar_f)
{
}


double Refiner::getDistance()
{
    double result = 0;
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int p = 0; p < _m; ++p)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                result += _cplex.getValue(_bar_f[chr][p][s]);
            }
        }
    }
    return result;
}


void Refiner::buildConstraints()
{
    CArchitect::buildConstraints();

    IloExpr sum_bar_f(_env);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int p = 0; p < _m; ++p)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                sum_bar_f += _bar_f[chr][p][s];
            }
        }
    }
     _model.add(sum_bar_f <= _bound_bar_f);
}


void Refiner::buildObjective()
{
    _obj = IloExpr(_env);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    _obj += _bar_a[chr][i][remapped_j][l] + _bar_d[chr][i][remapped_j][l];
                }
            }
        }
    }
    _model.add(IloMinimize(_env, _obj));
}
