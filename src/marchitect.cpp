#include "marchitect.h"
#include <lemon/time_measure.h>

MArchitect::MArchitect(const InputInstance& inputInstance,
                       const Int3Array& C,
                       const unsigned int k,
                       const IntMatrix &e)
    : _inputInstance(inputInstance)
    , _C(C)
    , _k(k)
    , _e(e)
    ,  _F(inputInstance.F())
    ,  _numChr(inputInstance.numChr())
    ,  _m(inputInstance.m())
    ,  _n(inputInstance.n())
    , _env()
    , _model(_env)
    , _cplex(_model)
    , _x()
    , _bar_f()
    , _obj()
    , _timer(0.0)
{
}

void MArchitect::init()
{
    buildVariables();
    buildConstraints();
    buildObjective();
}

bool MArchitect::solve(const int timeLimit, const int memoryLimit)
{
    if (g_verbosity >= VERBOSE_DEBUG)
    {
        _cplex.setOut(std::cerr);
        _cplex.setWarning(std::cerr);
        _cplex.setError(std::cerr);
        std::cerr << ">> Start CPLEX by M-step" << std::endl;
    }
    else
    {
        _cplex.setOut(_env.getNullStream());
        _cplex.setWarning(_env.getNullStream());
        _cplex.setError(_env.getNullStream());
    }

    if (timeLimit > 0)
    {
        _cplex.setParam(IloCplex::TiLim, timeLimit);
    }

    if (memoryLimit > 0)
    {
        _cplex.setParam(IloCplex::TreLim, memoryLimit);
    }
    
    _cplex.setParam(IloCplex::Threads, 1);

    lemon::Timer timer;
    bool res = _cplex.solve();
    _timer = timer.realTime();

    assert(res);
    if (res)
    {
        _cplex.out() << std::endl;
        _cplex.out() << "Solution status = " << _cplex.getStatus() << std::endl;
        _cplex.out() << "Solution LB = " << _cplex.getBestObjValue() << std::endl;
        _cplex.out() << "Solution UB = " << _cplex.getObjValue() << std::endl;
        _cplex.out() << "Runtime = " << _timer << " seconds" << std::endl;
    }
    else
    {
        _cplex.out() << "FAILED TO SOLVE:  " << _cplex.getStatus() << std::endl;
    }

    if (g_verbosity >= VERBOSE_DEBUG)
    {
        std::cerr << ">> End CPLEX" << std::endl;
    }

    return res;
}

DoubleMatrix MArchitect::getM()
{
    DoubleMatrix result(_m);
    for(unsigned int p = 0; p < _m; ++p)
    {
        result[p] = DoubleArray(_k);
        for(unsigned int i = 0; i < _k; ++i)
        {
            result[p][i] = _cplex.getValue(_x[p][i]);
        }
    }
    return result;
}

void MArchitect::buildVariables()
{
    _x = IloNumVarMatrix(_env, _m);
    for(unsigned int p = 0; p < _m; ++p)
    {
        _x[p] = IloNumVarArray(_env, _k);
        for(unsigned int i = 0; i < _k; ++i)
        {
            _x[p][i] = IloNumVar(_env, 0.0, 1.0);
        }
    }

    _bar_f = IloNumVar3Array(_env, _numChr);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        _bar_f[chr] = IloNumVarMatrix(_env, _m);
        for(unsigned int p = 0; p < _m; ++p)
        {
            _bar_f[chr][p] = IloNumVarArray(_env, _n[chr]);
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                _bar_f[chr][p][s] = IloNumVar(_env, 0, _e[chr][s]);
            }
        }
    }
}

void MArchitect::buildConstraints()
{
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int p = 0; p < _m; ++p)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                IloExpr sum_u_y(_env);

                for(unsigned int i = 0; i < _k; ++i)
                {
                    sum_u_y += _x[p][i] * _C[chr][i][s];
                }

                _model.add(_bar_f[chr][p][s] - _F[chr][p][s] + sum_u_y >= 0);
                _model.add(_bar_f[chr][p][s] + _F[chr][p][s] - sum_u_y >= 0);
            }
        }
    }

    for(unsigned int p = 0; p < _m; ++p)
    {
        IloExpr sum_x(_env);
        for(unsigned int i = 0; i < _k; ++i)
        {
            sum_x += _x[p][i];
        }
        _model.add(sum_x == 1.0);
    }
}

void MArchitect::buildObjective()
{
    _obj = IloExpr(_env);

    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int p = 0; p < _m; ++p)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                _obj += _bar_f[chr][p][s];
            }
        }
    }

    _model.add(IloMinimize(_env, _obj));
}
