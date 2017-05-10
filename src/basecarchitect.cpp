#include "basecarchitect.h"
#include <lemon/time_measure.h>


BaseCArchitect::BaseCArchitect(const InputInstance& inputInstance,
                               const DoubleMatrix& M,
                               const IntMatrix &e,
                               const int Z,
                               const unsigned int k,
                               const bool rootNotFixed,
                               const bool forceDiploid)
    : _inputInstance(inputInstance)
    , _M(M)
    , _e(e)
    , _Z(Z)
    , _k(k)
    , _forceDiploid(forceDiploid)
    , _rootNotFixed(rootNotFixed)
    , _F(inputInstance.F())
    , _numChr(inputInstance.numChr())
    , _m(inputInstance.m())
    , _n(inputInstance.n())
    , _num_vertices(2*_k - 1)
    , _env()
    , _model(_env)
    , _cplex(_model)
    , _T(_k, _numChr, _n)
    , _x()
    , _y()
    , _bar_f()
    , _obj(_env)
    , _completeHotStart(_env)
    , _partialHotStart(_env)
    , _timer(0.0)
{
}


void BaseCArchitect::init()
{
    buildVariables();
    buildConstraints();
    buildObjective();
}

bool BaseCArchitect::solve(const int timeLimit, const int memoryLimit, const int nrThreads)
{
    if (g_verbosity >= VERBOSE_DEBUG)
    {
        _cplex.setOut(std::cerr);
        _cplex.setWarning(std::cerr);
        _cplex.setError(std::cerr);
        std::cerr << ">> Start CPLEX by C-step" << std::endl;
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

    if (nrThreads > 0)
    {
        _cplex.setParam(IloCplex::Threads, nrThreads);
    }

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

        constructTree();
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

Int3Array BaseCArchitect::getC()
{
    Int3Array result(_numChr);
    for(unsigned chr = 0; chr < _numChr; ++chr)
    {
        result[chr] = IntMatrix(_k);
        for(unsigned int i = _k-1; i < _num_vertices; ++i)
        {
            unsigned int remapped_i = i - (_k-1);
            result[chr][remapped_i] = IntArray(_n[chr]);
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                result[chr][remapped_i][s] = _cplex.getIntValue(_y[chr][i][s]);
            }
        }
    }
    return result;
}

void BaseCArchitect::exportModel(const std::string& filename) const
{
    _cplex.exportModel(filename.c_str());
}

void BaseCArchitect::constructTree()
{
    // add and label arcs
    for(unsigned int i = 0; i < (_k - 1); ++i)
    {
        for(unsigned j = (i + 1); j < _num_vertices; ++j)
        {
            const int remapped_j = remap_j(i,j);
            if(_cplex.getIntValue(_x[i][remapped_j]))
            {
                _T.addArc(i, j);
            }
        }
    }

    // construct profiles
    for(unsigned int i = 0; i < _num_vertices; ++i)
    {
        CopyNumberTree::ProfileVector profile(_numChr);
        for (int chr = 0; chr < _numChr; ++chr)
        {
            profile[chr] = CopyNumberTree::Profile(_n[chr], 0);
            for (int s = 0; s < _n[chr]; ++s)
            {
                profile[chr][s] = _cplex.getIntValue(_y[chr][i][s]);
            }
        }
        _T.setProfile(i, profile);
    }
}

void BaseCArchitect::buildVariables()
{
    char buf[1024];

    _x = IloBoolVarMatrix(_env, _k - 1);
    for(unsigned int i = 0; i < (_k - 1); ++i)
    {
        _x[i] = IloBoolVarArray(_env, _num_vertices - (i + 1));
        for(unsigned j = (i + 1); j < _num_vertices; ++j)
        {
            snprintf(buf, 1024, "x_%d_%d", i, j);
            _x[i][remap_j(i, j)] = IloBoolVar(_env, buf);
            _completeHotStart.add(_x[i][remap_j(i, j)]);
            _partialHotStart.add(_x[i][remap_j(i, j)]);
        }
    }

    _y = IloIntVar3Array(_env, _numChr);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        _y[chr] = IloIntVarMatrix(_env, _num_vertices);
        for(unsigned int i = 0; i < _num_vertices; ++i)
        {
            _y[chr][i] = IloIntVarArray(_env, _n[chr]);
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                snprintf(buf, 1024, "y_%d_%d_%d", chr, i, s);
                _y[chr][i][s] = IloIntVar(_env, 0, _e[chr][s], buf);
                _completeHotStart.add(_y[chr][i][s]);
                _partialHotStart.add(_y[chr][i][s]);
            }
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
                snprintf(buf, 1024, "bar_f_%d_%d_%d", chr, p, s);
                _bar_f[chr][p][s] = IloNumVar(_env, 0, _e[chr][s], buf);
            }
        }
    }
}

void BaseCArchitect::buildConstraints()
{
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int p = 0; p < _m; ++p)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                IloExpr sum_u_y(_env);

                for(unsigned int i = (_k - 1); i < _num_vertices; ++i)
                {
                    int remapped_i = i - (_k - 1);
                    sum_u_y += _M[p][remapped_i] * _y[chr][i][s];
                }

                _model.add(_bar_f[chr][p][s] - _F[chr][p][s] + sum_u_y >= 0);
                _model.add(_bar_f[chr][p][s] + _F[chr][p][s] - sum_u_y >= 0);
            }
        }
    }

    char buf[1024];
    /**
     \sum_{i \in \delta^-(j)} x_{i,j} = 1 for 1 < j <= 2k-1
     **/
    for(unsigned int j = 1; j < _num_vertices; ++j)
    {
        IloExpr sum_x_minus(_env);
        unsigned int lim = std::min(j, _k - 1);

        for(unsigned i = 0; i < lim; ++i)
        {
            sum_x_minus += _x[i][remap_j(i, j)];
        }

        IloConstraint cons(sum_x_minus == 1);
        snprintf(buf, 1024, "in_deg_%d", j);
        cons.setName(buf);
        _model.add(cons);
    }

    /**
     \sum_{j \in \delta^+(i)} x_{i,j} = 2 for 1 <= i < k
     **/
    for(unsigned int i = 0; i < (_k - 1); ++i)
    {
        IloExpr sum_x_plus(_env);

        for(unsigned int j = (i + 1); j < _num_vertices; ++j)
        {
            sum_x_plus += _x[i][remap_j(i,j)];
        }

        IloConstraint cons(sum_x_plus == 2);
        snprintf(buf, 1024, "out_deg_%d", i);
        cons.setName(buf);
        _model.add(cons);
    }

    /**
     y_{1,s} = 2 for 1 <= s <= n
     **/
    if (!_rootNotFixed)
    {
        for(unsigned int chr = 0; chr < _numChr; ++chr)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                IloConstraint cons(_y[chr][0][s] == 2);
                snprintf(buf, 1024, "root_%d", s);
                cons.setName(buf);
                _model.add(cons);
            }
        }
    }

    if(_forceDiploid)
    {
        for(unsigned int chr = 0; chr < _numChr; ++chr)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                IloConstraint cons(_y[chr][_k-1][s] == 2);
                snprintf(buf, 1024, "force_diploid_%d", s);
                cons.setName(buf);
                _model.add(cons);
            }
        }

        _model.add(_x[0][remap_j(0, 1)] == 1);
        _model.add(_x[0][remap_j(0, _k-1)] == 1);
    }
}

void BaseCArchitect::buildObjective()
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

HotStart BaseCArchitect::getCompleteHotStart()
{
    int size =  _completeHotStart.getSize();
    IloNumArray values(_env, size);
    _cplex.getValues(_completeHotStart, values);
    HotStart result(size, 0);
    for(unsigned int i = 0; i < size; ++i)
    {
        result[i] = std::floor(values[i] + 0.5);
    }
    values.end();
    return result;
}

HotStart BaseCArchitect::getPartialHotStart()
{
    int size =  _partialHotStart.getSize();
    IloNumArray values(_env, size);
    _cplex.getValues(_partialHotStart, values);
    HotStart result(size, 0);
    for(unsigned int i = 0; i < size; ++i)
    {
        result[i] = std::floor(values[i] + 0.5);
    }
    values.end();
    return result;
}

int BaseCArchitect::addCompleteHotStart(const HotStart& completeHotStart)
{
    int size =  _completeHotStart.getSize();
    if(completeHotStart.size() != size)
    {
        throw std::runtime_error("ERROR: Complete Hot Start with an inconsistent number of elements");
    }
    IloIntArray cast(_env, size);
    for(unsigned int i = 0; i < size; ++i)
    {
        cast[i] = completeHotStart[i];
    }
    int res = _cplex.addMIPStart(_completeHotStart, cast.toNumArray(),
                                 IloCplex::MIPStartEffort::MIPStartSolveFixed, "Complete Hot Start");
    cast.end();
    return res;
}

int BaseCArchitect::addPartialHotStart(const HotStart& partialHotStart)
{
    int size =  _partialHotStart.getSize();
    if(partialHotStart.size() != size)
    {
        throw std::runtime_error("ERROR: Partial Hot Start with an inconsistent number of elements");
    }
    IloIntArray cast(_env, size);
    for(unsigned int i = 0; i < size; ++i)
    {
        cast[i] = partialHotStart[i];
    }
    int res = _cplex.addMIPStart(_partialHotStart, cast.toNumArray(),
                                 IloCplex::MIPStartEffort::MIPStartSolveMIP, "Partial Hot Start");
    cast.end();
    return res;
}


