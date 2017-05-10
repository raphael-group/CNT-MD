#include "carchitect.h"


CArchitect::CArchitect(const InputInstance& inputInstance,
                       const DoubleMatrix& M,
                       const IntMatrix &e,
                       const int Z,
                       const unsigned int k,
                       const bool rootNotFixed,
                       const bool forceDiploid)
    : BaseCArchitect(inputInstance, M, e, Z, k, rootNotFixed, forceDiploid)
    , _bar_y()
    , _num_z()
    , _z()
    , _a()
    , _d()
    , _bar_a()
    , _bar_d()
{
}

int CArchitect::getDelta()
{
    double result = 0.0;
    
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);

                if(_cplex.getIntValue(_x[i][remapped_j]) > 0)
                {
                    for(unsigned int l = 0; l < _n[chr]; ++l)
                    {
                        result += _cplex.getIntValue(_bar_a[chr][i][remapped_j][l]);
                        result += _cplex.getIntValue(_bar_d[chr][i][remapped_j][l]);
                    }
                }
            }
        }
    }
    return result;
}

void CArchitect::constructTree()
{
    // construct tree topology
    BaseCArchitect::constructTree();

    //collect the number of amplifications and deletions per segment
    Int4Array amplifications(_numChr);
    Int4Array deletions(_numChr);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        amplifications[chr] = Int3Array(_k - 1);
        deletions[chr] = Int3Array(_k - 1);
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            amplifications[chr][i] = IntMatrix(_num_vertices - (i + 1));
            deletions[chr][i] = IntMatrix(_num_vertices - (i + 1));
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                amplifications[chr][i][remapped_j] = IntArray(_n[chr]);
                deletions[chr][i][remapped_j] = IntArray(_n[chr]);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    amplifications[chr][i][remapped_j][l] = _cplex.getIntValue(_a[chr][i][remapped_j][l]);
                    deletions[chr][i][remapped_j][l] = _cplex.getIntValue(_d[chr][i][remapped_j][l]);
                }
            }
        }
    }

    // construct events
    const CopyNumberTree::Digraph& T = _T.T();
    for (CopyNumberTree::ArcIt a_ij(T); a_ij != lemon::INVALID; ++a_ij)
    {
        CopyNumberTree::Node v_i = T.source(a_ij);
        CopyNumberTree::Node v_j = T.target(a_ij);
        const int i = _T.index(v_i);
        const int j = _T.index(v_j);
        const int remapped_j = remap_j(i, j);
      
        for (unsigned int chr = 0; chr < _numChr; ++chr)
        {
            unsigned int counter = 0;
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                while(deletions[chr][i][remapped_j][s] != 0)
                {
                    int t = s;
                    while (t+1 < _n[chr] && deletions[chr][i][remapped_j][t+1] != 0) ++t;

                    _T.addEvent(chr, i, j, CopyNumberTree::Event(chr, s, t, -1));
                    ++counter;

                    for (int l = s; l <= t; ++l)
                        --deletions[chr][i][remapped_j][l];
                }
            }

            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                while(amplifications[chr][i][remapped_j][s] != 0)
                {
                    int t = s;
                    while (t+1 < _n[chr] && amplifications[chr][i][remapped_j][t+1] != 0) ++t;

                    _T.addEvent(chr, i, j, CopyNumberTree::Event(chr, s, t, 1));
                    ++counter;

                    for (int l = s; l <= t; ++l)
                        --amplifications[chr][i][remapped_j][l];
                }
            }

            unsigned int w = 0;
            for(unsigned int l = 0; l < _n[chr]; ++l)
            {
                w += _cplex.getIntValue(_bar_a[chr][i][remapped_j][l]) + _cplex.getIntValue(_bar_d[chr][i][remapped_j][l]);
            }

            // number of events is at most w, because we don't minimize for it
            assert(counter <= w);
        }
    }
}

void CArchitect::buildVariables()
{
    BaseCArchitect::buildVariables();

    char buf[1024];

    _bar_y = IloBoolVar3Array(_env, _numChr);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        _bar_y[chr] = IloBoolVarMatrix(_env, _num_vertices);
        for(unsigned int i = 0; i < _num_vertices; ++i)
        {
            _bar_y[chr][i] = IloBoolVarArray(_env, _n[chr]);
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                snprintf(buf, 1024, "bar_y_%d_%d_%d", chr, i, s);
                _bar_y[chr][i][s] = IloBoolVar(_env, buf);
                _completeHotStart.add(_bar_y[chr][i][s]);
            }
        }
    }

    _num_z = IntMatrix(_numChr);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        _num_z[chr] = IntArray(_n[chr]);
        for(unsigned int s = 0; s < _n[chr]; ++s)
        {
            _num_z[chr][s] = floor(log(_e[chr][s]) / log(2.0)) + 1;
        }
    }

    _z = IloBoolVar4Array(_env, _numChr);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        _z[chr] = IloBoolVar3Array(_env, _num_vertices);
        for(unsigned int i = 0; i < _num_vertices; ++i)
        {
            _z[chr][i] = IloBoolVarMatrix(_env, _n[chr]);
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                _z[chr][i][s] = IloBoolVarArray(_env, _num_z[chr][s]);
                for(unsigned int q = 0; q < _num_z[chr][s]; ++q)
                {
                    snprintf(buf, 1024, "z_%d_%d_%d_%d", chr, i, s, q);
                    _z[chr][i][s][q] = IloBoolVar(_env, buf);
                    _completeHotStart.add(_z[chr][i][s][q]);
                }
            }
        }
    }

    _a = IloIntVar4Array(_env, _numChr);
    _d = IloIntVar4Array(_env, _numChr);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        const int e = *max_element(_e[chr].begin(), _e[chr].end());
        _a[chr] = IloIntVar3Array(_env, _k - 1);
        _d[chr] = IloIntVar3Array(_env, _k - 1);
        for(unsigned int i = 0; i < _k - 1; ++i)
        {
            _a[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            _d[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            for(unsigned int j = i + 1; j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                _a[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                _d[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                for(unsigned int s = 0; s < _n[chr]; ++s)
                {
                    _a[chr][i][remapped_j][s] = IloIntVar(_env, 0, e);
                    _completeHotStart.add(_a[chr][i][remapped_j][s]);

                    _d[chr][i][remapped_j][s] = IloIntVar(_env, 0, e);
                    _completeHotStart.add(_d[chr][i][remapped_j][s]);
                }
            }
        }
    }

    _bar_a = IloIntVar4Array(_env, _numChr);
    _bar_d = IloIntVar4Array(_env, _numChr);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        _bar_a[chr] = IloIntVar3Array(_env, _k - 1);
        _bar_d[chr] = IloIntVar3Array(_env, _k - 1);
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            _bar_a[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            _bar_d[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                _bar_a[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                _bar_d[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    _bar_a[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr][l]);
                    _completeHotStart.add(_bar_a[chr][i][remapped_j][l]);

                    _bar_d[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr][l]);
                    _completeHotStart.add(_bar_d[chr][i][remapped_j][l]);
                }
            }
        }
    }
}

void CArchitect::buildConstraints()
{
    BaseCArchitect::buildConstraints();

    char buf[1024];

    IloExpr sum_bars(_env);
    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    sum_bars += _bar_a[chr][i][remapped_j][l] + _bar_d[chr][i][remapped_j][l];
                }
            }
        }
    }
    _model.add(sum_bars <= _Z);

    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        for(unsigned int i = 0; i < _num_vertices; ++i)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                IloExpr sum_z_2(_env);
                IloExpr sum_z(_env);

                for(unsigned int q = 0; q < _num_z[chr][s]; ++q)
                {
                    sum_z_2 += (1 << q) * _z[chr][i][s][q];
                    sum_z += _z[chr][i][s][q];

                    IloConstraint cons(_bar_y[chr][i][s] - _z[chr][i][s][q] >= 0);
                    snprintf(buf, 1024, "nonzero_lb_%d_%d_%d_%d", chr, i, s, q);
                    cons.setName(buf);
                    _model.add(cons);
                }

                IloConstraint cons2(_y[chr][i][s] - sum_z_2 == 0);
                snprintf(buf, 1024, "binary_%d_%d_%d", chr, i, s);
                cons2.setName(buf);
                _model.add(cons2);

                IloConstraint cons3(_bar_y[chr][i][s] - sum_z <= 0);
                snprintf(buf, 1024, "nonzero_ub_%d_%d_%d", chr, i, s);
                cons3.setName(buf);
                _model.add(cons3);
            }
        }
    }

    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        const int e = *max_element(_e[chr].begin(), _e[chr].end());
        for(unsigned int i = 0; i < _k - 1; ++i)
        {
            for(unsigned int j = i + 1; j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int s = 0; s < _n[chr]; ++s)
                {
                    _model.add(_bar_y[chr][i][s]
                               - _bar_y[chr][j][s]
                               + 1
                               - _x[i][remapped_j] >= 0);
                    _model.add(_y[chr][i][s]
                               - _d[chr][i][remapped_j][s]
                               - e
                               + e * _bar_y[chr][i][s]
                               - e * _bar_y[chr][j][s]
                               - e
                               + e * _x[i][remapped_j] <= 0);
                    _model.add(_y[chr][j][s]
                               - _y[chr][i][s]
                               + _d[chr][i][remapped_j][s]
                               - _a[chr][i][remapped_j][s]
                               - 4 * e
                               + (2 * e) * _bar_y[chr][i][s]
                               + (2 * e) * _bar_y[chr][j][s]
                               - 2 * e
                               + (2 * e) * _x[i][remapped_j] <= 0);
                    _model.add(_y[chr][j][s]
                               - _y[chr][i][s]
                               + _d[chr][i][remapped_j][s]
                               - _a[chr][i][remapped_j][s]
                               + 4 * e
                               - (2 * e) * _bar_y[chr][i][s]
                               - (2 * e) * _bar_y[chr][j][s]
                               + 2 * e
                               - (2 * e) * _x[i][remapped_j] >= 0);
                    _model.add(_d[chr][i][remapped_j][s]
                               - _y[chr][i][s]
                               + 1
                               - (e + 1) * 2
                               + (e + 1) * _bar_y[chr][i][s]
                               + (e + 1) * _bar_y[chr][j][s] <= 0);
                }
            }
        }
    }

    for(unsigned int chr = 0; chr < _numChr; ++chr)
    {
        const int e = *max_element(_e[chr].begin(), _e[chr].end());
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    if(l == 0)
                    {
                        _model.add(_bar_a[chr][i][remapped_j][l] - _a[chr][i][remapped_j][l] >= 0);
                        _model.add(_bar_d[chr][i][remapped_j][l] - _d[chr][i][remapped_j][l] >= 0);
                    }
                    else
                    {
                        _model.add(_bar_a[chr][i][remapped_j][l] - _a[chr][i][remapped_j][l] + _a[chr][i][remapped_j][l-1] >= 0);
                        _model.add(_bar_d[chr][i][remapped_j][l] - _d[chr][i][remapped_j][l] + _d[chr][i][remapped_j][l-1] >= 0);
                    }

                    _model.add(_a[chr][i][remapped_j][l]
                               - e * _x[i][remapped_j] <= 0);
                    _model.add(_d[chr][i][remapped_j][l]
                               - e * _x[i][remapped_j] <= 0);
                    _model.add(_bar_a[chr][i][remapped_j][l]
                               - _e[chr][l] * _x[i][remapped_j] <= 0);
                    _model.add(_bar_d[chr][i][remapped_j][l]
                               - _e[chr][l] * _x[i][remapped_j] <= 0);
                }
            }
        }
    }
}

HotStart CArchitect::firstCompleteHotStart(const InputInstance& inputInstance,
                                           const IntMatrix& e,
                                           const unsigned int k)
{
    DoubleMatrix M(inputInstance.m(), DoubleArray(k, 0.0));
    for(unsigned int p = 0; p < inputInstance.m(); ++p)
    {
        M[p][0] = 1.0;
    }

    const int Z = 0;

    g_mutex.lock();
    CArchitect architect(inputInstance, M, e, Z, k, false, true);
    g_mutex.unlock();
    
    architect.init();
    architect.solve(0, 0, 1);

    return architect.getCompleteHotStart();
}
