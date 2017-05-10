#include "tripletarchitect.h"
#include <lemon/time_measure.h>

TripletArchitect::TripletArchitect(const IntMatrix& y1, const IntMatrix& y2)
    : _y1(y1)
    , _y2(y2)
    , _e(y1.size())
    , _num_chr(y1.size())
    , _k(2)
    , _num_vertices(3)
    , _n(_num_chr)
    , _env()
    , _model(_env)
    , _cplex(_model)
    , _T()
    , _x()
    , _y()
    , _bar_y()
    , _num_z()
    , _z()
    , _obj()
    , _a()
    , _d()
    , _A()
    , _D()
{
    // initialize _n
    for (int chr = 0; chr < _num_chr; ++chr)
    {
        _n[chr] = y1[chr].size();
    }
    
    // initialize _T
    _T = CopyNumberTree(_k, _num_chr, _n);
    
    // initialize e
    for (int chr = 0; chr < _num_chr; ++chr)
    {
        int tmp1 = *std::max_element(y1[chr].begin(), y1[chr].end());
        int tmp2 = *std::max_element(y2[chr].begin(), y2[chr].end());
        
        _e[chr] = std::max(tmp1, tmp2);
    }
}

int TripletArchitect::cost(const IntMatrix &y1, const IntMatrix &y2)
{
    TripletArchitect triplet(y1, y2);
    triplet.init();
    triplet.solve(-1, -1);
    return triplet.getObjValue();
}

void TripletArchitect::init()
{
    buildVariables();
    buildConstraints();
    buildObjective();
}

bool TripletArchitect::solve(const int timeLimit, const int memoryLimit)
{
    _cplex.setOut(std::cerr);
    _cplex.setWarning(std::cerr);
    _cplex.setError(std::cerr);
    
    // shut up cplex
    bool verbose = false;
    if (!verbose) {
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
    
    lemon::Timer timer;
    bool res = _cplex.solve();
    if (res)
    {
        _cplex.out() << std::endl;
        _cplex.out() << "Solution status = " << _cplex.getStatus() << std::endl;
        _cplex.out() << "Solution LB = " << _cplex.getBestObjValue() << std::endl;
        _cplex.out() << "Solution UB = " << _cplex.getObjValue() << std::endl;
        _cplex.out() << "Runtime = " << timer.realTime() << " seconds" << std::endl;
        
        constructTree();
    }
    else
    {
        _cplex.out() << "FAILED TO SOLVE:  " << _cplex.getStatus() << std::endl;
    }
    
    return res;
}

void TripletArchitect::exportModel(const std::string& filename) const
{
    _cplex.exportModel(filename.c_str());
}

void TripletArchitect::constructTree()
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
        CopyNumberTree::ProfileVector profile(_num_chr);
        for (int chr = 0; chr < _num_chr; ++chr)
        {
            profile[chr] = CopyNumberTree::Profile(_n[chr], 0);
            for (int s = 0; s < _n[chr]; ++s)
            {
                profile[chr][s] = _cplex.getIntValue(_y[chr][i][s]);
            }
        }
        _T.setProfile(i, profile);
    }
    
    //collect the number of amplifications and deletions per segment
    Int4Array amplifications(_num_chr);
    Int4Array deletions(_num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
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
        
        for (unsigned int chr = 0; chr < _num_chr; ++chr)
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
                w += _cplex.getIntValue(_A[chr][i][remapped_j][l]) + _cplex.getIntValue(_D[chr][i][remapped_j][l]);
            }
            
            assert(counter == w);
        }
    }
}

void TripletArchitect::buildVariables()
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
        }
    }
    
    _y = IloIntVar3Array(_env, _num_chr);
    _bar_y = IloBoolVar3Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _y[chr] = IloIntVarMatrix(_env, _num_vertices);
        _bar_y[chr] = IloBoolVarMatrix(_env, _num_vertices);
        for(unsigned int i = 0; i < _num_vertices; ++i)
        {
            _y[chr][i] = IloIntVarArray(_env, _n[chr]);
            _bar_y[chr][i] = IloBoolVarArray(_env, _n[chr]);
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                snprintf(buf, 1024, "y_%d_%d_%d", chr, i, s);
                _y[chr][i][s] = IloIntVar(_env, 0, _e[chr], buf);
                
                snprintf(buf, 1024, "bar_y_%d_%d_%d", chr, i, s);
                _bar_y[chr][i][s] = IloBoolVar(_env, buf);
            }
        }
    }
    
    _num_z = IntMatrix(_num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _num_z[chr] = IntArray(_n[chr]);
        for(unsigned int s = 0; s < _n[chr]; ++s)
        {
            _num_z[chr][s] = floor(log(_e[chr]) / log(2.0)) + 1;
        }
    }
    
    _z = IloBoolVar4Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
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
                }
            }
        }
    }
    
    _a = IloIntVar4Array(_env, _num_chr);
    _d = IloIntVar4Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        const int e = _e[chr];
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
                    _d[chr][i][remapped_j][s] = IloIntVar(_env, 0, e);
                }
            }
        }
    }
    
    _A = IloIntVar4Array(_env, _num_chr);
    _D = IloIntVar4Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _A[chr] = IloIntVar3Array(_env, _k - 1);
        _D[chr] = IloIntVar3Array(_env, _k - 1);
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            _A[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            _D[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                _A[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                _D[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    _A[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr]);
                    _D[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr]);
                }
            }
        }
    }
}

void TripletArchitect::buildConstraints()
{
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
     y_{i,s} = c_{i-k+2,s} for  k <= i <= 2k-1 and 1 <= s <= n
     **/
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = (_k - 1); i < _num_vertices; ++i)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                if (i == 1)
                {
                    IloConstraint cons(_y[chr][i][s] == _y1[chr][s]);
                    snprintf(buf, 1024, "leaf_%d_%d", i, s);
                    cons.setName(buf);
                    _model.add(cons);
                }
                else if (i == 2)
                {
                    IloConstraint cons(_y[chr][i][s] == _y2[chr][s]);
                    snprintf(buf, 1024, "leaf_%d_%d", i, s);
                    cons.setName(buf);
                    _model.add(cons);
                }
                else
                {
                    assert(false);
                }
            }
        }
    }
    
    /**
     1. y_{i,s} = \sum_{j=0}^{\lceil \log_2(e) \rceil} 2^j \cdot z_{i,s,j} for 1 <= i <= 2k-1 and 1 <= s <= n
     2. \bar{y}_{i,s}  <= \sum_{j=0}^{\lceil \log_2(e) \rceil} z_{i,s,j} for 1 <= i <= 2k-1 and 1 <= s <= n
     3. \bar{y}_{i,s}  \ge z_{i,s,j} for  1 <= i <= 2k-1 and 1 <= s <= n and 0 <= j <= (\lceil \log_2(e) \rceil)
     **/
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
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
    
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    if(l == 0)
                    {
                        _model.add(_A[chr][i][remapped_j][l] - _a[chr][i][remapped_j][l] >= 0);
                        _model.add(_D[chr][i][remapped_j][l] - _d[chr][i][remapped_j][l] >= 0);
                    }
                    else
                    {
                        _model.add(_A[chr][i][remapped_j][l] - _a[chr][i][remapped_j][l] + _a[chr][i][remapped_j][l-1] >= 0);
                        _model.add(_D[chr][i][remapped_j][l] - _d[chr][i][remapped_j][l] + _d[chr][i][remapped_j][l-1] >= 0);
                    }
                }
            }
        }
    }
    
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        const int e = _e[chr];
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
}

void TripletArchitect::buildObjective()
{
    _obj = IloExpr(_env);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    _obj += _A[chr][i][remapped_j][l] + _D[chr][i][remapped_j][l];
                }
            }
        }
    }
    _model.add(IloMinimize(_env, _obj));
}

