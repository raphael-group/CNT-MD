#include "comparison.h"
#include "tripletarchitect.h"

Comparison::Comparison(const FMCSolution& T1, const FMCSolution& T2)
    : _T1(T1)
    , _T2(T2)
    , _k1(_T1.getTree().k())
    , _k2(_T2.getTree().k())
    , _G()
    , _dummy(_G)
    , _redNodeToIdx(_G)
    , _blueNodeToIdx(_G)
    , _idxToRedNode(_k1, lemon::INVALID)
    , _idxToBlueNode(_k2, lemon::INVALID)
    , _weight(_G)
    , _mwm(_G, _weight)
{
}

bool Comparison::init()
{
    // check #chromosomes
    if (_T1.getTree().numChr() != _T2.getTree().numChr())
    {
        return false;
    }
    
    for (int i1 = 0; i1 < _k1; ++i1)
    {
        BpRedNode v = _G.addRedNode();
        _redNodeToIdx[v] = i1;
        _idxToRedNode[i1] = v;
        _dummy[v] = false;
    }
    
    for (int i2 = 0; i2 < _k2; ++i2)
    {
        BpBlueNode v = _G.addBlueNode();
        _blueNodeToIdx[v] = i2;
        _idxToBlueNode[i2] = v;
        _dummy[v] = false;
    }
    
    for (int i1 = 0; i1 < _k1; ++i1)
    {
        for (int i2 = 0; i2 < _k2; ++i2)
        {
            const CopyNumberTree::ProfileVector& y1 = _T1.getTree().profile(i1 + _k1 - 1);
            const CopyNumberTree::ProfileVector& y2 = _T2.getTree().profile(i2 + _k2 - 1);
            
//            std::cout << "i1 = " << i1 << ", y1 = " << CopyNumberTree::profileToString(y1) << std::endl;
//            std::cout << "i2 = " << i2 << ", y2 = " << CopyNumberTree::profileToString(y2) << std::endl;
            
            BpEdge edge = _G.addEdge(_idxToRedNode[i1], _idxToBlueNode[i2]);
            _weight[edge] = -TripletArchitect::cost(y1, y2);
            
//            std::cout << "cost: " << _weight[edge] << std::endl << std::endl;
        }
    }
    
    // add dummy nodes
    if (_k1 < _k2)
    {
        int count = _k2 - _k1;
        for (int i = 0; i < count; ++i)
        {
            BpRedNode red = _G.addRedNode();
            _dummy[red] = true;
            
            for (int i2 = 0; i2 < _k2; ++i2)
            {
                BpEdge edge = _G.addEdge(red, _idxToBlueNode[i2]);
                _weight[edge] = 0;
            }
        }
    }
    else if (_k2 < _k1)
    {
        int count = _k1 - _k2;
        for (int i = 0; i < count; ++i)
        {
            BpBlueNode blue = _G.addBlueNode();
            _dummy[blue] = true;
            
            for (int i1 = 0; i1 < _k1; ++i1)
            {
                BpEdge edge = _G.addEdge(_idxToRedNode[i1], blue);
                _weight[edge] = 0;
            }
        }
    }
    
    _mwm.init();
    _mwm.start();
    
    _leaf1ToLeaf2 = std::vector<int>(2*_k1 - 1, -1);
    _leaf2ToLeaf1 = std::vector<int>(2*_k2 - 1, -1);
    
    // construct leaf to leaf mapping, based on min cost matching
    for (BpEdgeIt edge(_G); edge != lemon::INVALID; ++edge)
    {
        if (_mwm.matching(edge))
        {
            BpRedNode red = _G.redNode(edge);
            BpBlueNode blue = _G.blueNode(edge);
            
            if (_dummy[red] || _dummy[blue])
                continue;
            
            int i1 = _redNodeToIdx[red];
            int i2 = _blueNodeToIdx[blue];
            
            int mapped_i1 = i1 + _k1 - 1;
            int mapped_i2 = i2 + _k2 - 1;
            
            _leaf1ToLeaf2[mapped_i1] = mapped_i2;
            _leaf2ToLeaf1[mapped_i2] = mapped_i1;
            
            _mappedLeaves1.insert(mapped_i1);
            _mappedLeaves2.insert(mapped_i2);
        }
    }
    
    // check #segments per chromosome
    const IntArray& n1 = _T1.getTree().n();
    const IntArray& n2 = _T2.getTree().n();
    for (int chr = 0; chr < _T1.getTree().numChr(); ++chr)
    {
        if (n1[chr] != n2[chr])
        {
            return false;
        }
    }
    
    return true;
}

void Comparison::printBpGraph(std::ostream &out) const
{
    out << "graph G {" << std::endl;
    
    out << "\t" << "subgraph red {" << std::endl;
    for (BpRedNodeIt v(_G); v != lemon::INVALID; ++v)
    {
        int idx = _redNodeToIdx[v];
        const CopyNumberTree::ProfileVector& profile = _T1.getTree().profile(idx + _k1 - 1);
        
        out << "\t\t" << _G.id((BpNode)v) << " [label=\"";
        
        if (_dummy[v])
        {
            out << "dummy\",style=dashed";
        }
        else
        {
            bool first = true;
            for (int chr = 0; chr < _T1.getTree().numChr(); ++chr)
            {
                if (first)
                    first = false;
                else
                    out << "|";
                
                bool first2 = true;
                for (int s = 0; s < _T1.getTree().n()[chr]; ++s)
                {
                    if (first2)
                        first2 = false;
                    else
                        out << " ";
                    
                    out << profile[chr][s];
                }
            }
            out << "\"";
        }
        out << "]" << std::endl;
    }
    out << "\t}" << std::endl;
    
    out << "\t" << "subgraph blue {" << std::endl;
    for (BpBlueNodeIt v(_G); v != lemon::INVALID; ++v)
    {
        int idx = _blueNodeToIdx[v];
        const CopyNumberTree::ProfileVector& profile = _T2.getTree().profile(idx + _k2 - 1);
        
        out << "\t\t" << _G.id((BpNode)v) << " [label=\"";
        if (_dummy[v])
        {
            out << "dummy\",style=dashed";
        }
        else
        {
            bool first = true;
            for (int chr = 0; chr < _T2.getTree().numChr(); ++chr)
            {
                if (first)
                    first = false;
                else
                    out << "|";
                
                bool first2 = true;
                for (int s = 0; s < _T2.getTree().n()[chr]; ++s)
                {
                    if (first2)
                        first2 = false;
                    else
                        out << " ";
                    
                    out << profile[chr][s];
                }
            }
            out << "\"";
        }
        out << "]" << std::endl;

    }
    out << "\t}" << std::endl;
    
    for (BpEdgeIt e(_G); e != lemon::INVALID; ++e)
    {
        out << "\t" << _G.id(_G.u(e)) << " -- " << _G.id(_G.v(e))
            << " [label=\"" << _weight[e] << "\"";
        
        if (_mwm.matching(e) && !_dummy[_G.u(e)] && !_dummy[_G.v(e)])
        {
            out << ",color=red";
        }
        out << "]" << std::endl;
    }
    
    out << "}" << std::endl;
}

double Comparison::robinsonFoulds() const
{
    std::set<IntSetPair> T1_splits;
    for (int i1 = 0; i1 < _k1 - 1; ++i1)
    {
        // Map leaves
        IntSetPair tmp = _T1.getTree().splits(i1);
        IntSetPair splitPair1;
        
        std::set_intersection(tmp.first.begin(), tmp.first.end(),
                              _mappedLeaves1.begin(), _mappedLeaves1.end(),
                              std::inserter(splitPair1.first, splitPair1.first.begin()));
        
        std::set_intersection(tmp.second.begin(), tmp.second.end(),
                              _mappedLeaves1.begin(), _mappedLeaves1.end(),
                              std::inserter(splitPair1.second, splitPair1.second.begin()));
        
        T1_splits.insert(splitPair1);
    }

    std::set<IntSetPair> T2_splits;
    for (int i2 = 0; i2 < _k2 - 1; ++i2)
    {
        IntSetPair tmp = _T2.getTree().splits(i2);
        IntSetPair splitPair2;
        
        std::set_intersection(tmp.first.begin(), tmp.first.end(),
                              _mappedLeaves2.begin(), _mappedLeaves2.end(),
                              std::inserter(splitPair2.first, splitPair2.first.begin()));
        
        std::set_intersection(tmp.second.begin(), tmp.second.end(),
                              _mappedLeaves2.begin(), _mappedLeaves2.end(),
                              std::inserter(splitPair2.second, splitPair2.second.begin()));
        
        splitPair2.first = mapToT1(splitPair2.first);
        splitPair2.second = mapToT1(splitPair2.second);
        
        T2_splits.insert(splitPair2);
    }
    
    std::set<IntSetPair> result;
    std::set_symmetric_difference(T1_splits.begin(), T1_splits.end(),
                                  T2_splits.begin(), T2_splits.end(),
                                  std::inserter(result, result.begin()));
    
    return result.size() / (double)(_k1 - 1 + _k2 - 1);
}

IntSet Comparison::mapToT1(const IntSet& leafSet) const
{
    IntSet result;
    
    for (int i : leafSet)
    {
        if (_leaf2ToLeaf1[i] != -1)
        {
            result.insert(_leaf2ToLeaf1[i]);
        }
    }
    
    return result;
}

double Comparison::deltaM() const
{
    assert(_mappedLeaves1.size() == _mappedLeaves2.size());
    const int nrLeaves = _mappedLeaves1.size();
    
    const DoubleMatrix& M1 = _T1.getM();
    const DoubleMatrix& M2 = _T2.getM();
    
    assert(M1.size() == M2.size());
    const int m = M1.size();
    
    double delta = 0;
    for (int i1 = _k1 - 1; i1 < 2*_k1 - 1; ++i1)
    {
        int mapped_i1 = i1 - (_k1 - 1);
        if (_leaf1ToLeaf2[i1] != -1)
        {
            int i2 = _leaf1ToLeaf2[i1];
            int mapped_i2 = i2 - (_k2 - 1);
            
            for (int p = 0; p < m; ++p)
            {
                delta += fabs(M1[p][mapped_i1] - M2[p][mapped_i2]);
            }
        }
    }
    
    return delta / (m*nrLeaves);
}
