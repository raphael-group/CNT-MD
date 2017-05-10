#ifndef _COMPARISON_H_
#define _COMPARISON_H_

#include "fmcsolution.h"
#include "copynumbertree.h"
#include <lemon/matching.h>

class Comparison
{
public:
    Comparison(const FMCSolution& T1, const FMCSolution& T2);
    
    bool init();
    
    int leaf1ToLeaf2(int i) const
    {
        return _leaf1ToLeaf2[i];
    }
    
    int leaf2ToLeaf1(int j) const
    {
        return _leaf2ToLeaf1[j];
    }
    
    double robinsonFoulds() const;
    
    IntSet mapToT1(const IntSet& leafSet) const;
    
    void printBpGraph(std::ostream& out) const;
    
    double deltaM() const;
    
    double leafConsistency() const
    {
        assert(_mappedLeaves1.size() == _mappedLeaves2.size());
        assert(_T1.getInput().numChr() == _T2.getInput().numChr());
        return (-1.0) * _mwm.matchingWeight() / (_mappedLeaves1.size() * _T1.getInput().numChr());
    }
    
private:
    typedef lemon::ListBpGraph BpGraph;
    typedef BpGraph::RedNodeMap<int> IntRedNodeMap;
    typedef BpGraph::BlueNodeMap<int> IntBlueNodeMap;
    typedef BpGraph::EdgeMap<int> IntBpEdgeMap;
    typedef BpGraph::NodeMap<bool> BoolBpNodeMap;
    
    typedef BpGraph::RedNode BpRedNode;
    typedef BpGraph::RedNodeIt BpRedNodeIt;
    typedef BpGraph::BlueNode BpBlueNode;
    typedef BpGraph::BlueNodeIt BpBlueNodeIt;
    typedef BpGraph::Edge BpEdge;
    typedef BpGraph::EdgeIt BpEdgeIt;
    typedef BpGraph::Node BpNode;
    typedef BpGraph::NodeIt BpNodeIt;
    
    typedef std::vector<lemon::ListBpGraph::RedNode> RedNodeVector;
    typedef std::vector<lemon::ListBpGraph::BlueNode> BlueNodeVector;
    
    typedef lemon::MaxWeightedPerfectMatching<BpGraph> MWM;
    
private:
    const FMCSolution& _T1;
    const FMCSolution& _T2;
    const int _k1;
    const int _k2;
    
    BpGraph _G;
    BoolBpNodeMap _dummy;
    IntRedNodeMap _redNodeToIdx;
    IntBlueNodeMap _blueNodeToIdx;
    RedNodeVector _idxToRedNode;
    BlueNodeVector _idxToBlueNode;
    IntBpEdgeMap _weight;
    MWM _mwm;
    
    std::vector<int> _leaf1ToLeaf2;
    std::vector<int> _leaf2ToLeaf1;
    std::set<int> _mappedLeaves1;
    std::set<int> _mappedLeaves2;
};

#endif // _COMPARISON_H_
