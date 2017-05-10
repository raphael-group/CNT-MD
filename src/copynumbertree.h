#ifndef COPYNUMBERTREE_H
#define COPYNUMBERTREE_H

#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <vector>
#include <set>
#include <cassert>
#include <ostream>
#include "basic_types.h"

class CopyNumberTree
{
public:
    typedef lemon::ListDigraph Digraph;
    DIGRAPH_TYPEDEFS(Digraph);

    typedef std::vector<Node> NodeVector;
    typedef std::vector<int> Profile;
    typedef std::vector<Profile> ProfileVector;
    typedef Digraph::NodeMap<ProfileVector> ProfileNodeMap;
    typedef lemon::DynArcLookUp<Digraph> ArcLookUp;
    
    CopyNumberTree();

    CopyNumberTree(const CopyNumberTree& other);

    CopyNumberTree& operator=(const CopyNumberTree& other);
    
    CopyNumberTree(int k, int num_chr, const IntArray& n);
    
    struct Event
    {
        Event()
            : _chr(-1)
            , _s(-1)
            , _t(-1)
            , _b(-1)
        {
        }
        
        Event(int chr, int s, int t, int b)
            : _chr(chr)
            , _s(s)
            , _t(t)
            , _b(b)
        {
        }
        
        int _chr;
        int _s;
        int _t;
        int _b;
    };
    
    typedef std::vector<Event> EventVector;
    typedef Digraph::ArcMap<EventVector> EventVectorArcMap;
    
    bool isLeaf(int i) const;
    bool isArc(int i, int j) const;
    int parent(int i) const;
    int leftChild(int i) const;
    int rightChild(int i) const;
    
    void addArc(int i, int j);
    void setProfile(int i, const ProfileVector& y_i);
    void addEvent(int chr, int i, int j, const Event& event);
    
    int index(Node v_i) const
    {
        assert(v_i != lemon::INVALID);
        return _node2idx[v_i];
    }
    
    const ProfileVector& profile(int i) const
    {
        assert(0 <= i && i < _numVertices);
        Node v_i = _idx2node[i];
        assert(v_i != lemon::INVALID);
        return _profile[v_i];
    }

    const Int3Array getLeafProfiles() const;
    
    void init();
    
    void writeDOT(std::ostream& out) const;
    
    const IntArray& n() const
    {
        return _n;
    }
    
    const Digraph& T() const
    {
        return _T;
    }
    
    int k() const
    {
        return _k;
    }
    
    int numChr() const
    {
        return _num_chr;
    }
    
    int cost() const
    {
        int res = 0;
        for (ArcIt a_ij(_T); a_ij != lemon::INVALID; ++a_ij)
        {
            for (const Event& event : _events[a_ij])
            {
                res += abs(event._b);
            }
        }
        assert(res == _numEvents);
        return res;
    }
    
    IntSetPair splits(int i) const;
    
    static std::string profileToString(const ProfileVector& y);
    
    friend std::ostream& operator<<(std::ostream& out, const CopyNumberTree& T);
    friend std::istream& operator>>(std::istream& in, CopyNumberTree& T);
    
private:
    int _k;
    int _num_chr;
    IntArray _n;
    int _numVertices;
    Digraph _T;
    ArcLookUp _arcLookUp;
    IntNodeMap _node2idx;
    NodeVector _idx2node;
    ProfileNodeMap _profile;
    EventVectorArcMap _events;
    int _numEvents;
    
    void leafSet(const Node v_i, IntSet& leafSet) const;
};

std::ostream& operator<<(std::ostream& out, const CopyNumberTree& T);
std::istream& operator>>(std::istream& in, CopyNumberTree& T);

#endif //COPYNUMBERTREE_H
