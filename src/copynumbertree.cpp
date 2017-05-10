#include "copynumbertree.h"
#include <stdexcept>

CopyNumberTree::CopyNumberTree()
    : _k(0)
    , _num_chr(0)
    , _n(0)
    , _numVertices(0)
    , _T()
    , _arcLookUp(_T)
    , _node2idx(_T)
    , _idx2node(_numVertices, lemon::INVALID)
    , _profile(_T)
    , _events(_T)
    , _numEvents(0)
{
}

CopyNumberTree::CopyNumberTree(int k, int num_chr, const IntArray& n)
    : _k(k)
    , _num_chr(num_chr)
    , _n(n)
    , _numVertices(2*k - 1)
    , _T()
    , _arcLookUp(_T)
    , _node2idx(_T)
    , _idx2node(_numVertices, lemon::INVALID)
    , _profile(_T)
    , _events(_T)
    , _numEvents(0)
{
    init();
}

CopyNumberTree::CopyNumberTree(const CopyNumberTree &other)
    : _k(other._k)
    , _num_chr(other._num_chr)
    , _n(other._n)
    , _numVertices(other._numVertices)
    , _T()
    , _arcLookUp(_T)
    , _node2idx(_T)
    , _idx2node(other._idx2node.size(), lemon::INVALID)
    , _profile(_T)
    , _events(_T)
    , _numEvents(0)
{
    lemon::digraphCopy(other._T, _T)
            .nodeMap(other._node2idx, _node2idx)
            .nodeMap(other._profile, _profile)
            .arcMap(other._events, _events)
            .run();
    _numEvents = other._numEvents;

    for (NodeIt v(_T); v != lemon::INVALID; ++v)
    {
        int idx = _node2idx[v];
        _idx2node[idx] = v;
    }
}

CopyNumberTree& CopyNumberTree::operator=(const CopyNumberTree& other)
{
    if (&other != this)
    {
        _k = other._k;
        _num_chr = other._num_chr;
        _n = other._n;
        _numVertices = other._numVertices;
        _T.clear();
        lemon::digraphCopy(other._T, _T)
                .nodeMap(other._node2idx, _node2idx)
                .nodeMap(other._profile, _profile)
                .arcMap(other._events, _events)
                .run();
        _numEvents = other._numEvents;

        _idx2node = NodeVector(_numVertices, lemon::INVALID);
        for (NodeIt v(_T); v != lemon::INVALID; ++v)
        {
            int idx = _node2idx[v];
            _idx2node[idx] = v;
        }
    }
    return *this;
}

const Int3Array CopyNumberTree::getLeafProfiles() const
{
    Int3Array C(_num_chr);
    for(unsigned int c = 0; c < _num_chr; ++c)
    {
        for (int i = 0; i < _numVertices; ++i)
        {
            if(isLeaf(i))
            {
                const ProfileVector& y_i = profile(i);
                C[c].push_back(y_i[c]);
            }
        }
    }

    return C;
}

void CopyNumberTree::init()
{
    _T.clear();
    _numVertices = 2*_k - 1;
    _idx2node = NodeVector(_numVertices, lemon::INVALID);
    
    // 2k-1 vertices
    for (int i = 0; i < _numVertices; ++i)
    {
        Node v_i = _T.addNode();
        _idx2node[i] = v_i;
        _node2idx[v_i] = i;
        _profile[v_i] = ProfileVector(_num_chr);
        for (int chr = 0; chr < _num_chr; ++chr)
        {
            _profile[v_i][chr] = Profile(_n[chr], 0);
        }
    }
}

std::ostream& operator<<(std::ostream& out, const CopyNumberTree& T)
{
    out << "#PARAMS" << std::endl;
    out << T._num_chr << " #number of chromosomes" << std::endl;
    out << T._k << " #number of leaves" << std::endl;
    for (int n : T._n)
    {
        out << n << " ";
    }
    out << "#number of segments per chromosome" << std::endl;
    out << "#PROFILES" << std::endl;
    
    for (int i = 0; i < T._numVertices; ++i)
    {
        const CopyNumberTree::ProfileVector& y_i = T.profile(i);
        out << i << " :";
        for (int chr = 0; chr < T._num_chr; ++chr)
        {
            if (chr != 0)
            {
                out << " |";
            }
            for (int s = 0; s < T._n[chr]; ++s)
            {
                out << " " << y_i[chr][s];
            }
        }
        out << std::endl;
    }
    
    out << "#EDGES" << std::endl;
    for (CopyNumberTree::ArcIt a_ij(T._T); a_ij != lemon::INVALID; ++a_ij)
    {
        CopyNumberTree::Node v_i = T._T.source(a_ij);
        CopyNumberTree::Node v_j = T._T.target(a_ij);
        out << T._node2idx[v_i] << " -> " << T._node2idx[v_j] << std::endl;
    }
    
    out << "#EVENTS" << std::endl;
    for (CopyNumberTree::ArcIt a_ij(T._T); a_ij != lemon::INVALID; ++a_ij)
    {
        CopyNumberTree::Node v_i = T._T.source(a_ij);
        CopyNumberTree::Node v_j = T._T.target(a_ij);
        
        const CopyNumberTree::EventVector& events_ij = T._events[a_ij];
        for (const CopyNumberTree::Event& event : events_ij)
        {
            out << event._chr << " "
                << T._node2idx[v_i] << " " << T._node2idx[v_j] << " "
                << event._s << " "
                << event._t << " " << event._b << std::endl;
        }
    }
    
    return out;
}

std::istream& operator>>(std::istream& in, CopyNumberTree& T)
{
    T._n.clear();
    
    std::string line;
    std::string value;
    
    std::getline(in, line, '\n'); //Skip the first line "#PARAMS"
    
    /** Read the number of chromosomes **/
    std::getline(in, line, '\n');
    std::stringstream chromosomeline(line);
    chromosomeline >> T._num_chr;
    
    /** Read the number of vertices **/
    std::getline(in, line, '\n');
    std::stringstream leafline(line);
    leafline >> T._k;
    
    /** Read the number of segments for each chromosome  **/
    std::getline(in, line, '\n');
    std::stringstream segline(line);
    std::getline(segline, line, '#');
    std::stringstream split_seg(line);
    while (!split_seg.eof())
    {
        std::getline(split_seg, value, ' ');
        if(!value.empty())
        {
            T._n.push_back((unsigned int)atoi(value.c_str()));
        }
    }
    
    if(T._n.size() != T._num_chr)
    {
        throw std::runtime_error(std::string("ERROR: inconsistent numbers of segments in #PARAMS"));
    }
    
    T.init();
    
    getline(in, line, '\n'); //Skip the line "#PROFILES"
    
    /** Begin the reading of the profiles **/
    unsigned int counter_chromosomes = 0;
    unsigned int counter_vertices = 0;
    unsigned int counter_seg = 0;
    
    while (counter_vertices < T._numVertices && in.good())
    {
        getline(in, line, '\n');
        if(!line.empty() && !(line == "#EDGES"))
        {
            std::stringstream csline(line);
            getline(csline, line, ':'); //Skip the name Li of the leaf
            getline(csline, line, '\n');
            std::stringstream sline(line);
            
            counter_chromosomes = 0;
            while(!sline.eof())
            {
                std::string leaf;
                std::getline(sline, leaf, '|');
                std::stringstream sleaf(leaf);
                
                counter_seg = 0;
                while(!sleaf.eof())
                {
                    getline(sleaf, value, ' ');
                    if(!value.empty())
                    {
                        if(counter_vertices >= T._numVertices
                           || counter_chromosomes >= T._num_chr
                           || counter_seg >= T._n[counter_chromosomes])
                        {
                            throw std::runtime_error("ERROR: the input format is wrong or the number of chromosomes/leaves/corresponding segments is inconsistent with the specified PARAMS");
                        }
                        else
                        {
                            CopyNumberTree::Node v_i = T._idx2node[counter_vertices];
                            T._profile[v_i][counter_chromosomes][counter_seg] = atoi(value.c_str());
                            
                        }
                        ++counter_seg;
                    }
                }
                
                if(counter_seg != T._n[counter_chromosomes])
                {
                    std::stringstream tmp;
                    tmp << "ERROR: inconsistent number of segments for chromosome " << counter_chromosomes;
                    
                    throw std::runtime_error(tmp.str());
                }
                
                ++counter_chromosomes;
            }
            
            if(counter_chromosomes != T._num_chr)
            {
                throw std::runtime_error("ERROR: inconsistent number of chromosomes");
            }
            
            ++counter_vertices;
        }
    }
    
    if(counter_vertices != T._numVertices)
    {
        throw std::runtime_error("ERROR: inconsistent number of leaves");
    }
    
    getline(in, line, '\n'); //Skip the line "#EDGES"
    
    unsigned int counter_edges = 0;
    while (counter_edges < T._numVertices - 1 && in.good())
    {
        getline(in, line, '\n');
        if(!line.empty() && !(line == "#EVENTS"))
        {
            std::stringstream csline(line);
            std::string tmp;
            int i = -1, j = -1;
            csline >> i >> tmp >> j;
            
            if (!(0 <= i && i < T._numVertices))
            {
                throw std::runtime_error("ERROR: incorrect arc (i)");
            }
            if (!(0 <= j && j < T._numVertices))
            {
                throw std::runtime_error("ERROR: incorrect arc (j)");
            }
            
            T.addArc(i, j);
            ++counter_edges;
        }
    }
    
    if (counter_edges != T._numVertices - 1)
    {
        throw std::runtime_error("ERROR: inconsistent number of edges");
    }
    
    getline(in, line, '\n'); //Skip the line "#EVENTS"
    while (in.good())
    {
        getline(in, line, '\n');
        if (line.empty())
            continue;
        if (line == "#MIX-PARAMS")
            break;
        
        std::stringstream csline(line);
        
        int i = -1, j = -1, chr = -1, s = -1, t = -1, b = 0;
        csline >> chr >> i >> j >> s >> t >> b;
        
        if (!(0 <= i && i < T._numVertices))
        {
            throw std::runtime_error("ERROR: incorrect arc (i) for event");
        }
        if (!(0 <= j && j < T._numVertices))
        {
            throw std::runtime_error("ERROR: incorrect arc (j) for event");
        }
        
        if (!(0 <= chr && chr < T._num_chr))
        {
            throw std::runtime_error("ERROR: incorrect chromosome for event");
        }
        
        if (!(0 <= s && s <= t && t < T._n[chr]))
        {
            throw std::runtime_error("ERROR: incorrect positions for event");
        }
        
        T.addEvent(chr, i, j, CopyNumberTree::Event(chr, s, t, b));
    }
    
    return in;
}

bool CopyNumberTree::isLeaf(int i) const
{
    assert(0 <= i && i < _numVertices);
    return i >= _k - 1;
}

bool CopyNumberTree::isArc(int i, int j) const
{
    assert(0 <= i && i < _numVertices);
    assert(0 <= j && j < _numVertices);
    Node v_i = _idx2node[i];
    Node v_j = _idx2node[j];
    
    return _arcLookUp(v_i, v_j) != lemon::INVALID;
}

int CopyNumberTree::parent(int i) const
{
    assert(0 <= i && i < _numVertices);

    Node v_i = _idx2node[i];
    assert(v_i != lemon::INVALID);
    
    InArcIt a_ij(_T, v_i);
    if (a_ij == lemon::INVALID)
    {
        // this can only be the case for root
        assert(i == 0);
        return -1;
    }
    else
    {
        Node v_pi_i = _T.source(a_ij);
        return _node2idx[v_pi_i];
    }
}

int CopyNumberTree::leftChild(int i) const
{
    assert(0 <= i && i < _numVertices);
    if (isLeaf(i))
    {
        return -1;
    }
    else
    {
        Node v_i = _idx2node[i];
        assert(v_i != lemon::INVALID);
        OutArcIt a_ij(_T, v_i);
        assert(a_ij != lemon::INVALID);
        return _node2idx[_T.target(a_ij)];
    }
}

int CopyNumberTree::rightChild(int i) const
{
    assert(0 <= i && i < _numVertices);
    if (isLeaf(i))
    {
        return -1;
    }
    else
    {
        Node v_i = _idx2node[i];
        assert(v_i != lemon::INVALID);
        OutArcIt a_ij(_T, v_i);
        assert(a_ij != lemon::INVALID);
        ++a_ij;
        assert(a_ij != lemon::INVALID);
        return _node2idx[_T.target(a_ij)];
    }
}

void CopyNumberTree::addArc(int i, int j)
{
    assert(0 <= i && i < _numVertices);
    assert(0 <= j && j < _numVertices);
    Node v_i = _idx2node[i];
    Node v_j = _idx2node[j];
    
    assert(v_i != lemon::INVALID);
    assert(v_j != lemon::INVALID);
    
    _T.addArc(v_i, v_j);
}

void CopyNumberTree::setProfile(int i, const ProfileVector& y_i)
{
    assert(0 <= i && i < _numVertices);
    
    Node v_i = _idx2node[i];
    assert(v_i != lemon::INVALID);
    
    _profile[v_i] = y_i;
}

void CopyNumberTree::addEvent(int chr, int i, int j, const Event& event)
{
    assert(0 <= i && i < _numVertices);
    assert(0 <= j && j < _numVertices);
    assert(0 <= chr && chr < _num_chr);
    
    assert(0 <= event._s && event._s <= event._t && event._t < _n[chr]);
    assert(!isLeaf(i));
    
    Node v_i = _idx2node[i];
    Node v_j = _idx2node[j];
    
    assert(v_i != lemon::INVALID);
    assert(v_j != lemon::INVALID);
    
    Arc a_ij = _arcLookUp(v_i, v_j);
    assert(a_ij != lemon::INVALID);
    
    _events[a_ij].push_back(event);

    _numEvents += abs(event._b);
}

void CopyNumberTree::writeDOT(std::ostream &out) const
{
    out << "digraph T {" << std::endl;
    
    for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
    {
        const int i = _node2idx[v_i];
        const ProfileVector& y_i = _profile[v_i];
        
        out << "\t" << i << " [label=\"";
        
        bool first = true;
        for (int chr = 0; chr < _num_chr; ++chr)
        {
            if (first)
                first = false;
            else
                out << "|";
            
            bool first2 = true;
            for (int s = 0; s < _n[chr]; ++s)
            {
                if (first2)
                    first2 = false;
                else
                    out << " ";
                
                out << y_i[chr][s];
            }
        }
        out << "\"]" << std::endl;
    }
    
    for (ArcIt a_ij(_T); a_ij != lemon::INVALID; ++a_ij)
    {
        Node v_i = _T.source(a_ij);
        Node v_j = _T.target(a_ij);
        
        const int i = _node2idx[v_i];
        const int j = _node2idx[v_j];
        
        out << "\t" << i << " -> " << j << " [label=\"";
        
        bool first = true;
        for (const auto& event : _events[a_ij])
        {
            if (first)
                first = false;
            else
                out << "\n";
            
            out << "(" << event._chr << "," << event._s << "," << event._t << "," << event._b << ")";
        }
        out << "\"]" << std::endl;
    }
    out << "}" << std::endl;
}

IntSetPair CopyNumberTree::splits(int i) const
{
    IntSetPair splitSetPair;
    
    leafSet(_idx2node[i], splitSetPair.first);
    for (int i = _k - 1; i < 2 * _k - 1; ++i)
    {
        if (splitSetPair.first.count(i) == 0)
        {
            splitSetPair.second.insert(i);
        }
    }
    
    return splitSetPair;
}

void CopyNumberTree::leafSet(const Node v_i, IntSet& set) const
{
    const int i = _node2idx[v_i];
    if (isLeaf(i))
    {
        set.insert(i);
    }
    else
    {
        for (OutArcIt a_ij(_T, v_i); a_ij != lemon::INVALID; ++a_ij)
        {
            Node v_j = _T.target(a_ij);
            leafSet(v_j, set);
        }
    }
}

std::string CopyNumberTree::profileToString(const ProfileVector &y)
{
    std::stringstream ss;
 
    bool first = true;
    for (int chr = 0; chr < y.size(); ++chr)
    {
        if (first)
            first = false;
        else
            ss << "|";
        
        bool first2 = true;
        for (int s = 0; s < y[chr].size(); ++s)
        {
            if (first2)
                first2 = false;
            else
                ss << " ";
            
            ss << y[chr][s];
        }
    }

    return ss.str();
}
