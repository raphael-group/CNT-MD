#include "manager.h"

Manager::Manager(const InputInstance &inputInstance,
                 const unsigned int k,
                 const IntMatrix &e,
                 const unsigned int LZ,
                 const unsigned int UZ,
                 const bool forceDiploid,
                 const bool rootNotFixed,
                 const bool deactiveRefinement,
                 const int size_bubbles,
                 const unsigned int iterConvergence,
                 const unsigned int maxIter,
                 const int nrSeeds,
                 const int nrWorkers,
                 const int nrILPthreads,
                 const int timeLimit,
                 const int memoryLimit,
                 const double eps)
    : _inputInstance(inputInstance)
    , _k(k)
    , _e(e)
    , _LZ(LZ)
    , _UZ(UZ)
    , _forceDiploid(forceDiploid)
    , _rootNotFixed(rootNotFixed)
    , _deactiveRefinement(deactiveRefinement)
    , _size_bubbles(size_bubbles)
    , _iterConvergence(iterConvergence)
    , _maxIter(maxIter)
    , _nrSeeds(nrSeeds)
    , _nrWorkers(nrWorkers)
    , _nrILPthreads(nrILPthreads)
    , _allM0()
    , _timeLimit(timeLimit)
    , _memoryLimit(memoryLimit)
    , _eps(eps)
    , _sem(nrWorkers)
    , _threadGroup()
    , _mutex()
    , _bestObjValue()
    , _bestT()
    , _bestC()
    , _bestM()
    , _bestZ(-1)
    , _firstCompleteHotStart()
    , _lastCompleteHotStart()
    , _isComputed()
    , _norm(0)
    , _diploidCompleteHotStart()
    , _refinedTree()
    , _refinedObjValue()
{
    for (int i = 0; i < _nrSeeds; ++i)
    {
        _allM0.push_back(build_random_M(_inputInstance.m(), _k, _size_bubbles));
    }

    //_LZ = std::min(_LZ, _UZ);

    _norm = countElements(_inputInstance.F());
    assert(_norm == (sum_of_elements(_inputInstance.n()) * inputInstance.m()));

    _diploidCompleteHotStart = CArchitect::firstCompleteHotStart(_inputInstance, _e, _k);
}


void Manager::runBinarySearch()
{
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() << "Run in Binary Search mode" << std::endl;

    initialize();

    int lb = _LZ;
    int ub = _UZ;

    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() << "Computing the distance for the right-most value UZ=" << _UZ << " in the starting interval [LZ=" << lb << ", UZ=" << ub << "]" << std::endl;

    computeDistance(_UZ);

    while((ub - lb) > 1)
    {        
        int mid = std::ceil(((double)ub + (double)lb)/2.0);
        assert(mid != lb);
        assert(mid != ub);

        if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
            std::cerr << timestamp() << "Considering interval [LZ=" << lb << ", UZ=" << ub << "] with MZ=" << mid << std::endl;

        if(!_isComputed[mid])
            computeDistance(mid);

        if(isImproving(_bestObjValue[mid], _bestObjValue[_UZ]))
        {
            assert(g_tol.less(_bestObjValue[_UZ], _bestObjValue[mid]));
            if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
                std::cerr << timestamp() << "The distance d(" << mid << ")=" << _bestObjValue[mid] << " is improved by " << "d(" << _UZ << ")="
                          << _bestObjValue[_UZ] << " with a total tolerance of " << (_eps * _norm) << std::endl;

            lb = mid;
            ub = ub;
        } else {
            if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
                std::cerr << timestamp() << "The distance d(" << mid << ")=" << _bestObjValue[mid] << " is not improved by " << "d(" << _UZ << ")="
                          << _bestObjValue[_UZ] << " with a total tolerance of " << (_eps * _norm) << std::endl;

            lb = lb;
            ub = mid;
        }

        assert((lb < mid & ub == mid) | (ub > mid & lb == mid));
    }

    if((ub - lb) == 1)
    {
        if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
            std::cerr << timestamp() << "Considering interval [LZ=" << lb << ", UZ=" << ub << "]" << std::endl;

        if(!_isComputed[lb])
            computeDistance(lb);

        if(isImproving(_bestObjValue[lb], _bestObjValue[_UZ]))
        {
            if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
                std::cerr << timestamp() << "The distance d(" << lb << ")=" << _bestObjValue[lb] << " is improved by " << "d(" << _UZ << ")="
                          << _bestObjValue[_UZ] << " with a total tolerance of " << (_eps * _norm) << std::endl;

            if(!_isComputed[ub])
                computeDistance(ub);

            _bestZ = ub;

        } else {
            if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
                std::cerr << timestamp() << "The distance d(" << lb << ")=" << _bestObjValue[lb] << " is not improved by " << "d(" << _UZ << ")="
                          << _bestObjValue[_UZ] << " with a total tolerance of " << (_eps * _norm) << std::endl;

            _bestZ = lb;
        }
    } else {
        if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
            std::cerr << timestamp() << "Considering interval [LZ=" << lb << ", UZ=" << ub << "]" << std::endl;

        assert(ub == lb);
        if(!_isComputed[lb])
            computeDistance(lb);

        _bestZ = lb;
    }

    assert(_isComputed[_bestZ]);
    assert(_bestZ >= _LZ & _bestZ <= _UZ);
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() << "The selected slope is " << _bestZ << std::endl;

    refinement();
}


void Manager::runReverse()
{
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() << "Run in Reverse mode" << std::endl;

    initialize();

    int Z = _UZ;
    computeDistance(Z);

    if(Z > _LZ)
    {
        do {
            std::cerr << timestamp() << "The distance d(" << Z << ")=" << _bestObjValue[Z] << " is not improved by " << "d(" << _UZ << ")="
                      << _bestObjValue[_UZ] << " with a total tolerance of " << (_eps * _norm) << std::endl;

            --Z;
            computeDistance(Z, _diploidCompleteHotStart);

        } while(!isImproving(_bestObjValue[Z], _bestObjValue[_UZ]) && Z > _LZ);
        _bestZ = Z + 1;
    } else {
        _bestZ = Z;
    }

    assert(_bestZ >= 0 & _bestZ <= _UZ);
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() << "The selected slope is " << _bestZ << std::endl;

    refinement();
}


void Manager::runIterative()
{
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() << "Run in Iterative mode" << std::endl;

    initialize();

    for (int Z = _LZ; Z <= _UZ; ++Z)
    {
        computeDistance(Z);
    }

    int Z = _UZ;
    while(!isImproving(_bestObjValue[Z-1], _bestObjValue[_UZ]) && Z > _LZ)
    {
        std::cerr << timestamp() << "The distance d(" << Z-1 << ")=" << _bestObjValue[Z-1] << " is not improved by " << "d(" << _UZ << ")="
                  << _bestObjValue[_UZ] << " with a total tolerance of " << (_eps * _norm) << std::endl;

        --Z;
    }
    _bestZ = Z;

    assert(_bestZ >= 0 & _bestZ <= _UZ);
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() << "The selected slope is " << _bestZ << std::endl;

    refinement();
}


DoubleMatrix Manager::build_random_M(const int num_samples, const int num_leaves, const int size_bubbles)
{
    DoubleMatrix result;

    for(unsigned int p = 0; p < num_samples; ++p)
    {
        std::uniform_int_distribution<int> uni(1, num_leaves);

        //const int num_parts = generator();
        const int num_parts = std::max(uni(g_rng), uni(g_rng));
        result.push_back(build_partition_vector(num_leaves, num_parts, size_bubbles));
    }

    return result;
}


DoubleArray Manager::build_partition_vector(const int num_leaves, const int num_parts, const int size_bubbles)
{
    IntArray positions(num_leaves);
    std::iota(positions.begin(), positions.end(), 0);
    std::shuffle(positions.begin(), positions.end(), g_rng);

    //std::uniform_int_distribution<int> uni(1, size_bubbles);
    std::uniform_int_distribution<> dist(0, size_bubbles);

    DoubleArray bubbles;
    for(unsigned int i = 0; i < num_parts-1; ++i)
    {
        bubbles.push_back((double)dist(g_rng)/(double)size_bubbles);
    }
    std::sort(bubbles.begin(), bubbles.end());

    DoubleArray result(num_leaves, 0.0);
    if(num_parts > 1)
    {
        for(unsigned int i = 0; i < num_parts; ++i)
        {
            if(i == 0)
            {
                result[positions[i]] = bubbles[i];
            }
            else if(i == (num_parts-1))
            {
                result[positions[i]] = 1.0 - bubbles[i-1];
            }
            else
            {
                result[positions[i]] = bubbles[i] - bubbles[i-1];
            }
        }
    }
    else {
        result[positions[0]] = 1.0;
    }

    return result;
}


/**
DoubleMatrix Manager::build_random_M(int m, int k, double alpha)
{
    DoubleMatrix M(m, DoubleArray(k));
    
    boost::gamma_distribution<> dist(alpha);
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > generator(g_rng, dist);
    
    for(unsigned int p = 0; p < m; ++p)
    {
        double row_sum = 0.0;
        for(unsigned int i = 0; i < k; ++i)
        {
            double value = generator();
            M[p][i] = value;
            row_sum += value;
        }
        
        for(unsigned int i = 0; i < k; ++i)
        {
            M[p][i] = M[p][i] / row_sum;
        }
    }

    return M;
}
**/


void Manager::runInstance(const int Z, const int seedIdx, const HotStart &inputCompleteHotStart)
{
    Worker worker(_inputInstance, _k, _e, Z,
                  _forceDiploid, _rootNotFixed,
                  _iterConvergence, _maxIter,
                  _timeLimit, _memoryLimit, _nrILPthreads,
                  _allM0[seedIdx], seedIdx,
                  inputCompleteHotStart);
    double objValue = worker.solve();

    {
        boost::interprocess::scoped_lock<boost::mutex> lock(_mutex);

        if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
            std::cerr << ":s(" << seedIdx << ")=" << objValue << ":";

        _firstCompleteHotStart[Z][seedIdx] = worker.getFirstCompleteHotStart();
        if (g_tol.less(objValue, _bestObjValue[Z]))
        {
            _bestObjValue[Z] = objValue;
            _bestT[Z] = worker.getT();
            _bestC[Z] = worker.getC();
            _bestM[Z] = worker.getM();
            _lastCompleteHotStart[Z] = worker.getLastCompleteHotStart();
            _isComputed[Z] = true;
        }
    }
    
    _sem.post();
}


void Manager::computeDistance(const int Z)
{
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() <<"Computing distance with " << Z << " maximum number of events ";
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << "{";
    assert(!_isComputed[Z]);

    initializeZwithPrevious(Z);

    for (int i = 0; i < _nrSeeds; ++i)
    {
        _sem.wait();
        _threadGroup.create_thread(boost::bind(&Manager::runInstance, this, Z, i, previousCompleteHotStart(Z, i)));
    }
    _threadGroup.join_all();

    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << "} ";
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << std::endl << timestamp() << "Distance with " << Z <<" maximum number of events is " << _bestObjValue[Z] << std::endl;
}


void Manager::computeDistance(const int Z, const HotStart &inputCompleteHotStart)
{
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << timestamp() << "Computing distance with " << Z << " maximum number of events ";
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << "{";
    assert(!_isComputed[Z]);

    initializeZwithPrevious(Z);

    for (int i = 0; i < _nrSeeds; ++i)
    {
        _sem.wait();
        _threadGroup.create_thread(boost::bind(&Manager::runInstance, this, Z, i, inputCompleteHotStart));
    }
    _threadGroup.join_all();

    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << "} ";
    if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
        std::cerr << std::endl << timestamp() << "Distance with " << Z <<" maximum number of events is " << _bestObjValue[Z] << std::endl;
}


const HotStart& Manager::previousCompleteHotStart(const int Z, const int seedIdx)
{
    int ctr = Z;

    do {
        --ctr;
    } while(ctr >= 0 && !_isComputed[ctr]);

    if(ctr >= 0)
    {
        assert(_isComputed[ctr]);
        return _firstCompleteHotStart[ctr][seedIdx];
    } else {
        assert(ctr == -1);
        return _diploidCompleteHotStart;
    }
}


inline bool Manager::isImproving(const double previous, const double successive) const
{
    return g_tol.less(_eps, (double)(previous - successive) / _norm);
}


void Manager::initialize()
{
    const unsigned int size = _UZ + 1;

    _bestObjValue = DoubleArray(size, std::numeric_limits<double>::max());
    _bestT = std::vector<CopyNumberTree>(size);
    _bestC = Int4Array(size);
    _bestM = Double3Array(size);
    _firstCompleteHotStart = std::vector<std::vector<HotStart> >(size, std::vector<HotStart> (_nrSeeds));
    _lastCompleteHotStart = std::vector<HotStart>(size);
    _isComputed = std::vector<bool>(size, false);

    if (g_verbosity >= VerbosityLevel::VERBOSE_NON_ESSENTIAL)
        std::cout << "k" << "\t" << "Z" << "\t" << "seed" << "\t"
                  << "iteration" << "\t" << "step" << "\t" << "LB"
                  << "\t" << "UB" << "\t" << "runtime"
                  << "\t" << "delta(T)" << "\t" << "ObjValue" << std::endl;
}


void Manager::initializeZwithPrevious(const int Z)
{
    int ctr = Z;

    do {
        --ctr;
    } while(ctr >= 0 && !_isComputed[ctr]);

    if(ctr >= 0)
    {
        assert(_isComputed[ctr]);

        _bestObjValue[Z] = _bestObjValue[ctr];
        _bestT[Z] = _bestT[ctr];
        _bestC[Z] = _bestC[ctr];
        _bestM[Z] = _bestM[ctr];
        _lastCompleteHotStart[Z] = _lastCompleteHotStart[ctr];
        _isComputed[Z] = true;

        if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
            std::cerr << "previousBest=" << _bestObjValue[Z] << ":";
    }
}


void Manager::refinement()
{
    if(!_deactiveRefinement)
    {
        if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
            std::cerr << timestamp() << "Starting refinement with " << _bestT[_bestZ].cost() << " events in the tree and distance of " << _bestObjValue[_bestZ] << std::endl;

        Refiner refiner(_inputInstance, _bestM[_bestZ], _e, _bestZ, _k, _rootNotFixed,
                        _forceDiploid, _bestObjValue[_bestZ] + _norm * _eps);

        bool status = false;
        try {
            refiner.init();
            refiner.addCompleteHotStart(_lastCompleteHotStart[_bestZ]);
            status = refiner.solve(_timeLimit * _maxIter, _memoryLimit * _nrWorkers, _nrILPthreads * _nrWorkers);
        } catch (IloException &e) {
            std::cerr << "ILOG exception: "<< e.getMessage() << std::endl;
            e.end();
            abort();
        }

        assert(status);
        assert(refiner.getObjValue() <= _bestZ);
        assert(g_tol.less(refiner.getDistance(), _bestObjValue[_bestZ] + _norm * _eps) ||
               !g_tol.different(refiner.getDistance(), _bestObjValue[_bestZ] + _norm * _eps));

        if(status)
        {
            _refinedObjValue = refiner.getDistance();
            _refinedTree = refiner.getTree();
        } else {
            _refinedObjValue = _bestObjValue[_bestZ];
            _refinedTree = _bestT[_bestZ];
        }

        if(g_verbosity >= VerbosityLevel::VERBOSE_ESSENTIAL)
            std::cerr << timestamp() << "The refined tree has " << _refinedTree.cost() << " events and the distance is " << _refinedObjValue << std::endl;

    } else {
        _refinedObjValue = _bestObjValue[_bestZ];
        _refinedTree = _bestT[_bestZ];
    }
}
