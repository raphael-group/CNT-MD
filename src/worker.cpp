#include "worker.h"

Worker::Worker(const InputInstance &inputInstance,
               const int k,
               const IntMatrix& e,
               const int Z,
               const bool forceDiploid,
               const bool rootNotFixed,
               const unsigned int iterConvergence,
               const unsigned int maxIter,
               const int timeLimit,
               const int memoryLimit,
               const int nrThreads,
               const DoubleMatrix& M0,
               const int seedIndex,
               const HotStart& inputCompleteHotStart)
    : _inputInstance(inputInstance)
    , _k(k)
    , _e(e)
    , _Z(Z)
    , _forceDiploid(forceDiploid)
    , _rootNotFixed(rootNotFixed)
    , _iterConvergence(iterConvergence)
    , _maxIter(maxIter)
    , _timeLimit(timeLimit)
    , _memoryLimit(memoryLimit)
    , _nrThreads(nrThreads)
    , _M0(M0)
    , _allC()
    , _allM()
    , _allObjC()
    , _allObjM()
    , _allTrees()
    , _seedIndex(seedIndex)
    , _inputCompleteHotStart(inputCompleteHotStart)
    , _firstCompleteHotStart()
    , _lastCompleteHotStart()
{
}


Worker::Worker(const InputInstance &inputInstance,
               const int k,
               const IntMatrix& e,
               const int Z,
               const bool forceDiploid,
               const bool rootNotFixed,
               const unsigned int iterConvergence,
               const unsigned int maxIter,
               const int timeLimit,
               const int memoryLimit,
               const int nrThreads,
               const DoubleMatrix& M0,
               const int seedIndex)
    : Worker(inputInstance, k, e, Z, forceDiploid, rootNotFixed, iterConvergence,
             maxIter, timeLimit, memoryLimit, nrThreads, M0, seedIndex, CArchitect::firstCompleteHotStart(inputInstance, e, k))
{
}


double Worker::solve()
{
    unsigned int iter_convergence = 0;
    unsigned int iter = 0;
    bool first = true;
    
    HotStart completeHotStart = _inputCompleteHotStart;
    
    while((iter_convergence < _iterConvergence) && (iter < _maxIter))
    {
        g_mutex.lock();
        CArchitect carch(_inputInstance,
                         _allM.empty() ? _M0 : _allM.back(),
                         _e, _Z, _k, _rootNotFixed, _forceDiploid);
        g_mutex.unlock();
        
        try {
            carch.init();
            carch.addCompleteHotStart(completeHotStart);
            bool status = carch.solve(_timeLimit, _memoryLimit, _nrThreads);
            assert(status);
        } catch (IloException &e) {
            std::cerr << "ILOG exception: "<< e.getMessage() << std::endl;
            e.end();
            abort();
        }
        assert(first || (g_tol.less(carch.getObjValue(), _allObjM.back()) | !g_tol.different(carch.getObjValue(), _allObjM.back())));
        assert(first || (g_tol.less(carch.getObjValue(), _allObjC.back()) | !g_tol.different(carch.getObjValue(), _allObjC.back())));

        _allObjC.push_back(carch.getObjValue());
        _allC.push_back(carch.getC());
        _allTrees.push_back(carch.getTree());
        completeHotStart = carch.getCompleteHotStart();

        if(first)
        {
            _firstCompleteHotStart = completeHotStart;
            first = false;
        }
        assert(!first);

        if (g_verbosity >= VerbosityLevel::VERBOSE_NON_ESSENTIAL)
        {
            g_output_mutex.lock();
            std::cout << _k << "\t" << _Z << "\t" << _seedIndex << "\t"
                      << iter << "\t" << "C" << "\t" << carch.getLB()
                      << "\t" << carch.getUB() << "\t" << carch.getTime()
                      << "\t" << carch.getDelta() << "\t" << carch.getObjValue() << std::endl;
            g_output_mutex.unlock();
        }
        
        g_mutex.lock();
        MArchitect march(_inputInstance, _allC.back(), _k, _e);
        g_mutex.unlock();

        try {
            march.init();
            bool status = march.solve(_timeLimit, _memoryLimit);
            assert(status);
        } catch (IloException &e) {
            std::cerr << "ILOG exception: "<< e.getMessage() << std::endl;
            e.end();
            abort();
        }
        assert(g_tol.less(march.getObjValue(), _allObjC.back()) | !g_tol.different(carch.getObjValue(), _allObjC.back()));
        
        _allM.push_back(march.getM());
        _allObjM.push_back(march.getObjValue());
        
        if(g_tol.different(_allObjM.back(), _allObjC.back())) {
            iter_convergence = 0;
        } else {
            ++iter_convergence;
        }
        
        if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
        {
            g_output_mutex.lock();
            std::cout << _k << "\t" << _Z << "\t" << _seedIndex << "\t"
                      << iter << "\t" << "M" << "\t" << march.getLB()
                      << "\t" << march.getUB() << "\t" << march.getTime()
                      << "\t" << carch.getDelta() << "\t" << march.getObjValue() << std::endl;
            g_output_mutex.unlock();
        }
        
        ++iter;
    }
    assert(iter_convergence >= _iterConvergence | iter == _maxIter);
    
    _lastCompleteHotStart = completeHotStart;
    return _allObjM.back();
}




