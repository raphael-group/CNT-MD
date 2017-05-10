#ifndef _MANAGER_H_
#define _MANAGER_H_

#include "basic_types.h"
#include "fmcsolution.h"
#include "copynumbertree.h"
#include "worker.h"
#include "refiner.h"

#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/thread.hpp>


class Manager
{
public:
    Manager(const InputInstance &inputInstance,
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
            const double eps);
    
public:
    /// Run with the binary search method
    void runBinarySearch();
    /// Run with the reverse iterative mode from UB to LB until best is found
    void runReverse();
    /// Run with the iterative mode from LB to UB
    void runIterative();
    /// Construct random M
    static DoubleMatrix build_random_M(const int num_samples, const int num_leaves, const int size_bubbles);
    /// Construct random vector summing up to 1
    static DoubleArray build_partition_vector(const int num_leaves, const int num_parts, const int size_bubbles);
    /// Get best solution
    FMCSolution getSolution() const
    {
        return FMCSolution(_refinedTree, _bestM[_bestZ], _inputInstance);
    }
    /// Get best objective value
    double getObjValue() const
    {
        return _refinedObjValue;
    }
    /// Get the number of events selected as the best slope point by the running procedure
    int getSlopePoint() const
    {
        return _bestZ;
    }
    
private:
    /// Input instance
    const InputInstance& _inputInstance;
    /// Number of leaves
    const unsigned int _k;
    /// Maximum copy number per chromosome, per position
    const IntMatrix& _e;
    /// Lower bound on Z
    const int _LZ;
    /// Upper bound on Z
    const int _UZ;
    /// Force the presence of the normal diploid clone
    const bool _forceDiploid;
    /// Do not fix the root to the normal diploid
    const bool _rootNotFixed;
    /// Deactive refinement
    const bool _deactiveRefinement;
    /// Dirichlet parameter for generating initial M
    const int _size_bubbles;
    /// Number of iterations for checking convergence
    const unsigned int _iterConvergence;
    /// Maximum number of iterations for each seed
    const unsigned int _maxIter;
    /// Number of seeds
    const int _nrSeeds;
    /// Number of workers
    const int _nrWorkers;
    /// Number of LPthreads
    const int _nrILPthreads;
    /// All initial seeds
    std::vector<DoubleMatrix> _allM0;
    /// Time limit (seconds)
    const int _timeLimit;
    /// Memory limit (MB)
    const int _memoryLimit;
    /// Epsilon, the threshold estabilishing the tolerance on the objective value
    const double _eps;
    /// Semaphore
    boost::interprocess::interprocess_semaphore _sem;
    /// Thread group
    boost::thread_group _threadGroup;
    /// Mutex
    boost::mutex _mutex;
    /// Best objective value
    DoubleArray _bestObjValue;
    /// Best copy-number tree
    std::vector<CopyNumberTree> _bestT;
    /// Best C
    Int4Array _bestC;
    /// Best M
    Double3Array _bestM;
    /// Best number of events that has been found
    int _bestZ;
    /// Complete HotStarts from first iteration of each seed
    std::vector<std::vector<HotStart> > _firstCompleteHotStart;
    /// Complete HotStart related to best C
    std::vector<HotStart> _lastCompleteHotStart;
    /// Flags tracking whether the realtive value has already been computed
    std::vector<bool> _isComputed;
    /// Total number of values in F, used to normaize the objective value
    double _norm;
    /// Standard HotStart that can be always used
    HotStart _diploidCompleteHotStart;
    /// Refined tree
    CopyNumberTree _refinedTree;
    /// Objective value after refinement
    double _refinedObjValue;
    
    void runInstance(const int Z, const int seedIdx, const HotStart &inputCompleteHotStart);
    void computeDistance(const int Z);
    void computeDistance(const int Z, const HotStart &inputCompleteHotStart);
    const HotStart& previousCompleteHotStart(const int Z, const int seedIdx);
    inline bool isImproving(const double previous, const double successive) const;
    void initialize();
    void initializeZwithPrevious(const int Z);
    void refinement();
};

#endif // _MANAGER_H_
