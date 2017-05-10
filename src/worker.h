#ifndef _WORKER_H_
#define _WORKER_H_

#include <limits>
#include <fstream>
#include <ios>
#include <stdlib.h>
#include <cstdlib>
#include <algorithm>

#include "basic_types.h"
#include "carchitect.h"
#include "marchitect.h"
#include "inputinstance.h"

class Worker
{
public:
    Worker(const InputInstance &inputInstance,
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
           const HotStart& inputCompleteHotStart);

    Worker(const InputInstance &inputInstance,
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
           const int seedIndex);
    
    /// Solve for given seed M0
    double solve();
    
    const Int3Array& getC() const
    {
        return _allC.back();
    }
    
    const DoubleMatrix& getM() const
    {
        return _allM.back();
    }
    
    const CopyNumberTree& getT() const
    {
        return _allTrees.back();
    }

    const HotStart& getFirstCompleteHotStart()
    {
        return _firstCompleteHotStart;
    }

    const HotStart& getLastCompleteHotStart()
    {
        return _lastCompleteHotStart;
    }

private:
    /// Input instance
    const InputInstance& _inputInstance;
    /// Number of leaves
    const unsigned int _k;
    /// Maximum copy number per chromosome, per position
    const IntMatrix& _e;
    /// Maximum number of events per all chromosomes
    const int _Z;
    /// Force the presence of the normal diploid clone
    const bool _forceDiploid;
    /// Do not fix the root to the normal diploid
    const bool _rootNotFixed;
    /// Number of iterations for checking convergence
    const unsigned int _iterConvergence;
    /// Maximum number of iterations for each seed
    const unsigned int _maxIter;
    /// Time limit (seconds)
    const int _timeLimit;
    /// Memory limit (MB)
    const int _memoryLimit;
    /// Number of threads
    const int _nrThreads;
    /// Initial M0
    const DoubleMatrix& _M0;
    /// All C matrices
    std::vector<Int3Array> _allC;
    /// All M matrices
    std::vector<DoubleMatrix> _allM;
    /// All objective values for the C-steps
    DoubleArray _allObjC;
    /// All objective values for the M-steps
    DoubleArray _allObjM;
    /// All the trees
    std::vector<CopyNumberTree> _allTrees;
    /// Seed index
    const int _seedIndex;
    ///The complete HotStart that is used to hotstart the process
    const HotStart _inputCompleteHotStart;
    ///The complete HotStart that is computed after the first iteration
    HotStart _firstCompleteHotStart;
    ///The last complete HotStart
    HotStart _lastCompleteHotStart;
};


#endif //_WORKER_H_
