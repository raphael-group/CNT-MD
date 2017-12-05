#include <fstream>
#include <ios>
#include <stdlib.h>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include <lemon/arg_parser.h>

#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/gamma.hpp>

#include "basic_types.h"
#include "manager.h"

int main(int argc, char** argv)
{
    int timeLimit = -1;
    int memoryLimit = -1;
    const int size_bubbles = 10;
    std::string outputFilename;

    int maxCopyNumber = -1;
    int maxSizeTree = -1;
    int lbMaxSizeTree = 0;
    int k = 4;
    double eps = 0.0;
    bool forceDiploid = false;
    bool rootNotFixed = false;
    bool deactiveRefinement = false;

    int numStarts = 10;
    int numIterConvergence = 2;
    int maxIter = 7;
    int numILPThreads = 1;
    int numWorkers = 2;
    
    int seed = 0;
    int verbosityLevel = 1;
    int mode = 1;

    lemon::ArgParser ap(argc, argv);
    ap.refOption("o", "Output filename", outputFilename)
      .refOption("Z", "Maximum size of tree for all chromosomes", maxSizeTree, true)
      .refOption("lbZ", "Lower bound for maximum size of tree for all chromosomes (default: 0)", lbMaxSizeTree)
      .refOption("k", "Number of leaves", k, true)
      .refOption("t", "Epsilon, threshold level of tolerance for normalized distance (default: 0.0)", eps)
      .refOption("r", "Mode for searching parsimonious number of events: (1) Binary Search (2) Reverse Iterative (3) Full Iterative (default: 1)", mode)
      .refOption("s", "Time limit in seconds for each C-step (default: -1, disabled)", timeLimit)
      .refOption("ns", "Number of starting seeds (default: 10)", numStarts)
      .refOption("ni", "Number of iterations per seed (default: 7)", maxIter)
      .refOption("j", "Number of workers (default: 2)", numWorkers)
      .refOption("nt", "Number of ILP threads (default: 1)", numILPThreads)
      .refOption("m", "Memory limit in MB for each worker (default: -1, disabled)", memoryLimit)
      .refOption("e", "Maximum copy number (default: -1, inferred from leaves)", maxCopyNumber)
      .refOption("d", "Force one clone to be the normal diploid (default: false)", forceDiploid)
      .refOption("f", "Do not fix root to all 2s", rootNotFixed)
      .refOption("ss", "Random number seed (default: 0)", seed)
      .refOption("v", "Verbosity level from 0 to 4 (default: 1)", verbosityLevel)
      .refOption("dr", "Deactivate refinement", deactiveRefinement)
      .other("input", "Input file");
    ap.parse();
    g_rng = std::mt19937(seed);
    g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

    if (ap.files().size() == 0)
    {
        std::cerr << "ERROR: missing input file" << std::endl;
        return 1;
    }

    InputInstance inputInstance;
    std::ifstream inFile(ap.files()[0].c_str());

    if (!inFile.good())
    {
        std::cerr << "ERROR: could not open '" << ap.files()[0] << "' for reading" << std::endl;
        return EXIT_FAILURE;
    }

    inFile >> inputInstance;

    if (maxCopyNumber < 0)
    {
        maxCopyNumber = inputInstance.e();
        std::cerr << "e:  " << maxCopyNumber << endl;
    }
    IntMatrix e(inputInstance.numChr());
    for(unsigned int chr = 0; chr < inputInstance.numChr(); ++chr)
    {
        e[chr] = IntArray(inputInstance.n()[chr], maxCopyNumber);
    }
    
    Manager manager(inputInstance, k, e, lbMaxSizeTree, maxSizeTree,
                    forceDiploid, rootNotFixed, deactiveRefinement,
                    size_bubbles, numIterConvergence,
                    maxIter, numStarts, numWorkers, numILPThreads,
                    timeLimit, memoryLimit, eps);
    switch(mode)
    {
        case(1): manager.runBinarySearch(); break;
        case(2): manager.runReverse()   ; break;
        case(3): manager.runIterative(); break;
    }

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
        std::cerr << std::endl;
        std::cerr << "SOLUTION FOUND: " << manager.getObjValue() << std::endl;
        std::cerr << std::endl;
    }

    if(outputFilename.empty())
    {
        std::cout << manager.getSolution();
    }
    else
    {
        std::ofstream ofs(outputFilename.c_str());
        if (ofs.good())
        {
            ofs << manager.getSolution();
            ofs.close();
        }
        else
        {
            std::cerr << "Error: could not open '" << outputFilename << "' for writing" << std::endl;
            return EXIT_FAILURE;
        }
    }

    return 0;
}
