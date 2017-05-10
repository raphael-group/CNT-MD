#include "basic_check.h"
#include "manager.h"



int checkAllEqual();
const ReturnMessage testAllEqual(const InputInstance &inst, const unsigned int num_leaves,
                                 const unsigned int max_cn, const unsigned int max_events);

int checkComplete();
const ReturnMessage testComplete(const InputInstance &inst, const unsigned int num_leaves,
                                 const unsigned int max_cn, const unsigned int max_events);
const ReturnMessage testUntil(const InputInstance &inst, const unsigned int num_leaves,
                              const unsigned int max_cn, const unsigned int mode);

int checkFMCSolution(const FMCSolution &sol);



int main(int argc, char** argv)
{
    g_verbosity = VERBOSE_DEBUG;
    g_rng = boost::mt19937(1);

    std::cout << "CHECKING MANAGER" << std::endl;

    if(checkAllEqual() == EXIT_FAILURE)
        return EXIT_FAILURE;
    if(checkComplete() == EXIT_FAILURE)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}


int checkAllEqual()
{
    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,8,10,2.0),4,4,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,8,10,2.0),4,4,2*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,8,10,3.0),4,4,2*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,8,10,3.0),4,2,2*2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }
}


const ReturnMessage testAllEqual(const InputInstance &inst, const unsigned int num_leaves,
                                 const unsigned int max_cn, const unsigned int max_events)
{
    std::cout << "- Check all equal with " << inst.numChr() << " chromosomes, " << inst.m()
              << " samples, " << inst.n()[0] << " segments, " << inst.F()[0][0][0] << " value, "
              << num_leaves << " leaves, " << max_cn << " max cn, " << max_events << " max events : ";

    IntMatrix e(inst.numChr(), IntArray(inst.n()[0], max_cn));

    Manager manager(inst, num_leaves, e, std::max((int)max_events-2,0), max_events, false, false, false, 10, 1, 2, 1, 1, 0, 0, 0, 0.0);
    manager.runIterative();

    double obj = manager.getObjValue();
    FMCSolution sol = manager.getSolution();

    if(checkFMCSolution(sol) == EXIT_FAILURE)
        return ReturnMessage(ReturnType::FAILURE, "The FMC solution has been badly serialized");

    if(g_tol.different(obj, 0.0))
        return ReturnMessage(ReturnType::NO_ZERO, "The objective value is non-zero");

    DoubleMatrix M(sol.getM());
    for(unsigned int i = 0; i < M.size(); ++i)
    {
        double sum = 0.0;
        for(unsigned int j = 0; j < M[i].size(); ++j)
        {
            sum += M[i][j];
        }
        if(g_tol.different(sum, 1.0))
            return ReturnMessage(ReturnType::FAILURE, "The sum of the rows is not 1");
    }

    CopyNumberTree tree(sol.getTree());
    if(tree.cost() > max_events)
        return ReturnMessage(ReturnType::FAILURE, "The number of events in the tree is greater than the maximum");

    return ReturnMessage(ReturnType::SUCCESS);
}


int checkComplete()
{
    {
        ReturnMessage m(testComplete(makeCompleteFracInstance(1,6),4,4,3*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteFracInstance(2,6),3,4,3*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testUntil(makeCompleteIntInstance(1,3),4,4,3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS (UNTIL)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testUntil(makeCompleteFracInstance(2,1),3,3,3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS (UNTIL)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testUntil(makeCompleteIntInstance(1,3),4,4,2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS (UNTIL)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testUntil(makeCompleteFracInstance(2,1),3,3,2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS (UNTIL)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testUntil(makeCompleteIntInstance(1,3),4,4,1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS (UNTIL)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testUntil(makeCompleteFracInstance(2,1),3,3,1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS (UNTIL)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }
}


const ReturnMessage testComplete(const InputInstance &inst, const unsigned int num_leaves,
                                 const unsigned int max_cn, const unsigned int max_events)
{
    std::cout << "- Check complete example with " << inst.numChr() << " chromosomes, " << inst.m()
              << " samples, " << inst.n()[0] << " segments, "
              << num_leaves << " leaves, " << max_cn << " max cn, " << max_events << " max events : ";

    IntMatrix e(inst.numChr(), IntArray(inst.n()[0], max_cn));

    Manager manager(inst, num_leaves, e, std::max((int)max_events-1,0), max_events, false, false, false, 10, 1, 2, 2, 1, 0, 0, 0, 0.0);
    manager.runIterative();

    double obj = manager.getObjValue();
    FMCSolution sol = manager.getSolution();

    if(checkFMCSolution(sol) == EXIT_FAILURE)
        return ReturnMessage(ReturnType::FAILURE, "The FMC solution has been badly serialized");

    Int3Array C(sol.getTree().getLeafProfiles());
    for(unsigned int c = 0; c < inst.numChr(); ++c)
    {
        for(unsigned int i = 0; i < num_leaves; ++i)
        {
            for(unsigned int s = 0; s < inst.n()[c]; ++s)
            {
                if(C[c][i][s] < 0 || C[c][i][s] > max_cn)
                {
                    return ReturnMessage(ReturnType::FAILURE, "An element of a profile is out of bounds");
                }
            }
        }
    }

    DoubleMatrix M(sol.getM());
    for(unsigned int i = 0; i < M.size(); ++i)
    {
        double sum = 0.0;
        for(unsigned int j = 0; j < M[i].size(); ++j)
        {
            sum += M[i][j];
        }
        if(g_tol.different(sum, 1.0))
            return ReturnMessage(ReturnType::FAILURE, "The sum of the rows is not 1");
    }

    CopyNumberTree tree(sol.getTree());
    if(tree.cost() > max_events)
        return ReturnMessage(ReturnType::FAILURE, "The number of events in the tree is greater than the maximum");

    double result = 0.0;
    for(unsigned int c = 0; c < inst.numChr(); ++c)
    {
        for(unsigned int p = 0; p < inst.m(); ++p)
        {
            for(unsigned int s = 0; s < inst.n()[c]; ++s)
            {
                double collect = 0.0;
                for(unsigned int i = 0; i < num_leaves; ++i)
                {
                    collect += M[p][i] * C[c][i][s];
                }
                result += std::abs(inst.F()[c][p][s] - collect);
            }
        }
    }

    if(g_tol.different(result, obj))
        return ReturnMessage(ReturnType::FAILURE, "The distance computed is inconsistent with the extract obj value");

    if(g_tol.different(obj, 0.0))
        return ReturnMessage(ReturnType::NO_ZERO, "The objective value is non-zero");

    return ReturnMessage(ReturnType::SUCCESS);
}


const ReturnMessage testUntil(const InputInstance &inst, const unsigned int num_leaves,
                              const unsigned int max_cn, const unsigned int mode)
{
    string mode_str = "";
    switch(mode)
    {
        case(1): mode_str = "Binary Search mode"; break;
        case(2): mode_str = "Reverse Iterative mode"; break;
        case(3): mode_str = "Iterative mode"; break;
    }
    std::cout << "- Check until complete example on " << mode_str << " mode with " << inst.numChr() << " chromosomes, " << inst.m()
              << " samples, " << inst.n()[0] << " segments, "
              << num_leaves << " leaves, " << max_cn << " max cn : ";

    IntMatrix e(inst.numChr(), IntArray(inst.n()[0], max_cn));
    double prev_obj = inst.numChr() * inst.m() * inst.n()[0] * max_cn;

    unsigned int num_events = inst.n()[0] * inst.numChr();
    while(g_tol.nonZero(prev_obj))
    {
        Manager manager(inst, num_leaves, e, num_events, num_events*2, false, false, false, 10, 1, 3, 5, 1, 0, 0, 0, 0.0);
        switch(mode)
        {
            case(1): manager.runBinarySearch(); break;
            case(2): manager.runReverse(); break;
            case(3): manager.runIterative(); break;
        }

        double obj = manager.getObjValue();
        FMCSolution sol = manager.getSolution();

        if(checkFMCSolution(sol) == EXIT_FAILURE)
            return ReturnMessage(ReturnType::FAILURE, "The FMC solution has been badly serialized");

        Int3Array C(sol.getTree().getLeafProfiles());
        for(unsigned int c = 0; c < inst.numChr(); ++c)
        {
            for(unsigned int i = 0; i < num_leaves; ++i)
            {
                for(unsigned int s = 0; s < inst.n()[c]; ++s)
                {
                    if(C[c][i][s] < 0 || C[c][i][s] > max_cn)
                    {
                        return ReturnMessage(ReturnType::FAILURE, "An element of a profile is out of bounds");
                    }
                }
            }
        }

        DoubleMatrix M(sol.getM());
        for(unsigned int i = 0; i < M.size(); ++i)
        {
            double sum = 0.0;
            for(unsigned int j = 0; j < M[i].size(); ++j)
            {
                sum += M[i][j];
            }
            if(g_tol.different(sum, 1.0))
                return ReturnMessage(ReturnType::FAILURE, "The sum of the rows is not 1");
        }

        CopyNumberTree tree(sol.getTree());
        if(tree.cost() > num_events*2)
            return ReturnMessage(ReturnType::FAILURE, "The number of events in the tree is greater than the maximum");

        double result = 0.0;
        for(unsigned int c = 0; c < inst.numChr(); ++c)
        {
            for(unsigned int p = 0; p < inst.m(); ++p)
            {
                for(unsigned int s = 0; s < inst.n()[c]; ++s)
                {
                    double collect = 0.0;
                    for(unsigned int i = 0; i < num_leaves; ++i)
                    {
                        collect += M[p][i] * C[c][i][s];
                    }
                    result += std::abs(inst.F()[c][p][s] - collect);
                }
            }
        }

        if(g_tol.different(result, obj))
            return ReturnMessage(ReturnType::FAILURE, "The distance computed is inconsistent with the extract obj value");

        prev_obj = obj;
        num_events = num_events * 4;
    }

    if(g_tol.nonZero(prev_obj))
        return ReturnMessage(ReturnType::NO_ZERO, "The objective value is non-zero");

    return ReturnMessage(ReturnType::SUCCESS);
}


int checkFMCSolution(const FMCSolution &sol)
{


    return EXIT_SUCCESS;
}


