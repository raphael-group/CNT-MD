#include "basic_check.h"
#include "worker.h"



int checkAllEqual();
const ReturnMessage testAllEqual(const InputInstance &inst, const unsigned int num_leaves,
                                 const unsigned int max_cn, const unsigned int max_events);

int checkComplete();
const ReturnMessage testComplete(const InputInstance &inst, const unsigned int num_leaves,
                                 const unsigned int max_cn, const unsigned int max_events);



int main(int argc, char** argv)
{
    g_verbosity = VERBOSE_NONE;
    std::cout << "CHECKING WORKER" << std::endl;

    if(checkAllEqual() == EXIT_FAILURE)
        return EXIT_FAILURE;
    if(checkComplete() == EXIT_FAILURE)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}


int checkAllEqual()
{
    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(1,1,10,2.0),2,4,0*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(1,4,10,2.0),2,4,0*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(1,6,10,2.0),2,4,0*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,1,10,2.0),2,4,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,4,10,2.0),2,4,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),2,4,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),2,4,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),3,4,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),4,4,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),6,4,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),6,3,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),6,2,0*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),6,2,5*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,2.0),6,2,10*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,3.0),6,4,10*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,4.0),6,4,10*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,3.0),6,2,10*2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testAllEqual(makeAllEqualInstance(2,6,10,4.0),6,2,10*2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
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
              << " samples, " << inst.n()[0] << " segments, " << inst.F()[0][0][0] << " value,"
              << num_leaves << " leaves, " << max_cn << " max cn, " << max_events << " max events : ";

    DoubleMatrix M0(inst.m(), DoubleArray(num_leaves, 0.0));
    for(unsigned int p = 0; p < inst.m(); ++p)
    {
        M0[p][p%num_leaves] = 1.0;
    }

    IntMatrix e(inst.numChr(), IntArray(inst.n()[0], max_cn));

    Worker worker(inst, num_leaves, e, max_events, false, false, 1, 2, 0, 0, 0, M0, 0);
    double obj = worker.solve();

    if(g_tol.different(obj, 0.0))
        return ReturnMessage(ReturnType::NO_ZERO, "The objective value is non-zero");

    Int3Array C(worker.getC());
    for(unsigned int c = 0; c < inst.numChr(); ++c)
    {
        for(unsigned int i = 0; i < num_leaves; ++i)
        {
            for(unsigned int s = 0; s < inst.n()[c]; ++s)
            {
                if(g_tol.different(C[c][i][s], inst.F()[0][0][0]))
                {
                    return ReturnMessage(ReturnType::FAILURE, "Found a different value than the others expected");
                }
            }
        }
    }

    DoubleMatrix M(worker.getM());
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

    CopyNumberTree tree(worker.getT());
    if(tree.cost() > max_events)
        return ReturnMessage(ReturnType::FAILURE, "The number of events in the tree is greater than the maximum");

    return ReturnMessage(ReturnType::SUCCESS);
}


int checkComplete()
{
    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(1,2),4,4,3*1));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(1,2),4,4,4*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(1,4),4,4,3*1));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(1,4),4,4,8*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(1,8),4,4,3*1));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(1,8),4,4,8*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(2,2),4,4,3*2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(2,2),4,4,4*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(2,4),4,4,3*2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteIntInstance(2,4),4,4,8*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testComplete(makeCompleteFracInstance(1,4),4,4,8*1));
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
        ReturnMessage m(testComplete(makeCompleteFracInstance(2,4),4,4,3*2));
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
        ReturnMessage m(testComplete(makeCompleteFracInstance(2,4),4,4,8*2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO)" << std::endl; break;
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

    DoubleMatrix M0(inst.m(), DoubleArray(num_leaves, 0.0));
    for(unsigned int p = 0; p < inst.m(); ++p)
    {
        M0[p][p%num_leaves] = 1.0;
    }

    IntMatrix e(inst.numChr(), IntArray(inst.n()[0], max_cn));
    double prev_obj = inst.numChr() * inst.m() * inst.n()[0] * max_cn;

    for(unsigned int num_events = max_events/2; num_events <= max_events; ++num_events)
    {
        Worker worker(inst, num_leaves, e, num_events, false, false, 1, 2, 0, 0, 0, M0, 0);
        double obj = worker.solve();

        Int3Array C(worker.getC());
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

        DoubleMatrix M(worker.getM());
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

        CopyNumberTree tree(worker.getT());
        if(tree.cost() > num_events)
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

        if(!(obj <= prev_obj))
            return ReturnMessage(ReturnType::FAILURE, "The objective values are not monotically non-increasing for increasing number of events");

        prev_obj = obj;
    }

    if(g_tol.different(prev_obj, 0.0))
        return ReturnMessage(ReturnType::NO_ZERO, "The objective value is non-zero");

    return ReturnMessage(ReturnType::SUCCESS);
}
