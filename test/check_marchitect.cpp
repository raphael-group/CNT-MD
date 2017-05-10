#include "basic_check.h"
#include "marchitect.h"


const InputInstance makeInstance(Double3Array &F);
const InputInstance makeFAllDiploid(unsigned int num_chr, unsigned int num_sam);
const InputInstance makeFInt(unsigned int num_chr, unsigned int num_sam);
const InputInstance makeFDouble(unsigned int num_chr, unsigned int num_sam);

int checkSingle();
const ReturnMessage testSingle(const InputInstance &inst, const unsigned int max_cn);

int checkDouble();
const ReturnMessage testDouble(const InputInstance &inst, const unsigned int max_cn);

int checkFewLeaves();
const ReturnMessage testBetween(const InputInstance &inst, const unsigned int max_cn,
                               const unsigned int lb, const unsigned int ub);


int main(int argc, char** argv)
{
    g_verbosity = VERBOSE_NONE;
    std::cout << "CHECKING MARCHITECT" << std::endl;

    if(checkSingle() == EXIT_FAILURE)
        return EXIT_FAILURE;
    if(checkDouble() == EXIT_FAILURE)
        return EXIT_FAILURE;
    if(checkFewLeaves() == EXIT_FAILURE)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}


const InputInstance makeInstance(Double3Array &F, const unsigned int num_chr,
                                 const unsigned int num_sam)
{
    stringstream input;
    input << "#PARAMS" << std::endl;
    input << num_chr << " #number of chromosomes" << std::endl;
    input << num_sam << " #number of samples" << std::endl;
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        input << 10;
        if(c != num_chr-1)
            input << " ";
    }
    input << " #number of segments for each chromosome" << std::endl;
    input << "#SAMPLES" << std::endl;

    for(unsigned int p = 0; p < num_sam; ++p)
    {
        input << "S" << p << " :";
        for(unsigned int c = 0; c < num_chr; ++c)
        {
            for(unsigned int s = 0; s < 10; ++s)
            {
                input << " " << F[c][p][s];
            }

            if(c != num_chr-1)
            {
                input << " |";
            }
        }
        input << std::endl;
    }

    InputInstance inst;
    input >> inst;

    return inst;
}


const InputInstance makeFAllDiploid(unsigned int num_chr, unsigned int num_sam)
{
    Double3Array F(num_chr,
                   DoubleMatrix(num_sam,
                                DoubleArray(10, 0)));

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        for(unsigned int p = 0; p < num_sam; ++p)
        {
            F[c][p][0] = 2; F[c][p][1] = 2; F[c][p][2] = 2; F[c][p][3] = 2; F[c][p][4] = 2;
            F[c][p][5] = 2; F[c][p][6] = 2; F[c][p][7] = 2; F[c][p][8] = 2; F[c][p][9] = 2;
        }
    }

    return makeInstance(F, num_chr, num_sam);
}


const InputInstance makeFInt(unsigned int num_chr, unsigned int num_sam)
{
    Double3Array F(num_chr,
                   DoubleMatrix(num_sam,
                                DoubleArray(10, 0)));

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        for(unsigned int p = 0; p < num_sam; ++p)
        {
            if(p%4 == 0) {
                F[c][p][0] = 2; F[c][p][1] = 2; F[c][p][2] = 2; F[c][p][3] = 2; F[c][p][4] = 2;
                F[c][p][5] = 2; F[c][p][6] = 2; F[c][p][7] = 2; F[c][p][8] = 2; F[c][p][9] = 2;
            } else if (p%4 == 0) {
                F[c][p][0] = 3; F[c][p][1] = 3; F[c][p][2] = 3; F[c][p][3] = 3; F[c][p][4] = 3;
                F[c][p][5] = 3; F[c][p][6] = 3; F[c][p][7] = 3; F[c][p][8] = 3; F[c][p][9] = 3;
            } else if (p%4 == 0) {
                F[c][p][0] = 4; F[c][p][1] = 4; F[c][p][2] = 4; F[c][p][3] = 4; F[c][p][4] = 4;
                F[c][p][5] = 4; F[c][p][6] = 4; F[c][p][7] = 4; F[c][p][8] = 4; F[c][p][9] = 4;
            } else {
                F[c][p][0] = 5; F[c][p][1] = 5; F[c][p][2] = 5; F[c][p][3] = 5; F[c][p][4] = 5;
                F[c][p][5] = 5; F[c][p][6] = 5; F[c][p][7] = 5; F[c][p][8] = 5; F[c][p][9] = 5;
            }
        }
    }

    return makeInstance(F, num_chr, num_sam);
}


const InputInstance makeFDouble(unsigned int num_chr, unsigned int num_sam)
{
    Double3Array F(num_chr,
                   DoubleMatrix(num_sam,
                                DoubleArray(10, 0)));

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        for(unsigned int p = 0; p < num_sam; ++p)
        {
            if(p%4 == 0) {
                F[c][p][0] = 2.2; F[c][p][1] = 2.2; F[c][p][2] = 2.2; F[c][p][3] = 2.2; F[c][p][4] = 2.2;
                F[c][p][5] = 2.2; F[c][p][6] = 2.2; F[c][p][7] = 2.2; F[c][p][8] = 2.2; F[c][p][9] = 2.2;
            } else if (p%4 == 0) {
                F[c][p][0] = 1.3; F[c][p][1] = 1.3; F[c][p][2] = 1.3; F[c][p][3] = 1.3; F[c][p][4] = 1.3;
                F[c][p][5] = 1.3; F[c][p][6] = 1.3; F[c][p][7] = 1.3; F[c][p][8] = 1.3; F[c][p][9] = 1.3;
            } else if (p%4 == 0) {
                F[c][p][0] = 2.4; F[c][p][1] = 2.4; F[c][p][2] = 2.4; F[c][p][3] = 2.4; F[c][p][4] = 2.4;
                F[c][p][5] = 2.4; F[c][p][6] = 2.4; F[c][p][7] = 2.4; F[c][p][8] = 2.4; F[c][p][9] = 2.4;
            } else {
                F[c][p][0] = 1.5; F[c][p][1] = 1.5; F[c][p][2] = 1.5; F[c][p][3] = 1.5; F[c][p][4] = 1.5;
                F[c][p][5] = 1.5; F[c][p][6] = 1.5; F[c][p][7] = 1.5; F[c][p][8] = 1.5; F[c][p][9] = 1.5;
            }
        }
    }

    return makeInstance(F, num_chr, num_sam);
}

int checkSingle()
{
    {
        ReturnMessage m(testSingle(makeFAllDiploid(1, 2), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFAllDiploid(1, 4), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFAllDiploid(1, 6), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFAllDiploid(1, 8), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFAllDiploid(4, 2), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFAllDiploid(4, 4), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFAllDiploid(4, 6), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFAllDiploid(4, 6), 0));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }


    {
        ReturnMessage m(testSingle(makeFInt(1, 2), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFInt(1, 4), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFInt(1, 6), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFInt(1, 8), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFInt(4, 2), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFInt(4, 4), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFInt(4, 6), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFInt(4, 6), 0));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testSingle(makeFInt(4, 8), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}


const ReturnMessage testSingle(const InputInstance &inst, const unsigned int max_cn)
{
    std::cout << "- Check example with single clone per sample " << inst.numChr()
              << " chromosomes, " << inst.F()[0].size() << " samples : ";

    Int3Array C(inst.F().size(),
                IntMatrix(inst.F()[0].size(),
                          IntArray(inst.F()[0][0].size(), 0)));

    for(unsigned int c = 0; c < C.size(); ++c)
    {
        for(unsigned int p = 0; p < C[0].size(); ++p)
        {
            for(unsigned int s = 0; s < C[0][0].size(); ++s)
            {
                C[c][p][s] = inst.F()[c][p][s];
                if(inst.F()[c][p][s] != C[c][p][s])
                {
                    return ReturnMessage(ReturnType::FAILURE, "Wrong input to the test, F is not integer");
                }
            }
        }
    }

    IntMatrix e(C.size(), IntArray(C[0][0].size(), max_cn));

    MArchitect architect(inst, C, C[0].size(), e);
    architect.init();
    architect.solve(0, 0);

    DoubleMatrix M(architect.getM());

    for(unsigned int i = 0; i < M.size(); ++i)
    {
        double sum = 0.0;
        for(unsigned int j = 0; j < M[i].size(); ++j)
        {
            sum += M[i][j];
        }
        if(sum < 0.999 || sum > 1.001)
            return ReturnMessage(ReturnType::FAILURE, "The sum of the rows is not 1");
    }

    if(architect.getObjValue() > 0.00001)
        return ReturnMessage(ReturnType::NO_ZERO, "The objective function is non-zero");

    for(unsigned int c = 0; c < inst.numChr(); ++c)
    {
        for(unsigned int p = 0; p < inst.m(); ++p)
        {
            for(unsigned int s = 0; s < inst.n()[c]; ++s)
            {
                double res = 0.0;
                for(unsigned int i = 0; i < inst.m(); ++i)
                {
                    res += M[p][i] * C[c][i][s];
                }
                if(res < (inst.F()[c][p][s] - 0.0001) || res > (inst.F()[c][p][s] + 0.0001))
                {
                    return ReturnMessage(ReturnType::FAILURE, "Wrong estimation of F");
                }
            }
        }
    }

    return ReturnMessage(ReturnType::SUCCESS);
}


int checkDouble()
{
    {
        ReturnMessage m(testDouble(makeFAllDiploid(1, 2), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFAllDiploid(1, 4), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFAllDiploid(1, 6), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFAllDiploid(1, 8), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFAllDiploid(4, 2), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFAllDiploid(4, 4), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFAllDiploid(4, 6), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFAllDiploid(4, 6), 0));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }


    {
        ReturnMessage m(testDouble(makeFInt(1, 2), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFInt(1, 4), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFInt(1, 6), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFInt(1, 8), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFInt(4, 2), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFInt(4, 4), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFInt(4, 6), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFInt(4, 6), 0));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testDouble(makeFInt(4, 8), 2));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}


const ReturnMessage testDouble(const InputInstance &inst, const unsigned int max_cn)
{
    std::cout << "- Check example with double clone per sample " << inst.numChr()
              << " chromosomes, " << inst.F()[0].size() << " samples : ";

    unsigned int numChr = inst.F().size();
    unsigned int numSam = inst.F()[0].size();
    unsigned int numLea = numSam * 2;
    unsigned int numSeg = inst.F()[0][0].size();

    Int3Array C(numChr,
                IntMatrix(numLea,
                          IntArray(numSeg, 0)));

    for(unsigned int c = 0; c < numChr; ++c)
    {
        for(unsigned int p = 0; p < numSam; ++p)
        {
            for(unsigned int s = 0; s < numSeg; ++s)
            {
                C[c][p*2][s] = floor(inst.F()[c][p][s]);
                C[c][p*2+1][s] = ceil(inst.F()[c][p][s]);
            }
        }
    }

    IntMatrix e(numChr, IntArray(numSeg, max_cn));

    MArchitect architect(inst, C, numLea, e);
    architect.init();
    architect.solve(0, 0);

    DoubleMatrix M(architect.getM());

    for(unsigned int i = 0; i < M.size(); ++i)
    {
        double sum = 0.0;
        for(unsigned int j = 0; j < M[i].size(); ++j)
        {
            sum += M[i][j];
        }
        if(sum < 0.999 || sum > 1.001)
            return ReturnMessage(ReturnType::FAILURE, "The sum of the rows is not 1");
    }

    if(architect.getObjValue() > 0.00001)
        return ReturnMessage(ReturnType::NO_ZERO, "The objective function is non-zero");

    for(unsigned int c = 0; c < inst.numChr(); ++c)
    {
        for(unsigned int p = 0; p < inst.m(); ++p)
        {
            for(unsigned int s = 0; s < inst.n()[c]; ++s)
            {
                double res = 0.0;
                for(unsigned int i = 0; i < numLea; ++i)
                {
                    res += M[p][i] * C[c][i][s];
                }
                if(res < (inst.F()[c][p][s] - 0.0001) || res > (inst.F()[c][p][s] + 0.0001))
                {
                    return ReturnMessage(ReturnType::FAILURE, "Wrong estimation of F");
                }
            }
        }
    }

    return ReturnMessage(ReturnType::SUCCESS);
}


int checkFewLeaves()
{
    {
        ReturnMessage m(testBetween(makeFDouble(1, 1), 3, 1, 2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 1), 3, 2, 3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 1), 3, 1, 3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 2), 3, 1, 2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 2), 3, 2, 3));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 2), 3, 1, 3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 4), 3, 1, 2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 4), 3, 2, 3));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 4), 3, 1, 3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 6), 3, 1, 2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 6), 3, 2, 3));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 6), 3, 1, 3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 8), 3, 1, 2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 8), 3, 2, 3));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(1, 8), 3, 1, 3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(4, 8), 3, 1, 2));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(4, 8), 3, 2, 3));
        switch (m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(4, 8), 3, 1, 3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testBetween(makeFDouble(4, 8), 0, 1, 3));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }
}


const ReturnMessage testBetween(const InputInstance &inst, const unsigned int max_cn,
                               const unsigned int lb, const unsigned int ub)
{
    std::cout << "- Check between example per sample " << inst.numChr()
              << " chromosomes, " << inst.F()[0].size() << " samples with clones between "
              << lb << " and " << ub << " : ";

    unsigned int numChr = inst.F().size();
    unsigned int numLea = 2;
    unsigned int numSeg = inst.F()[0][0].size();

    Int3Array C(numChr,
                IntMatrix(numLea,
                          IntArray(numSeg, 0)));

    for(unsigned int c = 0; c < numChr; ++c)
    {
        for(unsigned int s = 0; s < numSeg; ++s)
        {
            C[c][0][s] = lb;
            C[c][1][s] = ub;
        }
    }

    IntMatrix e(numChr, IntArray(numSeg, max_cn));

    MArchitect architect(inst, C, numLea, e);
    architect.init();
    architect.solve(0, 0);

    DoubleMatrix M(architect.getM());

    for(unsigned int i = 0; i < M.size(); ++i)
    {
        double sum = 0.0;
        for(unsigned int j = 0; j < M[i].size(); ++j)
        {
            sum += M[i][j];
        }

        if(sum < 0.999 || sum > 1.001)
            return ReturnMessage(ReturnType::FAILURE, "The sum of the rows is not 1");
    }

    if(architect.getObjValue() > 0.00001)
        return ReturnMessage(ReturnType::NO_ZERO, "The objective function is non-zero");

    for(unsigned int c = 0; c < inst.numChr(); ++c)
    {
        for(unsigned int p = 0; p < inst.m(); ++p)
        {
            for(unsigned int s = 0; s < inst.n()[c]; ++s)
            {
                double res = 0.0;
                for(unsigned int i = 0; i < numLea; ++i)
                {
                    res += M[p][i] * C[c][i][s];
                }
                if(res < (inst.F()[c][p][s] - 0.0001) || res > (inst.F()[c][p][s] + 0.0001))
                {
                    return ReturnMessage(ReturnType::FAILURE, "Wrong estimation of F");
                }
            }
        }
    }

    return ReturnMessage(ReturnType::SUCCESS, "The objective value is zero");
}
