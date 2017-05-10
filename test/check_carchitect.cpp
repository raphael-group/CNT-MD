#include "basic_check.h"
#include "carchitect.h"

int checkAllDiploid();
const ReturnMessage testAllDiploid(const int num_chr, const int num_sam, const int num_seg,
                                   const int num_leaves, const int max_cn, const int max_e);
int checkFullExampleInd4Clones();
const ReturnMessage testFullExampleInd4Clones(const int num_chr, const int num_sam,
                                              const int max_cn, const int max_e);

int checkFullExampleHalf4Clones();
const ReturnMessage testFullExampleHalf4Clones(const int num_chr, const int max_cn,
                                               const int max_e);


int main(int argc, char** argv)
{
    g_verbosity = VERBOSE_NONE;
    std::cout << "CHECKING CARCHITECT" << std::endl;

    if(checkAllDiploid() == EXIT_FAILURE)
        return EXIT_FAILURE;
    if(checkFullExampleInd4Clones() == EXIT_FAILURE)
        return EXIT_FAILURE;
    if(checkFullExampleHalf4Clones() == EXIT_FAILURE)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}


int checkAllDiploid()
{
    //All Diploid

    //1 chromosome, 1 sample, 5 segments, 4 leaves, 2 max copy-number, 0 max events
    {
        ReturnMessage m(testAllDiploid(1,1,5,4,2,0*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    //1 chromosome, 4 sample, 5 segments, 4 leaves, 2 max copy-number, 0 max events
    {
        ReturnMessage m(testAllDiploid(1,4,5,4,2,0*1));
        switch (m.type)
        {
             case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
             default:
                 std::cout << "FAILED" << std::endl << m.message << std::endl;
                 return EXIT_FAILURE;
        }
    }

    //1 chromosome, 8 sample, 5 segments, 4 leaves, 2 max copy-number, 0 max events
    {
        ReturnMessage m(testAllDiploid(1,8,5,4,2,0*1));
        switch (m.type)
        {
           case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
           default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    //1 chromosome, 1 sample, 5 segments, 4 leaves, 4 max copy-number, 0 max events
    {
        ReturnMessage m(testAllDiploid(1,1,5,4,4,0*1));
        switch (m.type)
        {
           case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
           default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    //1 chromosome, 1 sample, 5 segments, 4 leaves, 2 max copy-number, 10 max events
    {
        ReturnMessage m(testAllDiploid(1,1,5,4,2,10*1));
        switch (m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    //4 chromosome, 1 sample, 5 segments, 4 leaves, 2 max copy-number, 0 max events
    {
        ReturnMessage m(testAllDiploid(4,1,5,4,2,0*4));
        switch (m.type)
        {
           case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
           default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    //4 chromosome, 4 sample, 5 segments, 4 leaves, 2 max copy-number, 0 max events
    {
        ReturnMessage m(testAllDiploid(4,4,5,4,2,0*4));
        switch (m.type)
        {
           case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
           default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    //4 chromosome, 8 sample, 5 segments, 4 leaves, 2 max copy-number, 0 max events
    {
        ReturnMessage m(testAllDiploid(4,8,5,4,2,0*4));
        switch (m.type)
        {
           case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
           default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    //4 chromosome, 1 sample, 5 segments, 4 leaves, 4 max copy-number, 0 max events
    {
        ReturnMessage m(testAllDiploid(4,1,5,4,4,0*4));
        switch (m.type)
        {
           case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
           default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    //2 chromosome, 1 sample, 5 segments, 4 leaves, 2 max copy-number, 10 max events
    {
        ReturnMessage m(testAllDiploid(2,1,5,4,2,10*2));
        switch (m.type)
        {
           case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
           default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }
}


const ReturnMessage testAllDiploid(const int num_chr, const int num_sam, const int num_seg,
                                   const int num_leaves, const int max_cn, const int max_e)
{
    std::cout << "- Check All Diploid with " << num_chr << " chromosomes, " << num_sam
              << " samples, " << num_leaves << " leaves, max cn of " << max_cn
              << ", and max events of " << max_e << " : ";

    stringstream input;
    input << "#PARAMS" << std::endl;
    input << num_chr << " #number of chromosomes" << std::endl;
    input << num_sam << " #number of samples" << std::endl;
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        input << num_seg;
        if(c != num_chr-1)
                input << " ";
    }
    input << " #number of segments for each chromosome" << std::endl;
    input << "#SAMPLES" << std::endl;
    for(unsigned int i = 0; i < num_sam; ++i)
    {
        input << "S" << i+1 << ":";
        for(unsigned int c = 0; c < num_chr; ++c)
        {
            for(unsigned int s = 0; s < num_seg; ++s)
            {
                input << " 2.0";
            }

            if (c != num_chr-1)
            {
                input << " |";
            }
        }
        input << std::endl;
    }

    InputInstance instance;
    input >> instance;

    DoubleMatrix M(num_sam, DoubleArray(num_leaves, 0.0));
    for(unsigned int p = 0; p < num_sam; ++p)
    {
        M[p][p%num_leaves] = 1.0;
    }

    /**
    std::cout << std::endl;
    for(unsigned int i = 0; i < M.size(); ++i)
    {
        for(unsigned int j = 0; j < M[i].size(); ++j)
        {
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
    **/

    IntMatrix MAX_CN(num_chr, IntArray(num_seg, max_cn));
    int MAX_E = max_e;

    CArchitect architect(instance, M, MAX_CN, MAX_E, num_leaves, false, false);
    architect.init();
    architect.solve(0,0,0);

    if(architect.getObjValue() > 0.00001)
        return ReturnMessage(ReturnType::NO_ZERO, "No zero objective value solution");

    Int3Array C(architect.getC());
    for(unsigned int c = 0; c < C.size(); ++c)
    {
        for(unsigned int i = 0; i < C[c].size(); ++i)
        {
            for(unsigned int j = 0; j < C[c][i].size(); ++j)
            {
                if(C[c][i][j] != 2 && i < num_sam)
                    return ReturnMessage(ReturnType::FAILURE, "A fixed profile containes an copy-number different than 2");
            }
        }
    }

    if(architect.getDelta() > max_e)
        return ReturnMessage(ReturnType::FAILURE, "The number of events is greater than the imposed upper bound");

    return ReturnMessage(ReturnType::SUCCESS);
}


int checkFullExampleInd4Clones()
{
    //Full Example with 4 Clones

    {
        ReturnMessage m(testFullExampleInd4Clones(1, 4, 4, 8*1));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 4, 4, 8*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(3, 4, 4, 8*3));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 6, 4, 8*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 8, 4, 8*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 12, 4, 8*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 4, 2, 8*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 4, 3, 8*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 6, 4, 0*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 6, 4, 4*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 6, 4, 7*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 6, 4, 9*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            case(ReturnType::FAILURE):
                if(m.message.substr(m.message.length() - 5) == "wrong")
                    std::cout << "SUCCESS (DIFFERENT INTERNALS)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleInd4Clones(2, 6, 4, 12*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            case(ReturnType::FAILURE):
                if(m.message.substr(m.message.length() - 5) == "wrong")
                    std::cout << "SUCCESS (DIFFERENT INTERNALS)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }
}


const ReturnMessage testFullExampleInd4Clones(const int num_chr, const int num_sam,
                                              const int max_cn, const int max_e)
{
    std::cout << "- Check Full Example of 4 Independent Clones with " << num_chr << " chromosomes, " << num_sam
              << " samples, " << "max cn of " << max_cn
              << ", and max events of " << max_e << " : ";

    //Build the input as InputInstance format
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

    //Build samples
    vector<string> samples;
    vector<string> test(7, "");
    string temp = "";

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        test[0] += "2 2 2 2 2 2 2 2 2 2";
        test[1] += "0 0 0 0 0 2 2 2 2 2";
        test[2] += "2 2 2 2 2 0 0 0 0 0";

        if(c < num_chr - 1)
        {
            test[0] += " | ";
            test[1] += " | ";
            test[2] += " | ";
        }
    }

    temp = "S1 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "0 0 0 0 0 2.0 2.0 2.0 3.0 3.0";
        test[3] += "0 0 0 0 0 2 2 2 3 3";
        if(c < num_chr - 1)
        {
            temp += " | ";
            test[3] += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S2 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "0 0 0 0 0 2.0 1.0 2.0 2.0 2.0";
        test[4] += "0 0 0 0 0 2 1 2 2 2";
        if(c < num_chr - 1)
        {
            temp += " | ";
            test[4] += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S3 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "4.0 4.0 4.0 4.0 4.0 0 0 0 0 0";
        test[5] += "4 4 4 4 4 0 0 0 0 0";
        if(c < num_chr - 1)
        {
            temp += " | ";
            test[5] += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S4 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "2.0 2.0 2.0 2.0 2.0 0 0 0 0 0";
        test[6] += "2 2 2 2 2 0 0 0 0 0";
        if(c < num_chr - 1)
        {
            temp += " | ";
            test[6] += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    input << "#SAMPLES" << std::endl;
    for(unsigned int i = 0; i < num_sam; ++i)
    {
        input << samples[i%4] << std::endl;
    }

    InputInstance instance;
    input >> instance;

    DoubleMatrix M(num_sam, DoubleArray(4, 0.0));
    for(unsigned int p = 0; p < num_sam; ++p)
    {
        M[p][p%4] = 1.0;
    }

    IntMatrix MAX_CN(num_chr, IntArray(10, max_cn));
    int MAX_E = max_e;

    CArchitect architect(instance, M, MAX_CN, MAX_E, 4, false, false);
    architect.init();
    architect.solve(0,0,0);

    if(architect.getObjValue() > 0.00001)
        return ReturnMessage(ReturnType::NO_ZERO, "No zero objective value solution");

    stringstream out;
    out << architect.getTree() << std::endl;
    string line;

    do {
        getline(out, line, '\n');
    } while(line != "#PROFILES");

    vector<string> profiles;
    do {
        getline(out, line, '\n');
        if(line != "#EDGES")
        {
            line = line.substr(line.find(':') + 2);
            profiles.push_back(line);
        }
    } while(line != "#EDGES");

    if(profiles[0] != test[0])
        return ReturnMessage(ReturnType::FAILURE, "The inferred root is not the diploid");

    if(!((profiles[1] == test[1] && profiles[2] == test[2]) ||
         (profiles[1] == test[2] && profiles[2] == test[1])))
    {
        return ReturnMessage(ReturnType::FAILURE, "The inferred internal vertices are wrong");
    }

    if(profiles[3] != test[3] && profiles[3] != test[4] &&
       profiles[3] != test[5] && profiles[3] != test[6])
        return ReturnMessage(ReturnType::FAILURE, "The inferred clone profile 3 is wrong");

    if(profiles[4] != test[3] && profiles[4] != test[4] &&
       profiles[4] != test[5] && profiles[4] != test[6])
        return ReturnMessage(ReturnType::FAILURE, "The inferred clone profile 4 is wrong");

    if(profiles[5] != test[3] && profiles[5] != test[4] &&
       profiles[5] != test[5] && profiles[5] != test[6])
        return ReturnMessage(ReturnType::FAILURE, "The inferred clone profile 5 is wrong");

    if(profiles[6] != test[3] && profiles[6] != test[4] &&
       profiles[6] != test[5] && profiles[6] != test[6])
        return ReturnMessage(ReturnType::FAILURE, "The inferred clone profile 6 is wrong");

    if(architect.getDelta() > max_e)
        return ReturnMessage(ReturnType::FAILURE, "The number of events is greater than the imposed upper bound");

    return ReturnMessage(ReturnType::SUCCESS);
}


int checkFullExampleHalf4Clones()
{
    {
        ReturnMessage m(testFullExampleHalf4Clones(1, 4, 8*1));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleHalf4Clones(2, 4, 8*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleHalf4Clones(3, 4, 8*3));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleHalf4Clones(2, 3, 8*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleHalf4Clones(2, 2, 8*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleHalf4Clones(2, 4, 7*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleHalf4Clones(2, 4, 0*2));
        switch(m.type)
        {
            case(ReturnType::NO_ZERO): std::cout << "SUCCESS (NON-ZERO OBJ)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleHalf4Clones(2, 4, 9*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            case(ReturnType::FAILURE):
                if(m.message.substr(m.message.length() - 5) == "wrong")
                    std::cout << "SUCCESS (DIFFERENT INTERNALS)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }

    {
        ReturnMessage m(testFullExampleHalf4Clones(2, 4, 20*2));
        switch(m.type)
        {
            case(ReturnType::SUCCESS): std::cout << "SUCCESS" << std::endl; break;
            case(ReturnType::FAILURE):
                if(m.message.substr(m.message.length() - 5) == "wrong")
                    std::cout << "SUCCESS (DIFFERENT INTERNALS)" << std::endl; break;
            default:
                std::cout << "FAILED" << std::endl << m.message << std::endl;
                return EXIT_FAILURE;
        }
    }
}


const ReturnMessage testFullExampleHalf4Clones(const int num_chr, const int max_cn,
                                               const int max_e)
{
    std::cout << "- Check Full Example of 4 Half Clones with " << num_chr
              << " chromosomes, " << "max cn of " << max_cn
              << ", and max events of " << max_e << " : ";

    //Build the input as InputInstance format
    stringstream input;
    input << "#PARAMS" << std::endl;
    input << num_chr << " #number of chromosomes" << std::endl;
    input << "8 #number of samples" << std::endl;
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        input << 10;
        if(c != num_chr-1)
                input << " ";
    }
    input << " #number of segments for each chromosome" << std::endl;

    //Build samples
    vector<string> samples;
    vector<string> test(7, "");

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        test[0] += "2 2 2 2 2 2 2 2 2 2";
        test[1] += "0 0 0 0 0 2 2 2 2 2";
        test[2] += "2 2 2 2 2 0 0 0 0 0";

        if(c < num_chr - 1)
        {
            test[0] += " | ";
            test[1] += " | ";
            test[2] += " | ";
        }
    }

    string temp1 = "S1 : ";
    string temp2 = "S2 : ";
    string temp3 = "S3 : ";
    string temp4 = "S4 : ";
    string temp5 = "S5 : ";
    string temp6 = "S6 : ";
    string temp7 = "S7 : ";
    string temp8 = "S8 : ";
    string temp9 = "S9 : ";
    string temp10 = "S10 : ";

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp1 += "0.0 0.0 0.0 0.0 0.0 2.0 1.1 2.0 2.1 2.1";
        temp2 += "0.0 0.0 0.0 0.0 0.0 2.0 1.6 2.0 2.6 2.6";
        temp3 += "0.8 0.8 0.8 0.8 0.8 1.4 1.4 1.4 2.1 2.1";
        temp4 += "0.4 0.4 0.4 0.4 0.4 1.6 1.2 1.6 2.0 2.0";
        temp5 += "3.2 3.2 3.2 3.2 3.2 0.4 0.2 0.4 0.4 0.4";
        temp6 += "2.0 2.0 2.0 2.0 2.0 1.0 0.5 1.0 1.0 1.0";
        temp7 += "1.0 1.0 1.0 1.0 1.0 1.0 0.5 1.0 1.0 1.0";
        temp8 += "3.0 3.0 3.0 3.0 3.0 0.0 0.0 0.0 0.0 0.0";
        temp9 += "0.0 0.0 0.0 0.0 0.0 2.0 2.0 2.0 3.0 3.0";
        temp10 += "2.0 2.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0";
        if(c < num_chr - 1)
        {
            temp1 += " | ";
            temp2 += " | ";
            temp3 += " | ";
            temp4 += " | ";
            temp5 += " | ";
            temp6 += " | ";
            temp7 += " | ";
            temp8 += " | ";
            temp9 += " | ";
            temp10 += " | ";
        }
    }
    samples.push_back(temp1);
    samples.push_back(temp2);
    samples.push_back(temp3);
    samples.push_back(temp4);
    samples.push_back(temp5);
    samples.push_back(temp6);
    samples.push_back(temp7);
    samples.push_back(temp8);
    samples.push_back(temp9);
    samples.push_back(temp10);

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        test[3] += "0 0 0 0 0 2 2 2 3 3";
        if(c < num_chr - 1)
        {
            test[3] += " | ";
        }
    }

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        test[4] += "0 0 0 0 0 2 1 2 2 2";
        if(c < num_chr - 1)
        {
            test[4] += " | ";
        }
    }

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        test[5] += "4 4 4 4 4 0 0 0 0 0";
        if(c < num_chr - 1)
        {
            test[5] += " | ";
        }
    }

    for(unsigned int c = 0; c < num_chr; ++c)
    {
        test[6] += "2 2 2 2 2 0 0 0 0 0";
        if(c < num_chr - 1)
        {
            test[6] += " | ";
        }
    }

    input << "#SAMPLES" << std::endl;
    for(unsigned int i = 0; i < 8; ++i)
    {
        input << samples[i] << std::endl;
    }

    InputInstance instance;
    input >> instance;

    DoubleMatrix M(10, DoubleArray(4, 0.0));
    M[0][0] = 0.1; M[0][1] = 0.9; M[0][2] = 0; M[0][3] = 0;
    M[1][0] = 0.6; M[1][1] = 0.4; M[1][2] = 0; M[1][3] = 0;
    M[2][0] = 0.7; M[2][1] = 0; M[2][2] = 0.1; M[2][3] = 0.2;
    M[3][0] = 0.4; M[3][1] = 0.4; M[3][2] = 0; M[3][3] = 0.2;
    M[4][0] = 0; M[4][1] = 0.2; M[4][2] = 0.8; M[4][3] = 0;
    M[5][0] = 0; M[5][1] = 0.5; M[5][2] = 0.5; M[5][3] = 0;
    M[6][0] = 0; M[6][1] = 0.5; M[6][2] = 0; M[6][3] = 0.5;
    M[7][0] = 0; M[7][1] = 0; M[7][2] = 0.5; M[7][3] = 0.5;
    M[8][0] = 1.0; M[8][1] = 0; M[8][2] = 0; M[8][3] = 0;
    M[9][0] = 0; M[9][1] = 0; M[9][2] = 0; M[9][3] = 1.0;

    IntMatrix MAX_CN(num_chr, IntArray(10, max_cn));
    int MAX_E = max_e;

    CArchitect architect(instance, M, MAX_CN, MAX_E, 4, false, false);
    architect.init();
    architect.solve(0,0,0);

    if(architect.getObjValue() > 0.00001)
        return ReturnMessage(ReturnType::NO_ZERO, "No zero objective value solution");

    stringstream out;
    out << architect.getTree() << std::endl;
    string line;

    do {
        getline(out, line, '\n');
    } while(line != "#PROFILES");

    vector<string> profiles;
    do {
        getline(out, line, '\n');
        if(line != "#EDGES")
        {
            line = line.substr(line.find(':') + 2);
            profiles.push_back(line);
        }
    } while(line != "#EDGES");

    if(profiles[0] != test[0])
        return ReturnMessage(ReturnType::FAILURE, "The inferred root is not the diploid");

    if(!((profiles[1] == test[1] && profiles[2] == test[2]) ||
         (profiles[1] == test[2] && profiles[2] == test[1])))
    {
        return ReturnMessage(ReturnType::FAILURE, "The inferred internal vertices are wrong");
    }

    if(profiles[3] != test[3] && profiles[3] != test[4] &&
       profiles[3] != test[5] && profiles[3] != test[6])
        return ReturnMessage(ReturnType::FAILURE, "The inferred clone profile 3 is wrong");

    if(profiles[4] != test[3] && profiles[4] != test[4] &&
       profiles[4] != test[5] && profiles[4] != test[6])
        return ReturnMessage(ReturnType::FAILURE, "The inferred clone profile 4 is wrong");

    if(profiles[5] != test[3] && profiles[5] != test[4] &&
       profiles[5] != test[5] && profiles[5] != test[6])
        return ReturnMessage(ReturnType::FAILURE, "The inferred clone profile 5 is wrong");

    if(profiles[6] != test[3] && profiles[6] != test[4] &&
       profiles[6] != test[5] && profiles[6] != test[6])
        return ReturnMessage(ReturnType::FAILURE, "The inferred clone profile 6 is wrong");

    if(architect.getDelta() > max_e)
        return ReturnMessage(ReturnType::FAILURE, "The number of events is greater than the imposed upper bound");

    return ReturnMessage(ReturnType::SUCCESS);
}

