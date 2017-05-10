#include "basic_check.h"

const InputInstance makeAllEqualInstance(const unsigned int num_chr, const unsigned int num_sam,
                                        const unsigned int num_seg, const double value)
{
    std::stringstream input;
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
                input << " " << value;
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

    return instance;
}


const InputInstance makeCompleteIntInstance(const unsigned int num_chr, const unsigned int num_sam)
{
    std::stringstream input;
    input << "#PARAMS" << std::endl;
    input << num_chr << " #number of chromosomes" << std::endl;
    input << num_sam << " #number of samples" << std::endl;
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        input << 6;
        if(c != num_chr-1)
            input << " ";
    }
    input << " #number of segments for each chromosome" << std::endl;

    //Build samples
    std::vector<std::string> samples;
    std::string temp = "";

    temp = "S1 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "0 2.0 2.0 2.0 3.0 3.0";
        if(c < num_chr - 1)
        {
            temp += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S2 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "0 2.0 1.0 2.0 2.0 2.0";
        if(c < num_chr - 1)
        {
            temp += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S3 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "4.0 0 0 0 0 0";
        if(c < num_chr - 1)
        {
            temp += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S4 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "2.0 0 0 0 0 0";
        if(c < num_chr - 1)
        {
            temp += " | ";
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

    return instance;
}


const InputInstance makeCompleteFracInstance(const unsigned int num_chr, const unsigned int num_sam)
{
    std::stringstream input;
    input << "#PARAMS" << std::endl;
    input << num_chr << " #number of chromosomes" << std::endl;
    input << num_sam << " #number of samples" << std::endl;
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        input << 6;
        if(c != num_chr-1)
            input << " ";
    }
    input << " #number of segments for each chromosome" << std::endl;

    //Build samples
    std::vector<std::string> samples;
    std::string temp = "";

    temp = "S1 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "0.3 2.2 2.2 2.8 2.2 2.1";
        if(c < num_chr - 1)
        {
            temp += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S2 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "0.9 2.2 1.1 2.2 2.3 2.2";
        if(c < num_chr - 1)
        {
            temp += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S3 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "2.0 0.6 0.2 0.1 0.17 0.189";
        if(c < num_chr - 1)
        {
            temp += " | ";
        }
    }
    samples.push_back(temp);
    temp = "";

    temp = "S4 : ";
    for(unsigned int c = 0; c < num_chr; ++c)
    {
        temp += "2.17 0.18 0.11 0.25 0.11 0.11";
        if(c < num_chr - 1)
        {
            temp += " | ";
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

    return instance;
}
