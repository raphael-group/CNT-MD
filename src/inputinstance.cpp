#include "inputinstance.h"
#include <sstream>
#include <string>
#include <stdexcept>

InputInstance::InputInstance()
    : _F()
    , _m(-1)
    , _num_chr(-1)
    , _n()
{
}

int InputInstance::e() const
{
    int max_e = 0;
    for (int chr = 0; chr < _num_chr; ++chr)
    {
        for (int p = 0; p < _m; ++p)
        {
            for (int s = 0; s < _n[chr]; ++s)
            {
                if (_F[chr][p][s] > max_e)
                {
                    max_e = ceil(_F[chr][p][s]);
                }
            }
        }
    }
    return max_e;
}

IntArray InputInstance::z() const
{
    IntArray result(_num_chr);
    for(int chr = 0; chr < _num_chr; ++chr)
    {
        result[chr] = 0;
        for(int s = 0; s < _n[chr]+1; ++s)
        {
            double contribution_neg = 0;
            double contribution_pos = 0;
            for(int p = 0; p < _m; ++p)
            {
                double left = 0;
                double right = 0;
                
                if(s == 0)
                {
                    left = 2.0;
                    right = _F[chr][p][s];
                } else if(s == _n[chr])
                {
                    left = _F[chr][p][s-1];
                    right = 2.0;
                } else
                {
                    left = _F[chr][p][s-1];
                    right = _F[chr][p][s];
                }
                
                if(right - left >= 0)
                    contribution_pos = std::max(contribution_pos, ceil(fabs(right - left)));
                else
                    contribution_neg = std::max(contribution_neg, ceil(fabs(right - left)));
            }
            result[chr] += contribution_pos + contribution_neg;
        }
        result[chr] = ceil(result[chr] / 2);
    }
    return result;
}

std::ostream& operator<<(std::ostream& out, const InputInstance& instance)
{
    out << "#PARAMS" << std::endl;
    out << instance._num_chr << " #number of chromosomes" << std::endl;
    out << instance._m << " #number of samples" << std::endl;
    for (int n : instance._n)
    {
        out << n << " ";
    }
    out << "#number of segments per chromosome" << std::endl;
    out << "#SAMPLES" << std::endl;

    for (int p = 0; p < instance._m; ++p)
    {
        out << p << " :";
        for (int chr = 0; chr < instance._num_chr; ++chr)
        {
            if (chr != 0)
            {
                out << " |";
            }
            for (int s = 0; s < instance._n[chr]; ++s)
            {
                out << " " << instance._F[chr][p][s];
            }
        }
        out << std::endl;
    }
    return out;
}

std::istream& operator>>(std::istream& in, InputInstance& instance)
{
    std::string line;
    std::string value;
    
    instance._F.clear();
    instance._n.clear();
    
    std::getline(in, line, '\n'); //Skip the first line "#PARAMS"
    
    /** Read the number of chromosomes **/
    std::getline(in, line, '\n');
    std::stringstream chromosomeline(line);
    chromosomeline >> instance._num_chr;
    
    /** Read the number of samples **/
    std::getline(in, line, '\n');
    std::stringstream sampleline(line);
    sampleline >> instance._m;
    
    /** Read the number of segments for each chromosome  **/
    std::getline(in, line, '\n');
    std::stringstream segline(line);
    std::getline(segline, line, '#');
    std::stringstream split_seg(line);
    while (!split_seg.eof())
    {
        std::getline(split_seg, value, ' ');
        if(!value.empty())
        {
            instance._n.push_back((unsigned int)atoi(value.c_str()));
        }
    }
    
    if(instance._n.size() != instance._num_chr)
    {
        throw std::runtime_error(std::string("ERROR: inconsistent numbers of segments in #PARAMS"));
    }
    
    /** Allocate the cube containing the result **/
    instance._F = Double3Array(instance._num_chr, DoubleMatrix(instance._m));
    
    for(unsigned int chr = 0; chr < instance._num_chr; ++chr)
    {
        for(unsigned int sample = 0; sample < instance._m; ++sample)
        {
            instance._F[chr][sample] = DoubleArray(instance._n[chr]);
        }
    }
    
    getline(in, line, '\n'); //Skip the line "#PROFILES"
    
    /** Begin the reading of the profiles **/
    unsigned int counter_chromosomes = 0;
    unsigned int counter_samples = 0;
    unsigned int counter_seg = 0;
    
    while (counter_samples < instance._m && in.good())
    {
        getline(in, line, '\n');
        if(!line.empty() && !(line == "#EDGES"))
        {
            std::stringstream csline(line);
            getline(csline, line, ':'); //Skip the name Li of the leaf
            getline(csline, line, '\n');
            std::stringstream sline(line);
            
            counter_chromosomes = 0;
            while(!sline.eof())
            {
                std::string sample;
                std::getline(sline, sample, '|');
                std::stringstream ssample(sample);
                
                counter_seg = 0;
                while(!ssample.eof())
                {
                    getline(ssample, value, ' ');
                    if(!value.empty())
                    {
                        if(counter_samples >= instance._m
                           || counter_chromosomes >= instance._num_chr
                           || counter_seg >= instance._n[counter_chromosomes])
                        {
                            throw std::runtime_error("ERROR: the input C format is wrong or the number of chromosomes/leaves/corresponding segments is inconsistent with the specified PARAMS");
                        }
                        else
                        {
                            instance._F[counter_chromosomes][counter_samples][counter_seg] = atof(value.c_str());
                        }
                        ++counter_seg;
                    }
                }
                
                if(counter_seg != instance._n[counter_chromosomes])
                {
                    std::stringstream tmp;
                    tmp << "ERROR: inconsistent number of segments for chromosome " << counter_chromosomes;
                    
                    throw std::runtime_error(tmp.str());
                }
                
                ++counter_chromosomes;
            }
            
            if(counter_chromosomes != instance._num_chr)
            {
                throw std::runtime_error("ERROR: inconsistent number of chromosomes");
            }
            
            ++counter_samples;
        }
    }
    
    if(counter_samples != instance._m)
    {
        throw std::runtime_error("ERROR: inconsistent number of samples");
    }

    return in;
}
