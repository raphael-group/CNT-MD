#include "fmcsolution.h"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>



FMCSolution::FMCSolution()
    : _tree()
    , _M()
    , _input()
{
}


FMCSolution::FMCSolution(const CopyNumberTree& tree, const DoubleMatrix& M, const InputInstance& input)
    : _tree(tree)
    , _M(M)
    , _input(input)
{
}

void FMCSolution::writeDOT(std::ostream &out) const
{
    out << "digraph FMC {" << std::endl;
    std::vector<std::string> colors;
    colors.push_back("blue");
    colors.push_back("red");
    colors.push_back("green");
    colors.push_back("purple");
    colors.push_back("yellow");
    colors.push_back("brown");
    colors.push_back("pink");
    colors.push_back("grey");
    colors.push_back("cyan");

    std::stringstream tstream;
    _tree.writeDOT(tstream);
    std::string line;

    getline(tstream, line, '\n');
    out << "\t" << "subgraph T {" << std::endl;

    for(unsigned int i = 0; i < 2 * _tree.k()-1; ++i) {
        getline(tstream, line, '\n');
        out << "\t" << line << std::endl;
    }

    out << "\t" << "}" << std::endl;

    out << "\t" << "subgraph F {" << std::endl;
    for(unsigned int p = 0; p < _input.m(); ++p)
    {
        out << "\t\tS" << p << " [label=\"";
        for(unsigned int chr = 0; chr < _input.numChr(); ++chr)
        {
            for(unsigned int s = 0; s < _input.n()[chr]; ++s)
            {
                out << std::setprecision(3) <<_input.F()[chr][p][s];
                if(s < _input.n()[chr] - 1 || chr < _input.numChr() - 1)
                {
                    out << " ";
                }
            }
            if (chr < _input.numChr() - 1)
            {
                out << "| ";
            }
            else
            {
                out << "\",color=" << colors[p % colors.size()] << ",shape=box]" << std::endl;
            }
        }
    }
    out << "\t" << "}" << std::endl;

    do {
        getline(tstream, line, '\n');
        if(line != "}" && !line.empty())
            out << "" << line << std::endl;
    } while(!tstream.eof());

    for(unsigned int i = 0; i < _tree.k(); ++i)
    {
        int remapped_i = i + _tree.k() - 1;
        for(unsigned int p = 0; p < _input.m(); ++p)
        {
            out << "\t" << remapped_i << " -> S" << p << " [label=\""
                << std::setprecision(2) << _M[p][i] << "\",style=dotted,color= "
                << colors[p % colors.size()]<<"]" << std::endl;
        }
    }
    out << "}" << std::endl;
}

std::ostream& operator<<(std::ostream& out, const FMCSolution& instance)
{
    out << instance._tree;

    out << "#MIX-PARAMS" << std::endl;
    out << instance._input.numChr() << " #number of chromosomes" << std::endl;
    out << instance._input.m() << " #number of samples" << std::endl;
    for(unsigned int chr = 0; chr < instance._input.numChr(); ++chr)
    {
        out << instance._input.n()[chr] << " ";
    }
    out << "#number of segments for each chromosome" << std::endl;

    out << "#PROPORTIONS" << std::endl;
    for(unsigned int p = 0; p < instance._input.m(); ++p)
    {
        for(unsigned int i = 0; i < instance._tree.k(); ++i)
        {
            if(i < instance._tree.k() - 1)
            {
                out << instance._M[p][i] << " ";
            }
            else
            {
                out << instance._M[p][i] << std::endl;
            }
        }
    }

    out << "#SAMPLES" << std::endl;
    for(unsigned int p = 0; p < instance._input.m(); ++p)
    {
        out << p << " :";
        for(unsigned int chr = 0; chr < instance._input.numChr(); ++chr)
        {
            for(unsigned int s = 0; s < instance._input.n()[chr]; ++s)
            {
                out << " " << instance._input.F()[chr][p][s];
            }

            if(chr < instance._input.numChr() - 1)
            {
                out << " |";
            }
        }
        out << std::endl;
    }
    return out;
}

std::istream& operator>>(std::istream& in, FMCSolution& instance)
{
    std::string line;
    std::string value;
    unsigned int m = 0;
    std::stringstream tstream;

    in >> instance._tree;

    //store #PARAMS
    tstream << "#PARAMS" << std::endl;

    //store #number of chromosomes
    getline(in, line, '\n');
    tstream << line << std::endl;

    //store #number of samples
    getline(in, line, '\n');
    tstream << line << std::endl;
    std::stringstream tline(line);
    tline >> m;

    //store #number of segments
    getline(in, line, '\n');
    tstream << line << std::endl;

    //skip #PROPORTIONS
    getline(in, line, '\n');

    instance._M.clear();
    instance._M = DoubleMatrix(m, DoubleArray(instance._tree.k(), 0.0));

    unsigned int counter_samples = 0;

    while (counter_samples < m && in.good())
    {
        getline(in, line, '\n');
        
        if(!line.empty() && !(line.at(0) == '#'))
        {
            std::stringstream sline(line);
            unsigned int counter_leaves = 0;

            while(!sline.eof())
            {
                getline(sline, value, ' ');
                if(!value.empty())
                {
                    instance._M[counter_samples][counter_leaves] = atof(value.c_str());
                    ++counter_leaves;
                }
            }

            if (counter_leaves != instance.getTree().k())
            {
                throw std::runtime_error("ERROR: inconsistent number of leaves in #PROPORTIONS");
            }
            ++counter_samples;
        }
    }
    if (counter_samples != m)
    {
        throw std::runtime_error("ERROR: inconsistent number of samples in #PROPORTIONS");
    }

    do {
        getline(in, line, '\n');
        tstream << line << std::endl;
    } while(!in.eof());
    
    tstream >> instance._input;
    return in;
}

