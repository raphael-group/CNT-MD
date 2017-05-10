#include "basic_types.h"
#include "tripletarchitect.h"
#include <string>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage " << argv[0] << ": y1 y2" << std::endl;
        return 1;
    }
    
    IntMatrix y1(1);
    IntMatrix y2(1);
    
    std::vector<std::string> s1, s2;
    boost::split(s1, argv[1], boost::is_any_of(",;"));
    boost::split(s2, argv[2], boost::is_any_of(",;"));
    
    if (s1.size() != s2.size())
    {
        std::cerr << "Error: different sizes of y1 and y2" << std::endl;
        return 1;
    }
    
    const int n = s1.size();
    for (int i = 0; i < n; ++i)
    {
        y1[0].push_back(boost::lexical_cast<int>(s1[i]));
        y2[0].push_back(boost::lexical_cast<int>(s2[i]));
    }
    
    TripletArchitect triplet(y1, y2);
    triplet.init();
    triplet.solve(-1, -1);
    
    std::cout << triplet.getTree();
    triplet.getTree().writeDOT(std::cout);
    
    return 0;
}
