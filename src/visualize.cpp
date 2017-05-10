#include "basic_types.h"
#include "fmcsolution.h"
#include "copynumbertree.h"
#include <string>
#include <fstream>

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <FILE> where" << std::endl
                  << "  <FILE> is solution filename or - for STDIN" << std::endl;
        
        return 1;
    }
    
    std::string filename(argv[1]);

    FMCSolution sol;
    if (filename == "-")
    {
        std::cin >> sol;
    }
    else
    {
        std::ifstream inFile(filename.c_str());
        if (!inFile.good())
        {
            std::cerr << "ERROR: could not open '" << filename << "' for reading" << std::endl;
            return 1;
        }
        inFile >> sol;
        inFile.close();
    }
    
    sol.writeDOT(std::cout);
    
    return 0;
}
