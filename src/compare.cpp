#include "basic_types.h"
#include "comparison.h"
#include "fmcsolution.h"
#include <string>
#include <fstream>

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <FILE1> <FILE2> where" << std::endl
                  << "  <FILE1> is copy-number tree filename for the true tree" << std::endl
                  << "  <FILE2> is copy-number tree filename for the inferred tree" << std::endl;
        
        return 1;
    }
    
    std::string filename1(argv[1]);
    std::string filename2(argv[2]);
    
    FMCSolution T1, T2;

    std::ifstream inFile1(filename1.c_str());
    if (!inFile1.good())
    {
        std::cerr << "ERROR: could not open '" << filename1 << "' for reading" << std::endl;
        return 1;
    }
    inFile1 >> T1;
    inFile1.close();
    
    std::ifstream inFile2(filename2.c_str());
    if (!inFile2.good())
    {
        std::cerr << "ERROR: could not open '" << filename2 << "' for reading" << std::endl;
        return 1;
    }
    inFile2 >> T2;
    inFile2.close();
    
    Comparison comp(T1, T2);
    if (!comp.init())
    {
        std::cerr << "Incompatible trees" << std::endl;
        return 1;
    }
    
//    comp.printBpGraph(std::cout);
    
    const int k = T1.getTree().k();
    for (int i = k-1; i < 2*k-1; ++i)
    {
        std::cout << "Leaf " << i << " of T_1 corresponds to leaf "
                  << comp.leaf1ToLeaf2(i) << " of T_2" << std::endl;
    }
    
    std::cout << "Delta(T_1) = " << T1.getTree().cost() << std::endl;
    std::cout << "Delta(T_2) = " << T2.getTree().cost() << std::endl;
    std::cout << "RF = " << comp.robinsonFoulds() << std::endl;
    std::cout << "delta M = " << comp.deltaM() << std::endl;
    std::cout << "leaf consistency = " << comp.leafConsistency() << std::endl;
    
    return 0;
}
