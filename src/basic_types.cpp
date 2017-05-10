#include "basic_types.h"

lemon::Tolerance<double> g_tol(1e-4);

boost::mt19937 g_rng;

VerbosityLevel g_verbosity;

boost::mutex g_mutex;

boost::mutex g_output_mutex;

int countElements(const IntMatrix arg)
{
    int res = 0;
    for(IntMatrix::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += (*it).size();
    }
    return res;
}

int countElements(const Int3Array arg)
{
    int res = 0;
    for(Int3Array::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += countElements(*it);
    }
    return res;
}

int countElements(const Int4Array arg)
{
    int res = 0;
    for(Int4Array::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += countElements(*it);
    }
    return res;
}

int countElements(const DoubleMatrix arg)
{
    int res = 0;
    for(DoubleMatrix::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += (*it).size();
    }
    return res;
}

int countElements(const Double3Array arg)
{
    int res = 0;
    for(Double3Array::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += countElements(*it);
    }
    return res;
}

int countElements(const Double4Array arg)
{
    int res = 0;
    for(Double4Array::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += countElements(*it);
    }
    return res;
}

int sum_of_elements(const IntArray arg)
{
    int sum = 0;
    for(IntArray::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        sum += *it;
    }
    return sum;
}

std::string timestamp()
{
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    return "[" + boost::posix_time::to_simple_string(now) + "]";
    /**
    time_facet *facet = new time_facet("%d-%b-%Y %H:%M:%S");
    out.imbue(locale(out.getloc(), facet));
    out << "[" << second_clock::local_time() << "]";
    **/
}
