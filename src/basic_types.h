#ifndef _BASIC_H_
#define _BASIC_H_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <lemon/tolerance.h>
#include <random>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <limits>

extern lemon::Tolerance<double> g_tol;

/// Random number generator
extern std::mt19937 g_rng;

extern boost::mutex g_mutex;

extern boost::mutex g_output_mutex;

typedef enum
{
    VERBOSE_NONE,
    VERBOSE_ESSENTIAL,
    VERBOSE_NON_ESSENTIAL,
    VERBOSE_DEBUG
} VerbosityLevel;

extern VerbosityLevel g_verbosity;

typedef std::vector<int> IntArray;
typedef std::vector<IntArray> IntMatrix;
typedef std::vector<IntMatrix> Int3Array;
typedef std::vector<Int3Array> Int4Array;

typedef std::vector<double> DoubleArray;
typedef std::vector<DoubleArray> DoubleMatrix;
typedef std::vector<DoubleMatrix> Double3Array;
typedef std::vector<Double3Array> Double4Array;


typedef std::set<int> IntSet;
typedef std::pair<IntSet, IntSet> IntSetPair;

typedef std::vector<int> HotStart;

int countElements(const IntMatrix arg);
int countElements(const Int3Array arg);
int countElements(const Int4Array arg);

int countElements(const DoubleMatrix arg);
int countElements(const Double3Array arg);
int countElements(const Double4Array arg);

int sum_of_elements(const IntArray arg);

std::string timestamp();

#endif // _BASIC_H_
