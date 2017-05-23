# CNT-MD: Coordinate-descent algorithm for the Copy-Number Tree Mixture Deconvolution (CNTMD) Problem
============================================

[License](LICENSE.txt)

Reference publication:

[Simone Zaccaria, Mohammed El-Kebir, Gunnar W. Klau, Benjamin J. Raphael: The Copy-Number Tree Mixture Deconvolution Problem and Applications to Multi-sample Bulk Sequencing Tumor Data. RECOMB 2017: 318-335](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_20)


For any questions regarding this tool or its usage, please contact:
  - Simone Zaccaria:  zaccaria  [at]  princeton  [dot]  edu
  - Mohammed El-Kebir:  melkebir  [at]  princeton  [dot]  edu

## Contents

  1. [Dependencies](#dep)
  2. [Compilation instructions](#comp)
  3. [Usage instructions](#usage)
  4. [Available data](#data)
	5. [Implementation details](#imp)

## <a name="dep"></a>Dependencies

* [LEMON Graph library](http://lemon.cs.elte.hu) (>= 1.3)
* [ILOG CPLEX](http://www.ibm.com/developerworks/downloads/ws/ilogcplex/) (>= 12.0)
* Boost libraries

[Graphviz](http://www.graphviz.org) is required to visualize the resulting DOT files, but is not required for compilation.

## <a name="comp"></a>Compilation instructions

To compile `CNT-MD`, execute the following commands from the root of the repository:

    mkdir build
    cd build
    cmake ..
    make
    
Note: On Mac OS >= 10.9, you have to ensure that LEMON is linked against the old C++ standard library by editing LEMON's `CMakeLists.txt` file as follows.

	if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
	  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
	endif()
    
In case CMake fails to detect either CPLEX or LEMON, run the following command with adjusted paths:

	cmake \
	-DLIBLEMON_ROOT=~/lemon \
	-DCPLEX_INC_DIR=~/ILOG/cplex/include/ \
	-DCPLEX_LIB_DIR=~/ILOG/cplex/lib/x86-64_osx/static_pic \
	-DCONCERT_LIB_DIR=~/ILOG/concert/lib/x86-64_osx/static_pic \
	-DCONCERT_INC_DIR=~/ILOG/concert/include/ ..
	
or use the tool `ccmake` to set the same corresponding paths.

The compilation results in the following executables in the `build` directory:

* `mixcnp`
* `visualize`
* `compare`

## <a name="usage"></a>Usage instructions

### Input format

Example:

	#PARAMS
	2 #number of chromosomes
	4 #number of samples
	6 6 #number of segments for each chromosome
	#SAMPLES
	1 : 2.14 1.43 1.43 3.71 2.0 2.0 | 3.29 3.29 4.0 2.0 2.0 1.29
	2 : 1.0 2.0 2.0 1.0 2.0 2.0 | 1.0 0.0 0.0 1.0 2.0 2.0
	3 : 2.9 2.0 2.0 3.85 2.0 2.0 | 2.9 2.85 3.8 1.95 2.0 1.05
	4 : 1.55 1.99 1.99 1.83 2.0 2.0 | 1.56 0.84 1.11 1.28 2.0 1.72

The first line must be `#PARAMS`. The three subsequent lines indicate the number of chromosomes, the number of samples and the number of segments for each chromosome (separated by a space character). The fourth line must be `#SAMPLES`. Next, for each sample the fractional copy numbers are given, such that the fractions of segments on the same chromosome are separated by spaces, whereas the chromosomes are separated by '|'.
*IMPORTANT:* Each sample must contain the same set of segments per chromosome. Therefore if the fractional copy numbers have been called independently for each sample, a procedure is needed to ensure this.

### Output format

Example:

    #PARAMS
    2 #number of chromosomes
    3 #number of leaves
    6 6 #number of segments for each chromosome
    #PROFILES
    0 : 2 2 2 2 2 2 | 2 2 2 2 2 2
    1 : 1 1 1 3 2 2 | 3 3 3 2 2 2
    2 : 1 2 2 1 2 2 | 1 0 0 1 2 2
    3 : 0 0 0 3 2 2 | 4 4 4 2 2 2
    4 : 3 2 2 4 2 2 | 3 3 4 2 2 1
    #EDGES
    0 -> 2
    0 -> 1
    1 -> 3
    1 -> 4
    #EVENTS
    0 0 1 0 2 -1
    0 0 1 3 3 1
    0 0 2 0 3 -1
    0 0 2 1 2 1
    0 1 3 0 4 -1
    0 1 3 2 4 1
    0 1 4 0 0 1
    0 1 4 0 3 1
    1 0 1 0 1 1
    1 0 1 2 2 1
    1 0 2 0 2 -1
    1 0 2 1 3 -1
    1 1 3 0 2 1
    1 1 4 2 4 1
    1 1 4 3 5 -1
    #MIX-PARAMS
    2 #number of chromosomes
    4 #number of samples
    6 6 #number of segments for each chromosome
    #PROPORTIONS
    0.0 0.29 0.71
    1.0 0.0 0.0
    0.05 0.0 0.95
    0.72 0.0 0.28
    #SAMPLES
    1 : 2.14 1.43 1.43 3.71 2.0 2.0 | 3.29 3.29 4.0 2.0 2.0 1.29
    2 : 1.0 2.0 2.0 1.0 2.0 2.0 | 1.0 0.0 0.0 1.0 2.0 2.0
    3 : 2.9 2.0 2.0 3.85 2.0 2.0 | 2.9 2.85 3.8 1.95 2.0 1.05
    4 : 1.55 1.99 1.99 1.83 2.0 2.0 | 1.56 0.84 1.11 1.28 2.0 1.72

The first line is `#PARAMS`. 
The three subsequent lines indicate the number of chromosomes, the number *k* of leaves that have been inferred and the number of segments for each chromosome (separated by a space). 
The fourth line is `#PROFILES`. 
Next, for each node corresponding to a clone in the copy number tree, the copy number profile is given, such that the chromosomes are separated by the symbol `|` and the copy numbers on each chromosome are separated by spaces. 
Each clone is uniquely identified by a number at the beginning of the line followed by the copy number profile after `:`. 
The *k* leaves, corresponding to the *k* extant clones that we inferred, correspond to the last *k* clones. 
The next line is `#EDGES`.
Each subsequent line corresponds to an edge of the tree that is given in the form `i -> j` such that the edge goes from the clone with id `i` to the clone with id `j`. 
The next line is `#EVENTS`. 
Each subsequent line corresponds to an interval event that is given in the form of `c i j s t`, where `c` is the chromosome that is affected by the interval event, `(i, j)` identifies the edge labeled by the corresponding event, and `(s, t)` define the start and end positions of the interval event (both bounds are inclusive, and are 0-based). The remaining lines replicate the provided input instance, described in the previous section.

### Example execution

To run `mixcnp`:

	./mixcnp input.samples -Z 120 -d -j 2 -k 3 -m 4000 -ni 5 -ns 50 -nt 2 -o result.out -s 120 -t 0.001 -ss 12 &> log

`input.samples` is the file containing the input samples in the input format. 
The rightmost bound *Z* is set to 120 using `-Z` and `-lbZ` is left to the default value of 0, therefore the best value of Lambda_max is searched for in the interval *[0, 120]*. 
The option `-d` is used to force the presence of the normal diploid clone, and this is indeed needed to study the normal admixture. 
Given parameters `-j 2 -ni 5 -ns 50 -s 600`, we can estimate an upper bound on the total running time.
There are 2 workers (`-j 2`) running in parallel and each of these will run roughly half of the total number of 50 starting points (`-ns 50`). 
Our method runs an iterative procedure composed of at most 5 steps (`-ni 5`) for each starting point.
Each step of the iterative procedure is composed of two steps: a C-step and a U-step. such that the parameter `-s 600` specifies the time limit for each step of each iteration.
The U-step typically does not affect the total running time.
While, the parameter `-s 120` specifies the time limit 120 in terms of seconds for each C-step.
Hence, the running time from each starting point is at most of `*120*5 sec = 600 sec = 10 min`.
As such, the execution of all the starting points will last for up to `25*10 min = 250 min = 4,2 hr` since each worker runs 25 starting points.
Moreover, the size of the interval for Lambda_max is 100 (`-Z 100` and `-lbZ` left to default 0) and the number of iteration of the searching scheme is logarithmic, i.e. approximately 8 iterations, since we are using the default binary search.
Therefore the main algorithm will be repeated by the seatchin scheme for 8 times and the entire procedure will be executed for at most `8*4,2 hr = 34 hr` (using a higher number of workers equal to 10 we can already improve to the expected total running time to 7 hr). Since the number of threads of each worker is set to 2, a total of 4 cores will be used on the same machines and they must be available to guarantee the expected running time. Furthermore, each worker can use at most 4 GB (`-m`) and this guarantees that the memory of the machine will not be exhausted by using a total of `2*4 GB = 8 GB`.
The threshold value `-t 0.001` means that a difference between each observed fractional copy number and the result is tolerated up to 0.001 (below that threshold, drops in the objective are not considered as an improvement). The parameter `-ss 12` determines the random seed that will be used and this is important for replicating the execution. Lastly, the final result will be written in the output file `result.out` and the log is written by redirection in the file `log`.

To visualize the resulting copy-number tree:

    ./visualize result.txt > T.dot
    dot -Tpdf T.dot -o T.pdf

### Usage suggestions

1. The real actual running time strongly depends on the complexity of the input. To choose the best value for the time limit `-s` that allows to minimize the needed running time and obtain the best performance, we suggest the user to test the input instance with an increasing time limit `-s` (using for example an interval of length 1, i.e. `-Z 100` and `-lbZ 100`) and to analyze how the objective function varies. More careful analysis can be done by considering the gap of the ILP formulation using and higher level of verbosity `-v 4` that measures how far the provided solution is from the optimum. The performance of the method highly benefits from an execution on multiple machines. We can achieve this solution by running different starting points on different machines with the same parameters and different random seeds `-ss` or the interval *[L, R]* where to search Lmabda_max can be split into different machines.

2. The same procedure described above can be symmetrically applied to estimate a number of starting points `-ns` and number `-ni` of iterations that are large enough for estimating the value of the objective function.

3. For large instances, we suggest the user to follow a 2-step procedure: In the first step, multiple starting points with different parameters and a running time relatively low can be used to estimate a smaller interval *[L, R]* containing the best value of Lambda_max with a high probability. In the second step, a higher running time and a higher number of starting points can be used to better explore the objective function on the resulting smaller interval *[L, R]*. 


### Main tools

In the following section we describe the details of the usage and of the parameters for each of the main tools, each corresponding to an executable.

EXECUTABLE | DESCRIPTION
-----------|-------------
`mixcnp`   | Main algorithm
`visualize`| Visualize the resulting output as a dot file
`compare`  | Compare the result with a truth instance

#### mixcnp

This is the main tool that implements the coordinate-descent algorithm for solving the CNTMD problem. The input and output formats are described in the next sections here below. For more detailes about the algorithm, please refer to the reference pubblications. Instead, for more details about the implementations of the algorithm please refer to the Section 4 of this document. The algorithm takes in input a collections of fractional copy numbers obtained from multiple samples, and infer a set of *n* extant clones, the copy-number tree describing their evolution, and the corresponding proportions of the extant clones in the various samples.

     Usage:
        ./mixcnp [--help|-h|-help] -Z int [-d] [-dr] [-e int] [-f] [-j int] -k int
                 [-lbZ int] [-m int] [-ni int] [-ns int] [-nt int] [-o str] [-r int]
                 [-s int] [-ss int] [-t num] [-v int] input
     Where:
            input
              Input file
            --help|-h|-help
              Print a short help message
            -Z int
               Maximum cost of tree considering all the chromosomes. This corresponds to the rightmost bound R of the interval [L, R] where the value of the maximum cost Lambda_max is searched. We suggest to use a large value of R and this can be estimate depending on the total number of genomic segments in the input and the maximum copy number that is allowed in the profiles of the inferred clones.
            -d
              Force one clone to be the normal diploid (default: false)
            -dr
              Deactivate refinement. The refinement step is a final step that aims to minimize the number of events in the final resulting copy-number tree.
            -e int
               Maximum copy number (default: -1, inferred from leaves)
            -f
              Do not fix root to all 2s
            -j int
               Number of workers (default: 2). Each worker solves a starting point independently and in parallel to the other workers. Each worker uses a number of threads specified by the following parameter '-nt' 
            -k int
               Number of leaves, corresponding to the extant clones that we are inferring
            -lbZ int
                 Lower bound for maximum size of tree for all chromosomes (default: 0). This corresponds to the lower bound L of the interval [L, R] where we are searching for the best value of Lambda_max.
            -m int
               Memory limit in MB for each worker (default: -1, disabled)
            -ni int
                Number of iterations per seed (default: 7). This is the maximum number of iterations that are applied to each starting point even if the convergence has not been reached.
            -ns int
                Number of starting seeds (default: 10). This corresponds to the total number of starting points.
            -nt int
                Number of ILP threads (default: 1) that are used by each worker.
            -o str
               Output filename
            -r int
               Mode for searching parsimonious number of events: (1) Binary Search (2) Reverse Iterative (3) Full Iterative (default: 1)
            -s int
               Time limit in seconds for each C-step (default: -1, disabled)
            -ss int
                Random number seed (default: 0) for the generations of the starting points.
            -t num
               Epsilon, threshold level of tolerance for normalized distance (default: 0.0)
            -v int
               Verbosity level from 0 to 4 (default: 1)


#### visualize

This tools is used to visualize the output of the main algorithm. Given in input a file in the output format, it produces a dot file that can be used to draw the inferred results in a pdf using the `dot` program.

     Usage: ./visualize <FILE> where
            <FILE> is solution filename or - for STDIN

#### compare

This tool is used to compare to outputs following the metric used in the reference pubblication. Usually, this is used to compare an inferred result with the ground truth.

     Usage: ./compare <FILE1> <FILE2> where
            <FILE1> is copy-number tree filename for the true tree
            <FILE2> is copy-number tree filename for the inferred tree

## <a name="data"></a>Available data

### Simulated data

In this table we describe the simulated datasets that are made available in this own repository together with the source of the method.

Dataset           | Location of Input          | True Solutions            | DESCRIPTION
------------------|----------------------------|---------------------------|-------------
`1chr`            | data/sim/1chr-input        | data/sim/1chr-input       | Simulated dataset containing instances composed of a single chromosome, 20 segments per chromosome, varying true number of unknown leaves in {4,6,8} (identified by the corresponding labels in the file name {k4,k6,k8}), and varying number of samples in {2,5,10} (identified by the corresponding labels in the filename {m2,m5,m10})
`4chr`            | data/sim/4chr-input        | data/sim/4chr-input       | Simulated dataset containing instances composed of 4 chromosomes, 20 segments per chromosome, varying true number of unknown leaves in {4,6,8} (identified by the corresponding labels in the file name {k4,k6,k8}), and varying number of samples in {2,5,10} (identified by the corresponding labels in the filename {m2,m5,m10})
`full-genome`     | data/sim/full-genome-input | data/sim/full-genome-true | Simulated full-genome dataset containing instances composed of 22 chromosomes, varying number of segments in the chromosomes that reflect the length distribution of the prostate-cancer dataset presented in (Gundem et al., Nature, 2015), true number of leaves equal to 4, and varying number of samples in {2,5,10} (identified by the corresponding labels in the filename {m2,m5,m10}), and randomly generated with different random seeds identified by the following labels in the file name {ss11,ss12,ss13}.

### Real data

Moreover, in `data/real/gundem2015` we make available the input instances that we obtain from 4 patients `{A10, A22, A29, A31}` of the prostate-cancer dataset published in (Gundem et al., Nature, 2015). We obtain the instances by processing the fractional copy numbers in these datasets (that the authors obtained by using the Battenberg tool) and we project the segments on all the samples, such that each sample is defined over the same set of segments. Moreover, we clean the instances by removing outlying segments that are too short or have a unlikely-high fractional copy number. The corresponding instances are available in the following files {`A10-cleaned.samples`, `A22-cleaned.samples`, `A29-cleaned.samples`, `A31-cleaned.samples`} and these other files {`A10-resulting-starts-ends.segments`, `A22-resulting-starts-ends.segments`, `A29-resulting-starts-ends.segments`, `A31-resulting-starts-ends.segments`} contain the starting position (first line) and the ending positions (second line) of each resulting segment.


## <a name="imp"></a>Details on implementation

Our method implements the coordinate-descent algorithm described in the reference pubblication by having in input a number of clones *n* and an interval *[L, R]* containing the value of the desired maximum cost Lambda_max. A searching scheme is applied on the interval *[L, R]* for finding the desired value of the maximum cost Lambda_max. Therefore, the coordinate-descent algorithm is repeated considering different values of Lambda_max in the searching scheme.  In the next subsection "Hostarts" we describe the method used to guarantee the objective function to be monotonically non-increasing and in the subsection "Workers" we describe the parallel implementation of the coordinate descent algorithm. In the last subsection "Searching scheme" we describe the different searching scheme for Lambda_max that are implemented.

### Hotstarts

For each starting point, the coordinate-descent algorithm yields a sequence of pairs of matrices where for consecutive pairs *(C_{q, t}, U_{q, t}), (C_{q, t+1}, U_{q, t+1})* it holds that the objective function is monotonically non-increasing, i.e. *||F - C_{q, t}, U_{q, t}|| >= ||F - C_{q, t+1}, U_{q, t+1}||*. The C-step computing each matrix *C_{t+1}* is implemented as an ILP formulation, whereas the U-step computing *U_{t+1}* is implemented as an LP formulation, and they are both solved using CPLEX. Therefore, to guarantee that the function is monotonically non-increasing, we use the HotStart technology of CPLEX such that each C-step is initialized with the previously computed *C_{t}*, whereas each U-step is initialized with the previously computed *U_{t}*. This guarantees the fact that the objective function is monotonically non-increasing for each step of each starting point.

### Workers

For a given value of maximum cost Lambda_max, our method applies the coordinate-descent algorithm described in the reference pubblication by using different starting points. We can consider and execute each starting point independently. 
As such, we exploit the parallel structure of our method by executing multiple starting points in parallel, when possible. A worker is the basic computational unit that executes the algorithm starting from a single starting point. Therefore, the the total number of workers determine the number of starting points that are executed on parallel. When a worker concludes a job, it executes the algorithm considering the next available (if there is) starting point. The execution of the method ends when the workers executed all the starting points.

### Searching scheme for the best value of Lambda_max

Here, we describe the alternative search schemes that can be used by our method to find the best value of the maximum cost Lambda_max of the corresponding copy-number tree in a starting interval *[L, R]*. By default L is set to *0* and *R* must be chosen by the user. We suggest to choose a large value for *R* that can be estimated from the total number of segments and the maximum copy number in the input. Our method uses by default a search scheme based on binary search for finding the best value of maximum cost Lambda_max within an interval *[L, R]*. This approach is very efficient when considering very large intervals. However, when considering small intervals a linear search could be preferred. As such, we implemented different search scheme to find the best value of Lambda_max in the interval *[L, R]* and the user can choose the method depending the psecific usage. The advantages on each method are described in the following:

1. Binary search. The running time of this search scheme is logarithmic in the size of the interval and it is ideal for considering large intervals. This is the standard search scheme used by our method since the searching interval is usually large. At each iteration, the computation with this search scheme benefits from the fact that the best solution found with Lambda_max equal to L can be used as a starting solution.

2. Linear search from L to R . The running time of this search scheme is linear in the size of the interval and it is ideal when considering small intervals. In fact, the best solution found with a specific value Lambda_t of the maximum cost of Lambda_max can be used as a starting solution for computing the solution with Lambda_{t+1}.

3. Reverse Linear search from R to L . The running time of this search scheme is linear in the size of the interval. This method can be very fast when the best value of Lambda_max is close to the rightmost bound R. This method is therefore 
