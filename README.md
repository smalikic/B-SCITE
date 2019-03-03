# Description of B-SCITE and details of comparisons against ddClone, OncoNEM and SCITE
======================================================================================




**B-SCITE** is a software package to compute trees of tumor evolution by integrative use of single-cell and bulk sequencing data. Full details about the use of the software will be added in the following days.

## Software requirements and installation

**B-SCITE** is currently supported on Linux OS. Installation requires 
(i)  C++ compiler (we are currently using GNU Compiler Collection - gcc version 5.2.0)
(ii) CPLEX Optimization Studio Version 12.5.1 or later.
In makefile, located inside folder src, Variable CC encoding path to C++ compiler and variables related to the CPLEX Optimization Studio (CPLEX_DIRECTORY and CPLEX_BUILD) should be set accordingly. Currently, we also provide an example of how these variables are set in our system. Once their values are updated, running single command "make" from within source directory compiles the source-code producing the executable bscite.exe. Details of the input data and running bscite.exe are provided below.

## Input data
-------------

Input to B-SCITE consists of single-cell matrix, which is ternary matrix D of dimension nxm, where n denotes number of mutations and m denotes number of single cells obtained in sequencing experiment, and bulk sequencing derived matrix containing details about read counts for each of n mutations. Entries of single cell matrix are 0,1 and 3, coding respectively for absence, presence or missing value for mutation calls. Folder testExample contains an example input of bulk file (bulkFile-n_50.txt) and SC matrix (SCFile-n_50-m_100.txt). i-th row of SC file and (i+1)-th row of bulk file must correspond to the same mutation. 

## Running B-SCITE

Simple Python script run_B-SCITE.py for running B-SCITE is provided inside folder testExample. Description of the parameters is also provided inside this file. In order to run B-SCITE it suffices to adjust related parameters in run_B-SCITE.py and run command "python run_B-SCITE.py". Please note that the values of number of iterations and repetitions provided in run_B-SCITE.py are set to very small value in order to facilitate verifying the correctness of installation, variable assignments etc. For the real application, these numbers should be significantly higher (we recommend at least 3 repretitions and at least several hundred-thousands of iterations.

## Interpreting Output


Assuming that we keep value of "./exampleOutput" as prefix to the B-SCITE output files and run run_B-SCITE.py script, after completing running of the script, we can observe that B-SCITE reports several files, all with prefix "exampleOutput". Among these files, the best tree can be obtained from the parent vector reported in the line starting with PARENT_VECTOR_OPTIMAL_COMBINED_SCORE) in exampleOutput.matrices file. Assuming zero indexing, the first mutation listed in the bulk file is given index 0, the second mutation is given index 1 etc. In parent vector of a given mutation tree, i-th entry represents the index of parent of mutation with index i.

For example of interpreting parent vector, consider the first few entries in the parent vector reported in the exampleOutput.matrices file (this file is also available in this repository). Based on the bulk file (exampleInput.bulk) index 0 correspond to mut0, index 1 corresponds to mut1, ... , index 49 corresponds to mut49 and 50 corresponds to the root node. Consequently, 1 4 50 10 3 7 ... can be interpreted as follows: 

1 is at position 0 in the parent vector  => parent of mutation with index 0 is mutation with index 1 (i.e. mut1 is parent of mut0)

4 is at position 1 in the parent vector  => parent of mutation with index 1 is mutation with index 4 (i.e. mut4 is parent of mut1)

50 is at position 2 in the parent vector => as we have 50 mutations in the input, 50 corresponds to the root implying that parent of mutation with index 2 is root (i.e. root is parent of mut2)

10 is at position 3 in the parent vector => parent of mutation with index 3 is mutation with index 10 (i.e. mut10 is parent of mut3)

3 is at position 4 in the parent vector  => parent of mutation with index 4 is mutation with index 3 (i.e. mut3 is parent of mut4)

7 is at position 5 in the parent vector  => parent of mutation with index 5 is mutation with index 7 (i.e. mut7 is parent of mut5).
 
 
Alternative way of obtainining the tree reported by B-SCITE is as follows: exampleOutput.matrices file also stores n x n matrix A encoding ordering relations between two mutations in the best-scoring tree T reported by B-SCITE. i-th row of A corresponds to i-th mutation from the input data. Analogous applies to i-th column. In other words, A[i,j]=1 if and only if i = j or mutation i is placed as an ancestor of mutation j in T. Note that A is equivalent to the matrix obtained from ancestry matrix of T by removing row and column corresponding to the root node (please refer to the manuscript for the definition of ancestry matrix).

Note that, in addition to reporting the best, B-SCITE also reports the second best scoring tree.

We are currently working on the code for reporting the output in other formats (e.g. graphviz file) and adding optional number of top scoring trees to be reported (currently .matrices contains the information about two best scoring trees).



### Compression of mutation trees into clonal trees

Clustering is performed along the chains formed by nodes x that lie between nodes A and B such that each of A and B is different and is either root node or has at least two descendants (i.e. is node where branching occurs). Also, we assume that the path between A and B does not contain a node C having more than one descendants. Code for performing clustering of mutations along such chains is available in the folder VAFclusterEMpackage. It can be installed from the provided source code by calling install.packages(PATH_SOURCE_FOLDER, repos=NULL, type="source"), where PATH_SOURCE_FOLDER encodes the path to the content of VAFclusterEM folder (for example, ./VAFclusterEMpackage/VAFclusterEM/). Folder VAFclusterEMpackage contains the examples of input (example.input) and output (example.output) of running mutation clustering script. In the input file, mutations are identified via unique ID provided in the first row and corresponding entry in the second row represent VAFs of the mutation. In the output file, entries in the second row represent unique identifier of the cluster where corresponding mutation is assigned. 

## Details about running tools in benchmarking step
---------------------------------------------------

### ddClone 

ddClone was run according to specifications from (cite ddClone paper). True simulated purity value was given as the input and single-cell data input was obtained after preprocessing simulated single cell matrix D by using Single-Cell Genotyper (cite Single-Cell genotyper paper) that provides genotype matrix that is desired input for ddClone. Tool was run for 300 iterations. 

### OncoNEM

As per its specifications, OncoNEM was run after pre-processing single-cell input matrix D by removing single cells with reference genotypes. 

### SCITE

The input to SCITE is single-cell matrix used also as the input for B-SCITE. 3 repetitions were run, each with 200 000 iterations. 

### B-SCITE

B-SCITE was also run for 3 repetitions, each with 200 000 iterations, using the input data as specified above.


