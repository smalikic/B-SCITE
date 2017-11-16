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

Simple Python script run_B-SCITE.py for running B-SCITE is provided inside folder testExample. Description of the parameters is also provided inside this file. In order to run B-SCITE it suffices to adjust related parameters in run_B-SCITE.py and run command "python run_B-SCITE.py".

## Interpreting Output

Assuming that in the previous step, when running run_B-SCITE.py, we set prefix of B-SCITE output filenames to "./example", B-SCITE reports three output files (in this case stored inside folder testExample) with the following names and meaning:
#### (i) example.matrices file 
  This file stores n x n matrix A encoding ordering relations between two mutations in the best-scoring tree T reported by B-SCITE. i-th row of A corresponds to i-th mutation from the input data. Analogous applies to i-th column. In other words, A[i,j]=1 if and only if i = j or mutation i is placed as an ancestor of mutation j in T. Note that A is equivalent to the matrix obtained from ancestry matrix of T by removing row and column corresponding to the root node.
#### (ii) example.gv
  Best-scoring tree represented in Graphviz format. i-th mutation from the input data is labelled as i and root node is labelled as (n+1).
#### (iii) example.newick
  Best-scoring tree represented in Newick format. 

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


