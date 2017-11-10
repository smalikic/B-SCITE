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

Input to B-SCITE consists of single-cell matrix, which is ternary matrix D of dimension nxm, where n denotes number of mutations and m denotes number of single cells obtained in sequencing experiment, and bulk sequencing derived matrix containing details about read counts for each of n mutations. Entries of single cell matrix are 0,1 and 3, coding respectively for absence, presence or missing value for mutation calls. We assume that each mutation 

## Running B-SCITE

Coming soon.

## Interpreting Output

Coming soon.

#### Compression of mutation trees into clonal trees

Coming soon.


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


