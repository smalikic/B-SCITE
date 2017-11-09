# B-SCITE
========


## Content

Here we present description of B-SCITE software and detais about running ddClone, OncoNEM and SCITE used in the benchmarking.

## Description
--------------


**B-SCITE** is a software package to compute trees of tumor evolution by integrative use of single-cell and bulk sequencing data. Full details about the use of the software will be added in the following days.

## Input data
-------------

Input to B-SCITE consists of single-cell matrix, which is ternary matrix D of dimension nxm, where n denotes number of mutations and m denotes number of single cells obtained in sequencing experiment, and bulk sequencing derived matrix containing details about read counts for each of n mutations. Entries of single cell matrix are 0,1 and 3, coding respectively for absence, presence or missing value for mutation calls. We assume that each mutation 

## Running ddClone 
----------------------------------

ddClone was run according to specifications from (cite ddClone paper). True simulated purity value was given as the input and single-cell data input was obtained after preprocessing simulated single cell matrix D by using Single-Cell Genotyper (cite Single-Cell genotyper paper) that provides genotype matrix that is desired input for ddClone. Tool was run for 300 iterations. 

## Running OncoNEM

As per its specifications, OncoNEM was run after pre-processing single-cell input matrix D by removing single cells with reference genotypes. 

## Running SCITE

The input to SCITE is single-cell matrix used also as the input for B-SCITE. 3 repetitions were run, each with 200 000 iterations. 

## Running B-SCITE

B-SCITE was also run for 3 loops and 200 000 iterations using the input data specified above.
