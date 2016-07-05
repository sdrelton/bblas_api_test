# bblas_api_test

Recently it has become apparent that Batched BLAS (BBLAS) routines,
for computing multiple BLAS operations simultaneously,
are required to obtain maximal performance in a number of HPC applications.

Although there is some initial support for BBLAS in CuBLAS, MKL, and MAGMA,
there is no standard set of functions and calling sequence between these libraries.
Setting such a standard would allow these libraries to interoperate easily and
give users a standard set of functionality for use in their applications.

There is a draft standard available at http://eprints.ma.man.ac.uk/2464/
which is open to comments and suggestions from the community.

At the recent "Workshop on Batched, Reproducible, and Reduced Precision BLAS"
(http://www.netlib.org/utk/people/JackDongarra/WEB-PAGES/Batched-BLAS-2016/)
a number of alternative APIs were suggested by the participants.

## Inside this repository
The aim of this repository is to create sample code for the different APIs
analyzed in the technical report *(FILL IN URL HERE)*.
Please use the repository in combination with the tech. report in order
to fully understand the different approaches.

This repository contains implementations of three possible APIs for BBLAS routines.
The APIs we implement are:
- Separate fixed and variable batch API
- Flag-based API
- Group API

We also include a version of DGEMM using a strided memory layout, as opposed to the
pointer-to-pointer approach used in the other functions.

## Download and compilation
To clone the repository please use the following command.
`git clone git@github.com:sdrelton/bblas_api_test.git`

After this you will need to edit the make.inc file.
