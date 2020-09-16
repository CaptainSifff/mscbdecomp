# MSCBDECOMP - The minimum split checkerboard decomposition
This library is a Fortran implementation of the minimum split checkerboard decomposition(MSCBD) by
Che-Rung Lee in SIAM Vol 35, No 2, pp. C143-C171.
The MSCBD is useful to obtain sparse representations of the matrix exponential.
Stated more precisely given a sparse matrix M the MSCBD will find an exact decomposition

![equation](https://latex.codecogs.com/gif.latex?M%20%3D%5Csum_i%5EN%20M_i%20%5C%5C%20%5Ctext%7Bsuch%20that%20%7D%20e%5E%7B%5Cdelta%20t%20M%7D%20%3D%20%5Cprod_i%5EN%20e%5E%7B%5Cdelta%20t%20M_i%7D%20&plus;%20%5Cmathcal%7BO%7D%28%5Cdelta%20t%5E2%29)

The number of families N depends on the sparsity of M. The decisive property is that the exponentials can now be exactly evaluated in linear time.
With the help of splitting methods even higher order approximations can be built as outlined in the appendix of [arXiv:2009.04491](https://arxiv.org/abs/2009.04491)

## AUTHOR
Florian Goth, Universität Würzburg, SFB1170, Projekt Z03

## LICENSE
MIT License, see LICENSE file.

## LICENSE

MIT License

## PREREQUISITES/ENVIRONMENT

A standard Desktop System with a fortran 2003 compatible compiler.

## CONFIGURATION

## TESTING

## USAGE

If you want to use this library/set of routines, start in mscbdecomp.f90.
We provide the mat2vert function that takes a Fortran matrix and converts it to
our internal graph structure. The you can apply MvG_decomp which applies the 
Misra - van Gries algorithm to decompose the matrix.
If you want to use our helper Exponential Objects you can create them
via a call to createFullExponentialfromGraphData()

## Notes
This code utilizes the Misra and Gries edge coloring algorithm.
Info can be found here:
https://thorehusfeldt.files.wordpress.com/2010/08/gca.pdf

https://en.wikipedia.org/wiki/Misra_%26_Gries_edge_coloring_algorithm

Here's the original paper on the algorithm:
http://www.cs.utexas.edu/users/misra/psp.dir/vizing.pdf

It features an implementation of the Misra-van-Gries Edge coloring algorithm.
Stuff that can also be found is a path structure that mimics the basic behaviour of 
the C++ vector (due to a lack of templates only for two integers.)
as well as a quicksort on integers and a binary search on integers.
