# Homogeneous Closest Point Algorithm

2018-08-08, B. Galdrikian

## Introduction

This directory contains a reference implementation of the Homogeneous Closest
Point (HCP) algorithm (see [1]), as well as test programs.  These programs
include the benchmarks used to measure performance for the results presented
in [1].

The files contain two variants of the algorithm: a general version that works
for any number of spatial dimensions, and dimension-specific versions
specialized for 1, 2, 3, and 4 dimensions.  Depending on the compiler and
processor, the specialized versions can be up to several times faster than the
general one.

All of the code is written in C++ (uses C++ 11 features). The files needed
for either version of the algorithm are given below.


## Reference Implementations

### General Algorithm

#### Files Required

../common/types.h  
../common/math/discrete.h  
../common/math/linalg.h  
./src/hcp_defs.h  
./src/general/hcp.cpp  
./src/general/hcp.h  
./src/general/hcp_linalg.h

#### To Use

Include the file ./src/general/hcp.h, and compile ./src/general/hcp.cpp with
your project (or link it with your project in some way).

HCP is implemented in the functions hcp_solve_i and hcp_solve, defined in
hcp.h.  These function require an implementation of HCP_Halfspace_Set, an
abstract base class defined in ./src/hcp_defs.h.  Some common implementations
are given in

./src/general/sets/

These include a polytope, an ellipsoid, and a class which represents the union
of two sets, in hcp_shape_pair.h (this is useful for testing intersection of
two convex shapes).

#### Time of Intersection

For ToI, in addition to the files listed above, you need these two files:

./src/general/hcp_toi.cpp  
./src/general/hcp_toi.h

### Specializations for D = 1, 2, 3, and 4

#### Files Required

../common/types.h
./src/hcp_defs.h
./src/specialized/hcpd.h
./src/specialized/hcpd_math.h

In addition, you must have the files required for the dimensionality of the
problem to be solved.  For example, for solving D = 3 problems, you will
need:

./src/specialized/hcp3d.h

#### To Use

Include the dimension-specific files (./src/general/hcp3d.h for example).

Dimension-specialized HCP is implemented in the functions hcpd_solve_i<D> and
hcpd_solve<D>, defined in hcpd.h.  These function require an implementation of
HCP_Halfspace_Set, just as the general algorithm does.  Implementations which
are templated in D are given in

./src/specialized/sets/

The sets templated currently included are HCPD_Static_Polytope<D> in
hcpd_static_polytope.h, and HCPD_Shape_Pair<D> in hcpd_shape_pair.h, a class
which represents the union of two sets.

#### Time of Intersection

For ToI, in addition to the files listed above, you need this files:

./src/specialized/hcpd_toi.h

*NOTE* The ToI calculation in the reference implementations solves the
problem in one higher dimension, adding a temporal dimension internally.
This means that if you are solving a 3-dimensional ToI problem, you will
need to include

./src/specialized/hcp4d.h

in order for the code to compile.

## Examples

The directory

./examples/

contains minimal programs which demonstrate what is required to compile
and perform several types of calculations using the HCP algorithm.

See the readme:

./examples/readme.txt

for more.

## Test Programs

The directory

./tests/

contains test programs including a performance test, soak test, and unit
tests.  See the readme:

./tests/

for more.

## Licensing

The software in this directory is open-source, licensed under the MIT license
below:

```
// Copyright (c) 2014-2018 NVIDIA Corporation
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
```

## Citations

[1] "A Unified Closest Point and Linear Optimization Algorithm,"
B. Galdrikian.  (To be published.)
