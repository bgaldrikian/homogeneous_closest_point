# unit_tests.exe

These tests rely on a unit test framework defined in

../../../common/util/unit_test.h

The tests themselves are divided into three source files, organizing them into
categories:

./tests/discrete_tests.cpp
./tests/hcp_tests.cpp
./tests/linalg_tests.cpp
./tests/hcp_linalg_tests.cpp

## Tests

### Discrete Tests

The tests in discrete_tests.cpp help ensure that the iterators defined in

../../common/math/discrete.h

work as expected.

### Linear Algebra Tests

The tests in linalg_tests.cpp help ensure that orthogonal basis creation
using QR decomposition is accurate and stable for a range of dimensionalities.
This is an important core piece of the HCP algorithm.

### HCP Tests

The tests in hcp_tests.cpp are high-level function tests, setting up simple use
cases for HCP such as intersection and time-of-intersection tests, using simple
geometric shapes such as squares and triangles (in two dimensions) and cubes
and tetrahedrons (in three dimensions).

### HCP Linear Algebra Tests

The test in hcp_linalg_tests.cpp exercises a specialized QR decomposition and
orthogonal basis creation function that is used by HCP.  It ensures that the
determinant calculated by the QR decomposition and basis construction matches
that calculated by a more traditional means (LU decompostion).

## Building

Run make in the current folder.  The executable hcp_unit_tests.exe will be
created in the directory ../bin/

## Running

Execute the command

../bin/unit_tests.exe [substring_1] [substring_2] [...]

Any argument is interpreted as a substring filter to be applied to the test
name.  For example,

../bin/unit_tests.exe HCP_

will cause only the tests with names that contain "HCP_" to be run.

The filters have an "or" relationship, so multiple arguments will run the
tests that contain substring_1 or substring_2, etc.

If no arguments are given, all tests are run.
