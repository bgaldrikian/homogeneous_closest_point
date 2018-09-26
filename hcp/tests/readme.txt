# Tests

This folder contains subfolders, each containing a test program.  These
programs can all be built at once using the makefile in this directory, or
individually using the makefiles found in the separate program directories.
The programs are:

hcp_perf (in ./perf/)
hcp_soak (in ./soak/)
hcp_unit_tests (in ./unit_tests/)

Running make puts the executables into a ./bin/ folder (this folder will
be created if it doesn't already exist.)

A brief description of the tests is given below.  See the readme.txt in the
program directories for more detailed descriptions.

## hcp_perf

These tests time the execution of HCP with pre-calculated random arrangements
of polytopes, calculating mean test time and variance.  Timing measurements
are also performed on other algorithms for comparison with HCP.

## hcp_soak

This is a test which puts two polytopes in near proximity with one another,
using random arrangements of half-spaces.

This test ensures that HCP terminates within its maximum number of iterations.
It also compares the result with other algorithms and reports discrepancies.

This test has been used to ensure robustness of HCP.  Typically it is run for
one billion trials.

## hcp_unit_tests

This suite of unit tests check low- and high-level functions, used routinely
during code development.
