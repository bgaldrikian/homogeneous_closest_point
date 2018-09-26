# perf.exe

This program executes various timing tests for the purpose of measuring
scaling behavior, comparing with other algorithms, and regression testing.

This program has two main modes: N-sweep and D-sweep.  In each mode, four
different arrangement types are generated.  Three of them are "two-polytope"
arrangements, with polytopes formed from bounding planes tangent to two
spheres, each with uniformly-distributed tangent points.  The spheres' centers
are separated by three different distances: (a) zero, (b) a distance which
leads to intersection for approximately 50% of the arrangements, and (c) a
distance which separates the arrangements in all but rare cases.  The fourth
type of arrangement is one that ensures the half-spaces have an unbounded
intersection.

In the N-sweep mode, the total number of bounding planes N (divided equally
between the two polytopes) is swept through the values 10, 20, 50, 100, 200,
500, 1000, 2000, 5000, and 10000.  In this case the number of spatial dimen-
sions D is fixed at D = 3.  The generalized version of the homogeneous closest
point algorithm is timed, as well as the D = 3 specialized version.

In the D-sweep mode, the total number of bounding planes is fixed at 1000, and
the number of spatial dimensions D is swept through the values 2, 3, 5, 7, 10,
15, 20, 30, 50, 70, and 100.  In this case the only the generalized version of
the homogeneous closest point algorithm is used.

One other algorithm is tested as well for compaarison, Seidel's Linear
Programming algorithm (see ../common/seidel.h).

For each arrangement, two random objective vectors are used.  One represents a
point (for a closest-point calculation), and the other a direction (for a
linear programming problem).  (Note, for Seidel's LP algorithm, only direction
vector objectives are used.)

For each arrangement type, N and D value, and objective vector type, timing
measurements are made and statistics calculated.  Statistics are printed to
the standard output, as well as to separate files for each arrangement type if
the user supplies a directory name for these files (see the usage notes
below).  The filenames are of the form 'N_<TYPE>' and 'D_<TYPE>', where <TYPE>
is 'OVL', 'MRG', 'SEP', and 'UNB'.  These denote the four types of arrange-
ments described above: overlapping, marginal, separated, and unbounded.

The statistics output is given in pairs of columns, each pair giving the
median and median absolute deviation for the N or D value in the corresponding
row.  The headers for each column pair are are of the form '<ALG>/<OBJ>' and
'd(<ALG>/<OBJ>)'.  The values of <ALG> are 'HCP', 'HCP3' (for N sweeps only),
and 'Seid'.  These correspond to homogeneous closest point, D = 3 specialized
HCP, and Seidel's LP algorithm, respectively.  The values of <OBJ> are 'LP'
and (for HCP and HCP3), 'CP'.  These correspond to the Linear Programming
(direction) and Closest Point (point) objective vector types.

The individual measurements may be recorded into separate files (one file for
each statistic calculated), if the user supplies a directory name for those
files (again, see the usage notes below).  The filename format is
't_{N,D}SWEEP_<TYPE>_<ALG>_<OBJ>_N<CNT>_D<DIM>'.  Here, <TYPE>, <ALG>, and
<OBJ> have the same meaning as described above, <CNT> is the total number of
half-spaces in the arrangements measured, and <DIM> is the number of spatial
dimensions.

The units of all timing values and statistics are microseconds.

## Building

Run make in the current folder.  The executable perf.exe will be created in
the directory ../bin/

## Running

Execute the command

../bin/perf.exe type [m M] [g G] [x X] [o o_dir] [a a_dir] [v verbosity]

  type = 'N', 'D', or a test ID.  If 'N', the total number of half-spaces
    is swept through {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000} while
    the number of spatial dimensions is fixed at D = 3.  If 'D', the number of
    spatial dimensions is swept through {2, 3, 5, 7, 10, 15, 20, 30, 50, 70,
    100} while the total number of half-spaces is fixed at N = 1000.
    Otherwise, this argument needs to be a test ID.  A test ID is output when
    an error occurs during a run. It will specify a test, objective, D, N, and
    arrangement to run.

  M = number of arrangments per test group.  If 0, M is set to (D,N)-dependent
    default values.

  G = number of test groups for calculating statistics.  If 0, G is set to
    (D,N)-dependent default values.

  X = 0 or 1.  If 1, only HCP tests are performed.  Default = 0.

  o_dir = directory to write statistics files.  If not given, no files will be.
    written (stats are written to std out in any case).

  a_dir = directory to write files containing all measurments.  If not given,
    no such files will be written.

  verbosity = 0, 1, or 2. How much progress text is output.  Default = 1.

Use

../bin/perf.exe ?

to output help.
