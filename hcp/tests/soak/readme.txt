# soak.exe

This is a soak test for the homogeneous closest point algorithm meant to test
robustness.  It creates random arrangements of half-spaces which represent two
convex polytopes in chosen number of dimensions.  The number of dimensions D
can be any positive integer, however the O(D^2) memory requirements of HCP
will make execution unfeasible for large enough D.

There are two arrangement types: bounded and unbounded.  In the bounded mode,
bounding planes are chosen to be tangent to unit spheres, with tangent points
distributed uniformly on the surfaces of the spheres.  Note, this does not
guarantee a bounded polytope, but it makes boundedness very likely when the
number of planes is large enough.  The spheres' centers are separated by a
distance which gives intersection in approximately 50% of the arrangements,
when each polytope has 50 faces, for D <= 10.  These separations are arrived
at numerically using a separate program.

In the unbounded mode, a single polytope is created such that a ray with
endpoint at the origin, and direction given by a randomly-chosen direction,
is contained within the polytope.

There are also two different query types: closest point and support point.
In the closest point mode, a random query point is chosen within a sphere of
radius 10.  In support point mode, a random direction is chosen (uniformly
distributed).

Finally, there are two versions of the algorithm which can be chosen.  One is
the generalized algorithm which in theory works for any dimensionality D
(with the memory caveat given above).  The other is the dimension-specialized
version, defined for D = 1, 2, 3, and 4.

If the general algorithm is used, then the result is compared with Seidel's
Linear Programming algorithm (see ../common/seidel.h).  If a discrepancy is
found that is outside of numerical tolerance, it is reported.

If a D-specialized algorithm is used, the result is compared with the general
algorithm.  Again, if a discrepancy is found that is outside of numerical
tolerance, it is reported.

When an error is reported, the trial index is given.  This index may be used
as an optional input parameter (see the usage below), causing the program to
reproduce the arrangement for that trial only.  This way the user may debug
error and failure cases.

The trials are grouped and a progress indicator shows intersection statistics
at the end of each group.  The number of trials per group and the number of
groups is selectable by the user.  The default is one million trials per group
and one thousand groups.

## Building

Run make in the current folder.  The executable soak.exe will be created in
the directory ../bin/

## Running

Execute the command

../bin/soak.exe dim [s S] [w W] [u U] [n N] [r R] [t trialIndex] [m M] [g G] [v V] [outfile]

	dim = the dimensionality of the test.
	S = 0 or 1.  Set to 1 to test against dimension-specialized hcp.  Default is 0.
	W = 0 (linear optimization) or 1 (closest point test).  Default is 1.
	U = 0 or 1, determines if all intersections will be unbounded.  Default is 0.
	N = total number of halfspaces in tests.  Default is 100.
	R = override offset to use.  If 0, uses pre-calculated offsets.  Default is 0.
	trialIndex = the single trial number to run, if given.
	M = number of groups of trials to perform.  Default is 1000.
	G = number of trials per group.  Default is 1000000.
	V = verbosity.  0 = result only, 1 = progress output.  Default is 1.
	outfile = an optional output filename.  If none is given, the standard output is used.

Use

../bin/soak.exe ?

(or no arguments) to output help.
