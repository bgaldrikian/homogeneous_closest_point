// Closest-point example, specialized for D = 3.  Finds the closest point on a tetrahedron to a query point.

#include "src/specialized/hcp3d.h"
#include "src/specialized/sets/hcpd_static_polytope.h"

#include <stdio.h>

#define SQRT1_3  ((real)0.577350269189625764509)  // 1/sqrt(3)

int main()
{
	const real tetrahedron_planes[4][4] = { { -1, 0, 0, 0 },{ 0,-1, 0, 0 },{ 0, 0,-1, 0 },{ SQRT1_3, SQRT1_3, SQRT1_3,-1 } };
	HCPD_Static_Polytope<3> tetrahedron(4, tetrahedron_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes

	real p[4] = { 1, 2, 3, 1 };	// Query point (1, 2, 3)
	const int result = hcpd_solve<3>(tetrahedron, p);

	printf("result = %d, closest point = (%g %g %g %g)\n", result, p[0], p[1], p[2], p[3]);

	return 0;
}
