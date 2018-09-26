// Separation example.  Attempt to find the closest point to the intersection of two shapes that have no intersection.
// HCP returns 0 to indicate there is no intersection.  The returned point is meaningless.

#include "src/general/hcp.h"
#include "src/general/sets/hcp_static_polytope.h"
#include "src/general/sets/hcp_shape_pair.h"

#include <stdio.h>

#define SQRT1_3  ((real)0.577350269189625764509)  // 1/sqrt(3)

int main()
{
	const real tetrahedron_planes[4][4] = { { -1, 0, 0, 0 }, { 0,-1, 0, 0 }, { 0, 0,-1, 0 }, { SQRT1_3, SQRT1_3, SQRT1_3,-1 } };
	HCP_Static_Polytope tetrahedron(3, 4, tetrahedron_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes

	const real cube_planes[6][4] = { { -1, 0, 0, 1.25 }, { 1, 0, 0,-2.25 }, { 0,-1, 0, 1.25 }, { 0, 1, 0,-2.25 }, { 0, 0,-1,-0.5 }, { 0, 0, 1,-0.5 } };
	HCP_Static_Polytope cube(3, 6, cube_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes

	real p[4] = { 0, 0, 0, 1 };	// Query point (0, 0, 0)
	const int result = HCPA(3).solve(HCP_Shape_Pair(3, &tetrahedron, &cube), p);

	printf("result = %d\n", result);

	return 0;
}
