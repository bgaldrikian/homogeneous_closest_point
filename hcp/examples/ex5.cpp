// Separation example.  Attempt to find the closest point to the intersection of two shapes that have no intersection.
// HCP returns 0 to indicate there is no intersection.  The returned point is meaningless.

#include "src/general/hcp.h"
#include "src/general/sets/hcp_static_polytope.h"
#include "src/general/sets/hcp_shape_pair.h"

#include <stdio.h>

#define SQRT1_3  ((real)0.577350269189625764509)  // 1/sqrt(3)

int main()
{
	printf("Separation example.  Attempt to find the closest point to the origin in the\n");
	printf("intersection of a tetrhedron and a cube, however the shapes have no\n");
	printf("intersection.\n");
	printf("\n");
	printf("The tetrahedron has sides bounded by the planes :\n");
	printf("  yz plane\n");
	printf("  zx plane\n");
	printf("  xy plane\n");
	printf("  the plane with normal and point (1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3))\n");
	printf("\n");
	printf("The cube is axis-aligned with unit width, and centered at (7/4, 7/4, 0).\n");
	printf("\n");
	printf("The correct query return value is 0, signifying an infeasible set.\n");

	const real tetrahedron_planes[4][4] = { { -1, 0, 0, 0 }, { 0,-1, 0, 0 }, { 0, 0,-1, 0 }, { SQRT1_3, SQRT1_3, SQRT1_3,-1 } };
	HCP_Static_Polytope tetrahedron(3, 4, tetrahedron_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes

	const real cube_planes[6][4] = { { -1, 0, 0, 1.25 }, { 1, 0, 0,-2.25 }, { 0,-1, 0, 1.25 }, { 0, 1, 0,-2.25 }, { 0, 0,-1,-0.5 }, { 0, 0, 1,-0.5 } };
	HCP_Static_Polytope cube(3, 6, cube_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes

	real p[4] = { 0, 0, 0, 1 };	// Query point (0, 0, 0)
	const int result = HCPA(3).solve(HCP_Shape_Pair(3, &tetrahedron, &cube), p);

	printf("\nresult = %d\n\n", result);

	return 0;
}
