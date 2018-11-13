// Support point example.  Finds a point p in the given tetrahedron which maximizes p dot q.

#include "src/general/hcp.h"
#include "src/general/sets/hcp_static_polytope.h"

#include <stdio.h>

#define SQRT1_3  ((real)0.577350269189625764509)  // 1/sqrt(3)

int main()
{
	printf("Support point example.  Finds a point p in the given tetrahedron which\n");
	printf("maximizes p dot q, where q is the query direction.\n");
	printf("\n");
	printf("The tetrahedron has sides bounded by the planes :\n");
	printf("  yz plane\n");
	printf("  zx plane\n");
	printf("  xy plane\n");
	printf("  the plane with normal and point (1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3))\n");
	printf("\n");
	printf("The query direction q is (-1, 1,-1).\n");
	printf("\n");
	printf("The correct query return value is 1, signifying a feasible set.  The exact\n");
	printf("maximal point is (0, sqrt(3), 0).\n");

	const real tetrahedron_planes[4][4] = { { -1, 0, 0, 0 }, { 0,-1, 0, 0 }, { 0, 0,-1, 0 }, { SQRT1_3, SQRT1_3, SQRT1_3,-1 } };
	HCP_Static_Polytope tetrahedron(3, 4, tetrahedron_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes

	real p[4] = {-1, 1,-1, 0 };	// Query direction (-1, 1,-1)
	const int result = HCPA(3).solve(tetrahedron, p);

	printf("\nresult = %d, support point = (%g %g %g %g)\n\n", result, p[0], p[1], p[2], p[3]);

	return 0;
}
