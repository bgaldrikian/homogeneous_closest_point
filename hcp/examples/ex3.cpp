// Unbounded support point example.  Finds a point p in the given open pyramid which maximizes p dot q.
// In this case there is no bounded maximum, so the point is at infinity.  This is represented by
// a direction vector with last (projective) component equal to 0.

#include "src/general/hcp.h"
#include "src/general/sets/hcp_static_polytope.h"

#include <stdio.h>

#define SQRT1_3  ((real)0.577350269189625764509)  // 1/sqrt(3)

int main()
{
	const real pyramid_planes[3][4] = { { 0,-1, 0, 0 }, { 0, 0,-1, 0 }, { SQRT1_3, SQRT1_3, SQRT1_3,-1 } };
	HCP_Static_Polytope pyramid(3, 3, pyramid_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes

	real p[4] = {-1, 1,-1, 0 };	// Query direction (-1, 1,-1)
	const int result = HCPA(3).solve(pyramid, p);

	printf("result = %d, support direction = (%g %g %g %g)\n", result, p[0], p[1], p[2], p[3]);

	return 0;
}
