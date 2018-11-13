// Time of Intersection example.  The time when two moving shapes intersect.  The point of intersection is also returned.

#include "src/general/hcp_toi.h"
#include "src/general/sets/hcp_static_polytope.h"

#include <stdio.h>

#define SQRT1_3  ((real)0.577350269189625764509)  // 1/sqrt(3)

int main()
{
	printf("Time of Intersection example.  Calculates the time when a moving cube\n");
	printf("intersects a static tetrahedron, and returns the point of intersection.\n");
	printf("\n");
	printf("The tetrahedron has sides bounded by the planes :\n");
	printf("  yz plane\n");
	printf("  zx plane\n");
	printf("  xy plane\n");
	printf("  the plane with normal and point (1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3))\n");
	printf("\n");
	printf("The cube is axis-aligned with unit width, and centered at (5, 5, 5).  It has\n");
	printf("velocity (-1,-1,-1).\n");
	printf("\n");
	printf("The correct ToI return value is 1, signifying that there is a time when the\n");
	printf("sets intesect.  The exact time of intersection is 4 - 1/sqrt(3).  The exact\n");
	printf("point of intersection at that time is (1/sqrt(3), 1/sqrt(3), 1/sqrt(3)).\n");

	const real tetrahedron_planes[4][4] = { {-1, 0, 0, 0 }, { 0,-1, 0, 0 }, { 0, 0,-1, 0 }, { SQRT1_3, SQRT1_3, SQRT1_3,-1 } };
	HCP_Static_Polytope tetrahedron(3, 4, tetrahedron_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes
	const real tetrahedron_vel[3] = { 0, 0, 0 };

	const real cube_planes[6][4] = { {-1, 0, 0, 4 }, { 1, 0, 0,-6 }, { 0,-1, 0, 4 }, { 0, 1, 0,-6 }, { 0, 0,-1, 4 }, { 0, 0, 1,-6 } };
	HCP_Static_Polytope cube(3, 6, cube_planes[0]);	// Form a HCP_Halfspace_Set from bounding planes
	const real cube_vel[3] = {-1,-1,-1 };

	real toi;
	real p[4];
	const int result = HCPA_TOI(3).toi(toi, tetrahedron, tetrahedron_vel, cube, cube_vel, true, p);

	printf("\nresult = %d, ToI = %g, point of intersection = (%g %g %g %g)\n\n", result, toi, p[0], p[1], p[2], p[3]);

	return 0;
}
