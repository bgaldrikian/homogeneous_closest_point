// Copyright (c) 2014-2018 NVIDIA Corporation
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.


#include "util/unit_test.h"

#include "src/general/hcp.h"
#include "math/linalg.h"

#include "src/general/hcp_toi.h"

#include "src/general/sets/hcp_static_polytope.h"
#include "src/general/sets/hcp_polytope.h"
#include "src/general/sets/hcp_shape_pair.h"

#include "src/specialized/hcpd_toi.h"
#include "src/specialized/hcp2d.h"
#include "src/specialized/hcp3d.h"
#include "src/specialized/hcp4d.h"
#include "src/specialized/sets/hcpd_static_polytope.h"

#include "math/dist.h"

#include <vector>
#include <cmath>
#include <algorithm>


#define SQRT1_2  ((real)0.707106781186547524401)  // 1/sqrt(2)
#define SQRT1_3  ((real)0.577350269189625764509)  // 1/sqrt(3)


static const real s_void2d[3][3] =
{
	{-1, 0, 0 },
	{ 0,-1, 0 },
	{ SQRT1_2, SQRT1_2, 1 }
};


static const real s_triangle[3][3] =
{
	{-1, 0, 0 },
	{ 0,-1, 0 },
	{ SQRT1_2, SQRT1_2,-1 }
};


static const real s_square[4][3] =
{
	{-1, 0,-1 },
	{ 1, 0,-1 },
	{ 0,-1,-1 },
	{ 0, 1,-1 }
};


static const real s_square_translated_5_5[4][3] =
{
	{-1, 0, 4 },
	{ 1, 0,-6 },
	{ 0,-1, 4 },
	{ 0, 1,-6 }
};


static const real s_void3d[4][4] =
{
	{-1, 0, 0, 0 },
	{ 0,-1, 0, 0 },
	{ 0, 0,-1, 0 },
	{ SQRT1_3, SQRT1_3, SQRT1_3, 1 }
};


static const real s_tetrahedron[4][4] =
{
	{-1, 0, 0, 0 },
	{ 0,-1, 0, 0 },
	{ 0, 0,-1, 0 },
	{ SQRT1_3, SQRT1_3, SQRT1_3,-1 }
};


static const real s_cube[6][4] =
{
	{-1, 0, 0,-1 },
	{ 1, 0, 0,-1 },
	{ 0,-1, 0,-1 },
	{ 0, 1, 0,-1 },
	{ 0, 0,-1,-1 },
	{ 0, 0, 1,-1 }
};


static const real s_cube_translated_5_5_5[6][4] =
{
	{-1, 0, 0, 4 },
	{ 1, 0, 0,-6 },
	{ 0,-1, 0, 4 },
	{ 0, 1, 0,-6 },
	{ 0, 0,-1, 4 },
	{ 0, 0, 1,-6 }
};


UNIT_TEST(HCPLP_degenerate_void_2D)
{
	const real planes[2][3] = {{ 0, 1, 1 }, { 0,-1, 1 }};
	const real q[3] = { 1, 0, 0 };

	const HCP_Static_Polytope halfspace_set(2, 2, planes[0]);
	real p[3];
	la::cpy(p, q, 3);
	
	const int CP_result = HCPA(2).solve(halfspace_set);
	const int LP_result = HCPA(2).solve(halfspace_set, p);

	EXPECT(CP_result == 0);
	EXPECT(LP_result == 1);
	EXPECT(la::diff_norm_2_sq(p, q, 3) == 0);
}


UNIT_TEST(HCP_3D_dynamic_spatial_polytopes)
{
	// Create two HCP_Polytopes to represent the sets of half-spaces.
	HCP_Polytope poly1(3), poly2(3);

	// Point HCP_Polytopes to the arrays of bounding_planes
	poly1.set_halfspaces(6, s_cube[0]);
	poly2.set_halfspaces(4, s_tetrahedron[0]);

	// Using a temporal halfspace set for a static query.
	HCP_Shape_Pair halfspaces(3, &poly1, &poly2);

	// These shapes are overlapping.
	std::vector<char> scratch(hcp_solve_scratch_required(3));
	const int result_A = hcp_solve(3, halfspaces, scratch.data());
	EXPECT(result_A == 1);

	// Translate the cube, now the shapes are disjoint
	const real t[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 5, 5, 5, 1 };
	poly1.set_local_to_world_transform(3, t);
	const int result_B = hcp_solve(3, halfspaces, scratch.data());
	EXPECT(result_B == 0);
}


UNIT_TEST(HCP_temporal_and_spatial_2D_polytopes)
{
	// Create two HCP_Static_Polytopes to represent the sets of half-spaces.
	HCP_Static_Polytope poly1(2), poly2(2);

	// Point HCP_Static_Polytopes to the arrays of bounding_planes
	poly1.set_halfspaces(4, s_square_translated_5_5[0]);
	poly2.set_halfspaces(3, s_triangle[0]);

	// Give the square a velocity towards the triangle
	const real v1[2] = {-1,-1 };
	const real v2[2] = { 0, 0 };

	std::vector<char> scratch(std::max(hcp_toi_scratch_required(2), hcp_solve_scratch_required(2)));
	real toi;
	const int result_toi = hcp_toi(toi, 2, poly1, v1, poly2, v2, scratch.data(), true);
	EXPECT(result_toi == 1);

	const real expected_toi = 4 - SQRT1_2;
	const real error = toi - expected_toi;
	EXPECT(std::abs(error) <= hcp_default_tolerance()*expected_toi);

	const int result_find = hcp_solve(2, HCP_Shape_Pair(2, &poly1, &poly2), scratch.data());
	EXPECT(result_find == 0);
}


UNIT_TEST(HCP_temporal_and_spatial_3D_polytopes)
{
	// Create two HCP_Static_Polytopes to represent the sets of half-spaces.
	HCP_Static_Polytope poly1(3), poly2(3);

	// Point HCP_Static_Polytopes to the arrays of bounding_planes
	poly1.set_halfspaces(6, s_cube_translated_5_5_5[0]);
	poly2.set_halfspaces(4, s_tetrahedron[0]);

	// Give the cube a velocity towards the tetrahedron
	const real v1[3] = {-1,-1,-1 };
	const real v2[3] = { 0, 0, 0 };

	std::vector<char> scratch(std::max(hcp_toi_scratch_required(3), hcp_solve_scratch_required(3)));
	real toi;
	const int result_toi = hcp_toi(toi, 3, poly1, v1, poly2, v2, scratch.data(), true);
	EXPECT(result_toi == 1);

	const real expected_toi = 4 - SQRT1_3;
	const real error = toi - expected_toi;
	EXPECT(std::abs(error) <= hcp_default_tolerance()*expected_toi);

	const int result_find = hcp_solve(3, HCP_Shape_Pair(3, &poly1, &poly2), scratch.data());
	EXPECT(result_find == 0);
}


UNIT_TEST(HCP_temporal_and_spatial_2D_polytopes_specialized)
{
	// Point HCP2D_Static_Polytopes to the arrays of bounding_planes
	HCPD_Static_Polytope<2> poly1(4, s_square_translated_5_5[0]);
	HCPD_Static_Polytope<2> poly2(3, s_triangle[0]);

	// Give the square a velocity towards the triangle
	const real v1[2] = {-1,-1 };
	const real v2[2] = { 0, 0 };

	real toi;
	const int result_toi = hcpd_toi<2>(toi, poly1, v1, poly2, v2, true);
	EXPECT(result_toi == 1);

	const real expected_toi = 4 - SQRT1_2;
	const real error = toi - expected_toi;
	EXPECT(std::abs(error) <= hcp_default_tolerance()*expected_toi);

	const int result_find = hcpd_solve<2>(HCPD_Shape_Pair<2>(&poly1, &poly2));
	EXPECT(result_find == 0);
}


UNIT_TEST(HCP_temporal_and_spatial_3D_polytopes_specialized)
{
	// Point HCP2D_Static_Polytopes to the arrays of bounding_planes
	HCPD_Static_Polytope<3> poly1(6, s_cube_translated_5_5_5[0]);
	HCPD_Static_Polytope<3> poly2(4, s_tetrahedron[0]);

	// Give the cube a velocity towards the tetrahedron
	const real v1[3] = {-1,-1,-1 };
	const real v2[3] = { 0, 0, 0 };

	real toi;
	const int result_toi = hcpd_toi<3>(toi, poly1, v1, poly2, v2, true);
	EXPECT(result_toi == 1);

	const real expected_toi = 4 - SQRT1_3;
	const real error = toi - expected_toi;
	EXPECT(std::abs(error) <= hcp_default_tolerance()*expected_toi);

	const int result_find = hcpd_solve<3>(HCPD_Shape_Pair<3>(&poly1, &poly2));
	EXPECT(result_find == 0);
}


UNIT_TEST(HCP_toi_3D_dynamic_polytopes)
{
	// Create two HCP_Polytopes to represent the sets of half-spaces.
	HCP_Polytope poly1(3), poly2(3);

	// Point HCP_Polytopes to the arrays of bounding_planes
	poly1.set_halfspaces(6, s_cube[0]);
	poly2.set_halfspaces(4, s_tetrahedron[0]);

	// Translate the cube and give it a velocity towards the tetrahedron
	const real t[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 5, 5, 5, 1 };
	poly1.set_local_to_world_transform(3, t);
	const real v1[3] = {-1,-1,-1 };
	const real v2[3] = { 0, 0, 0 };

	std::vector<char> scratch(hcp_toi_scratch_required(3));
	real toi;
	const int result = hcp_toi(toi, 3, poly1, v1, poly2, v2, scratch.data(), true);
	EXPECT(result == 1);

	const real expected_toi = 4  - SQRT1_3;
	const real error = toi - expected_toi;
	EXPECT(std::abs(error) <= hcp_default_tolerance()*expected_toi);
}


UNIT_TEST(HCP_radius_estimate_3D)
{
	const real separation = 100;
	const uint D = 3;
	const uint N = 100;

	real v1[D];
	la::zero(v1, D);
	v1[0] = -1;
	real v2[D];
	la::zero(v2, D);
	v2[0] = 1;

	std::vector<char> hcp_toi_scratch(hcp_toi_scratch_required(D));

	Random<real> rnd;

	std::vector<real> planes(2 * N*(D + 1));
	HCP_Static_Polytope poly1(D, N, &planes[0]);
	HCP_Static_Polytope poly2(D, N, &planes[N*(D + 1)]);

	const uint trialCount = 10000;
	std::vector<real> R(trialCount);
	uint numR = 0;
	real sumR = 0;

	for (uint trial = 0; trial < trialCount; ++trial)
	{
		for (uint i = 0; i < 2 * N; ++i)
		{
			rnd.sphere(&planes[i*(D + 1)], D);
			planes[i*(D + 1) + D] = -1;
		}

		for (uint i = 0; i < N; ++i)
			planes[i*(D + 1) + D] -= planes[i*(D + 1)] * separation;
		for (uint i = N; i < 2 * N; ++i)
			planes[i*(D + 1) + D] += planes[i*(D + 1)] * separation;

		real toi;
		const int result = hcp_toi(toi, D, poly1, v1, poly2, v2, hcp_toi_scratch.data(), true);

		EXPECT(result == 1);
		R[numR] = separation - toi;
		sumR += R[numR++];
	}

	const real meanR = sumR / numR;
	real varR = 0;
	for (uint i = 0; i < numR; ++i) varR += pow(R[i] - meanR, 2);
	const real sigmaR = sqrt(varR / (numR - 1));

	EXPECT_NEAR((real)1.0315, meanR, (real)0.005);
	EXPECT_NEAR((real)0.0189, sigmaR, (real)0.0002);
}
