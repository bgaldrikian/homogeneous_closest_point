// Copyright (c) 2017-2018 NVIDIA Corporation
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


/*** Implementation of Time of Intersection using the homogeneous closest point algorithm ***/

#include "src/general/hcp_toi.h"
#include "src/general/hcp.h"
#include "math/linalg.h"

#include <limits>


/*
Uses a spatial HCP_Halfspace_Set implementation for its underlying query.

The vector format is:

     / \
    | t |  <=  1-dimensional temporal component
v = | x |  <=  D-dimensional spatial components
    | w |  <=  1-dimensional projective component
     \ /
*/
struct HCP_Temporal_Halfspace_Set : public HCP_Halfspace_Set
{
	HCP_Temporal_Halfspace_Set(uint D, const HCP_Halfspace_Set& halfspace_set, const real* lin_vel) : m_D(D), m_halfspace_set(halfspace_set), m_lin_vel(lin_vel) {}

	/* Implementation of HCP_Halfspace_Set interface */
	real
	farthest_halfspace(real* plane, const real* point) const override
	{
		la::madd(plane, m_lin_vel, -point[0], point + 1, m_D);	// Transform into moving frame of reference.  Using plane as temporary storage
		plane[m_D] = point[m_D + 1];	// Create spatial homogeneous vector

		// Arguments of farthest_halfspace_internal must be able to point to same memory.  Now plane will refer to a plane equation
		const real s = m_halfspace_set.farthest_halfspace(plane + 1, plane);
		plane[0] = -la::dot(plane + 1, m_lin_vel, m_D);	// Set plane d_dot

		la::normal_proj_plane(plane, plane, m_D + 1);	// Normalize plane
		return la::dot(plane, point, m_D + 2);
	}

protected:
	uint						m_D;
	const HCP_Halfspace_Set&	m_halfspace_set;
	const real*					m_lin_vel;
};


int
hcp_toi(real& toi, uint D, const HCP_Halfspace_Set& shape1, const real* v1, const HCP_Halfspace_Set& shape2, const real* v2, void* scratch, bool intersection, real* p, real* rho, real tol)
{
	char* shape_pair_scratch = (char*)scratch;
	real* q = (real*)(shape_pair_scratch + HCP_Shape_Pair_NA::storage_required(D + 1));
	scratch = q + (D + 2);

	const real q_t = intersection ? -(real)1 : (real)1;
	la::zero(q, D + 2);
	q[0] = q_t;

	HCP_Temporal_Halfspace_Set shape1_t(D, shape1, v1);
	HCP_Temporal_Halfspace_Set shape2_t(D, shape2, v2);
	const int result = hcp_solve(D + 1, HCP_Shape_Pair_NA(D + 1, shape_pair_scratch, &shape1_t, &shape2_t), scratch, q, rho, tol);

	if (result == 1)
	{
		toi = q[D + 1] != 0 ? q[0] / q[D + 1] : q_t*std::numeric_limits<real>().max();
		if (p) la::cpy(p, q + 1, D + 1);
	}

	return result;
}


