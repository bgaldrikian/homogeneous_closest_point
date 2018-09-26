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

#ifndef _HCPD_TOI_H_
#define _HCPD_TOI_H_


#include "src/specialized/hcpd.h"
#include "src/specialized/sets/hcpd_shape_pair.h"

#include <cmath>
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
template<uint D>
struct HCPD_Temporal_Halfspace_Set : public HCP_Halfspace_Set
{
	HCPD_Temporal_Halfspace_Set(const HCP_Halfspace_Set& halfspace_set, const real* lin_vel) : m_halfspace_set(halfspace_set), m_lin_vel(lin_vel) {}

	/* Implementation of HCP_Halfspace_Set interface */
	real
	farthest_halfspace(real* plane, const real* point) const override
	{
		madd<D>(plane, m_lin_vel, -point[0], point + 1);	// Transform into moving frame of reference.  Using plane as temporary storage
		plane[D] = point[D + 1];	// Create spatial homogeneous vector

		// Arguments of farthest_halfspace_internal must be able to point to same memory.  Now plane will refer to a plane equation
		const real s = m_halfspace_set.farthest_halfspace(plane + 1, plane);
		plane[0] = -dot<D>(plane + 1, m_lin_vel);	// Set plane d_dot

		// Normalize plane
		const real plane_norm2_sq = dot<D+1>(plane, plane);
		div<D+2>(plane, plane, sqrt(plane_norm2_sq));

		return dot<D+2>(plane, point);
	}

protected:
	const HCP_Halfspace_Set&	m_halfspace_set;
	const real*					m_lin_vel;
};


/*
Find the Time of Intersection of two shapes defined by half-space sets, with given linear velocities, using a (D+1)-specific specialization of the homogeneous closest point algorithm.

The corresponding (D+1)-specific specialization of the homogeneous closest point algorithm must be defined by including the appropriate header (hcp1d.h, hcp2d.h, hcp3d.h, or hcp4d.h).

For example, to use hcpd_toi<3>, the file hcp4d.h must be included.

\param[out]		toi				if this function returns 1 then toi is filled with the time of intersection or separation (depending on if 'intersection' is true or false, respectively).
\param[in]		shape1			a user-defined implementation of HCP_Halfspace_Set.
\param[in]		v1				points to an array of D real values representing the world velocity of shape1.
\param[in]		shape2			a user-defined implementation of HCP_Halfspace_Set.
\param[in]		v2				points to an array of D real values representing the world velocity of shape2.
\param[in]		intersection	if true, this function calculates time of intersection.  Otherwise it calculates time of separation.
\param[in,out]	p				(optional) if not NULL, and if this function returns 1, then p[0...D] is filled with the homogeneous position where intersection or separation occurs at the calculated toi.
								Default is NULL.
\param[out]		rho				(optional) if not NULL, *rho will be set to the absolute value of the determinant of the linear system solved for the returned solution.  0 < rho <= 1.  Default is NULL.
								NOTE: a crude approximation for the condition number is (1+sqrt(1-rho))/rho.
\param[in]		tol				(optional) the requested solution accuracy, relative to the size of the output vector.

Return values:
	 1: there is a time when the sets shape1 and shape2 intersect
	 0: there is no time when the sets shape1 and shape2 intersect
	-1: the function failed to find a result
*/
template<uint D>
inline int
hcpd_toi(real& toi, const HCP_Halfspace_Set& shape1, const real* v1, const HCP_Halfspace_Set& shape2, const real* v2, bool intersection, real* p = 0, real* rho = 0, real tol = hcp_default_tolerance())
{
	real q[D + 2];

	const real q_t = intersection ? -(real)1 : (real)1;
	zero<D + 2>(q);
	q[0] = q_t;

	HCPD_Temporal_Halfspace_Set<D> shape1_t(shape1, v1);
	HCPD_Temporal_Halfspace_Set<D> shape2_t(shape2, v2);
	const int result = hcpd_solve<D+1>(HCPD_Shape_Pair<D+1>(&shape1_t, &shape2_t), q, rho, tol);

	if (result == 1)
	{
		toi = q[D + 1] != 0 ? q[0] : q_t*std::numeric_limits<real>().max();
		if (p) cpy<D+1>(p, q + 1);
	}

	return result;
}


#endif // #ifndef _HCPD_TOI_H_
