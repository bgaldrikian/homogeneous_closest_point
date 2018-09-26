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

#ifndef _HCP_TOI_H_
#define _HCP_TOI_H_


/** Time of Intersection function using the homogeneous closest point algorithm. */

#include "src/general/hcp.h"
#include "src/general/sets/hcp_shape_pair.h"


/**
Find the Time of Intersection of two shapes defined by half-space sets, with given linear velocities, using the homogeneous closest point algorithm.

\param[out]		toi				if this function returns 1 then toi is filled with the time of intersection or separation (depending on if 'intersection' is true or false, respectively).
\param[in]		D				the number of spatial dimensions.
\param[in]		shape1			a user-defined implementation of HCP_Halfspace_Set.
\param[in]		v1				points to an array of D real values representing the world velocity of shape1.
\param[in]		shape2			a user-defined implementation of HCP_Halfspace_Set.
\param[in]		v2				points to an array of D real values representing the world velocity of shape2.
\param[in]		scratch			user-provided scratch space.  It must point to hcp_toi_scratch_required(D) bytes.
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
int hcp_toi
(
	real& toi,
	uint D,
	const HCP_Halfspace_Set& shape1,
	const real* v1,
	const HCP_Halfspace_Set& shape2,
	const real* v2,
	void* scratch,
	bool intersection,
	real* p = 0,
	real* rho = 0,
	real tol = hcp_default_tolerance()
);


/** Scratch space required for hcp_toi for a given D. */
inline size_t hcp_toi_scratch_required(uint D) { return HCP_Shape_Pair_NA::storage_required(D + 1) + (D + 2) * sizeof(real) + hcp_solve_scratch_required(D + 1); }


/**
This class calls hcp_toi with storage that it manages.
*/
class HCPA_TOI
{
public:
	/**
	ctor

	\param[in]	D	the number of spatial dimensions.
	*/
	HCPA_TOI(uint D) : m_D(D), m_storage(new char[hcp_toi_scratch_required(D)]) {}
	~HCPA_TOI() { delete[] m_storage; }

	/**
	Find the Time of Intersection of two shapes defined by half-space sets, with given linear velocities, using the homogeneous closest point algorithm.

	\param[out]	toi				if this function returns 1 then toi is filled with the time of intersection or separation (depending on if 'intersection' is true or false, respectively).
	\param[in]	shape1			a user-defined implementation of HCP_Halfspace_Set.
	\param[in]	v1				points to an array of D real values representing the world velocity of shape1.
	\param[in]	shape2			a user-defined implementation of HCP_Halfspace_Set.
	\param[in]	v2				points to an array of D real values representing the world velocity of shape2.
	\param[in]	intersection	if true, this function calculates time of intersection.  Otherwise it calculates time of separation.
	\param[out]	p				(optional) if not NULL, and if this function returns 1, then p[0...D] is filled with the homogeneous position where intersection or separation occurs at the calculated toi.
								Default is NULL.
	\param[out]	rho				(optional) if not NULL, *rho will be set to the absolute value of the determinant of the linear system solved for the returned solution.  0 < rho <= 1.  Default is NULL.
								NOTE: a crude approximation for the condition number is (1+sqrt(1-rho))/rho.
	\param[in]	tol				(optional) the requested solution accuracy, relative to the size of the output vector.

	Return values:
		 1: there is a time when the sets shape1 and shape2 intersect
		 0: there is no time when the sets shape1 and shape2 intersect
		-1: the function failed to find a result
	*/
	int
	toi(real& toi, const HCP_Halfspace_Set& shape1, const real* v1, const HCP_Halfspace_Set& shape2, const real* v2, bool intersection, real* p = 0, real* rho = 0, real tol = hcp_default_tolerance())
	{
		return hcp_toi(toi, m_D, shape1, v1, shape2, v2, m_storage, intersection, p, rho, tol);
	}

	/**
	\return the number of spatial dimensions with which this object is set (and allocated) to perform calculations.
	*/
	uint	get_D() const { return m_D; }

private:
	uint	m_D;		//!< The number of spatial dimensions
	char*	m_storage;	//!< A pointer to the allocated storage space
};


#endif // #ifndef _HCP_TOI_H_
