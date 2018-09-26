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

#ifndef _HCP_DEFS_H_
#define _HCP_DEFS_H_


/*** Public definitions for the homogeneous closest point algorithm. */

#include "types.h"
#include <limits>


/** Floating point tolerance for homogeneous closest point algorithm. */
inline constexpr real	hcp_default_tolerance() { return 8 * std::numeric_limits<real>().epsilon(); }


/** Abstraction of a set of half-spaces used by the homogeneous closest point algorithm. */
struct HCP_Halfspace_Set
{
	/**
	The implementation of farthest_halfspace should return the half-space "most below" the given point in D spatial
	dimensions.  The value of D is *implicitly defined* by the implementation.  The point is represented by a vector
	in projective coordinates, and its last element, point[D], will not necessarily equal 1.  However, point[D] will
	be non-negative.  If point[D] is zero, the vector is interpreted as a direction, or a point at infinity.

	The plane returned is the boundary of the half-space found (with outward-pointing normal), and is also
	represented as a vector in projective coordinates (the coefficients of the plane equation).  So the plane vector
	returned should have the greatest dot product with the input point, of any bounding plane for the half-spaces
	in the set.

	\param[out]	plane	the returned half-space boundary, a (D+1)-vector.
	\param[in]	point	the input point, a (D+1)-vector.

	\return the dot product of the point and returned plane.
	*/
	virtual	real	farthest_halfspace(real* plane, const real* point) const = 0;
};


#endif // #ifndef _HCP_DEFS_H_
