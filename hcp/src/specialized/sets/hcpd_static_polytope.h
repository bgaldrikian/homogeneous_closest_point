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

#ifndef _HCPD_STATIC_POLYTOPE_H_
#define _HCPD_STATIC_POLYTOPE_H_


#include "src/hcp_defs.h"
#include "src/specialized/hcpd_math.h"

#include <limits>
#include <stddef.h>


/**
Implementation of HCP_Halfspace_Set, an explicit representation of a finite number of half-spaces.

The shape represented is the interesection of the half-spaces.

Written for general D, but functions used within might only be specialized for particular values of D.
*/
template<uint D>
struct HCPD_Static_Polytope : HCP_Halfspace_Set
{
	size_t		halfspace_count;	//!< The number of half-spaes in the halfspaces array.
	const real*	halfspaces;			//!< Points to an array of groups of five floats, each representing the bounding plane of a half-space with outward-pointing normal.  Array size must be a multiple of D+1 times halfspace_count.

	/** Default ctor creates a blank object. */
	HCPD_Static_Polytope() : halfspace_count(0), halfspaces(nullptr)	{}

	/** ctor which initializes the array of half-spaces. */
	HCPD_Static_Polytope(size_t in_halfspace_count, const real* in_halfspaces) : halfspace_count(in_halfspace_count), halfspaces(in_halfspaces) {}

	/** Implementation of HCP_Halfspace_Set interface. */
	real
	farthest_halfspace(real* plane, const real* point) const override
	{
		real greatest_s = -std::numeric_limits<real>().max();
		const real* halfspace = halfspaces;
		const real* stop = halfspace + (D + 1)*halfspace_count;
		for (const real* test = halfspaces; test < stop; test += D+1)
		{
			const real s = dot<D + 1>(point, test);
			if (s > greatest_s)
			{
				greatest_s = s;
				halfspace = test;
			}
		}

		// Return results
		if (halfspace < stop) cpy<D + 1>(plane, halfspace);
		else { zero<D + 1>(plane); plane[D] = -(real)1; }	// Empty set represents all of space, so the plane is at infinity with outward-pointing normal

		return greatest_s;
	}
};


#endif // #ifndef _HCPD_STATIC_POLYTOPE_H_
