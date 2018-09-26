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

#ifndef _HCPD_SHAPE_PAIR_H_
#define _HCPD_SHAPE_PAIR_H_

#include "src/hcp_defs.h"
#include "src/specialized/hcpd_math.h"


/**
Implementation of HCP_Halfspace_Set, for fixed D, representing the union of two other half-space sets.

The shape defined is the intersection of the shapes defined by the two other half-space sets.
*/
template<unsigned D>
struct HCPD_Shape_Pair : HCP_Halfspace_Set
{
	HCP_Halfspace_Set* shape1;	//!< A pointer to a user-defined implementation of HCP_Halfspace_Set.
	HCP_Halfspace_Set* shape2;	//!< A pointer to a user-defined implementation of HCP_Halfspace_Set.

	/** Default ctor creates a blank object. */
	HCPD_Shape_Pair() : shape1(nullptr), shape2(nullptr) {}

	/** ctor which initializes the shapes. */
	HCPD_Shape_Pair(HCP_Halfspace_Set* in_shape1, HCP_Halfspace_Set* in_shape2) : shape1(in_shape1), shape2(in_shape2) {}

	/** Implementation of HCP_Halfspace_Set interface. */
	real
	farthest_halfspace(real* plane, const real* point) const override
	{
		// Choose the farthest halfspace from either set
		const real s1 = shape1->farthest_halfspace(plane, point);
		real plane2[D+1];
		const real s2 = shape2->farthest_halfspace(plane2, point);
		if (s2 < s1) return s1;
		cpy<D+1>(plane, plane2);
		return s2;
	}
};


#endif // #ifndef _HCPD_SHAPE_PAIR_H_
