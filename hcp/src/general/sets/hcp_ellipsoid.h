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

#ifndef _HCP_ELLIPSOID_H_
#define _HCP_ELLIPSOID_H_


/**
Implementation of HCP_Halfspace_Set, an implicit representation of an infinite number of half-spaces
bounding an ellipsoid, the ellipsoid being the interesection of the half-spaces.

This class derives from HCP_Mutable, giving it a transformation matrix that is effectivly applied
to a D-dimensional ball during the execution of farthest_halfspace.  Actually, this is transformation
is only approximately applied if it contains a non-uniform scale.  In that case, the returned
half-space may not be the one farthest from the given point when transformed.  However, upon iteration
of the homogeneous closest point algorithm, the planes returned by this function will lead to the
correct result.
*/


#include "src/hcp_defs.h"
#include "src/general/sets/hcp_mutable.h"
#include "src/general/sets/hcp_allocating_halfspace_set.h"


/**
This class does not perform any allocation, and requires storage be available it when it performs a calculation.

(The suffix "NA" stands for "Non-Allocating.")
*/
struct HCP_Ellipsoid_NA : public HCP_Mutable, public virtual HCP_Halfspace_Set
{
	/** ctor creates a blank object. */
	HCP_Ellipsoid_NA() : m_D(0), m_plane(nullptr) {}

	/** Implementation of HCP_Halfspace_Set interface. */
	real
	farthest_halfspace(real* plane, const real* point) const override
	{
		la::mul_t(m_plane, point, m_world_to_local_t, m_D + 1, m_D + 1);	// Transform spatial vector to local space, this will become the local plane direction
		real n = la::normal(m_plane, m_plane, m_D);
		if (n == 0) m_plane[0] = n = 1;	// Create plane normal
		m_plane[m_D] = -1;	// Circle displacement
		const real s = n - point[m_D];
		la::mul(plane, m_world_to_local_t, m_plane, m_D + 1, m_D + 1);
		return s / la::normal_proj_plane(plane, plane, m_D);
	}

	//// Functions required by HCP_Heap_Allocating_Halfspace_Set<typename>.  See the documentation for that class. ////
	static size_t	storage_required(uint D) { return HCP_Mutable::storage_required(D) + (D + 1) * sizeof(real); }

	void*
	init(uint D, void* storage)
	{
		m_D = D;
		storage = HCP_Mutable::init(D, storage);
		m_plane = (real*)storage;
		return m_plane + (D + 1);
	}

protected:
	uint	m_D;
	real*	m_plane;
};


/**
This ellipsoid object allocates and manages its own storage.
*/
typedef HCP_Heap_Allocating_Halfspace_Set<HCP_Ellipsoid_NA> HCP_Ellipsoid;


#endif // #ifndef _HCP_ELLIPSOID_H_
