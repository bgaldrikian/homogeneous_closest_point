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

#ifndef _HCP_POLYTOPE_H_
#define _HCP_POLYTOPE_H_


/**
Implementation of HCP_Halfspace_Set, an explicit representation of a finite number of half-spaces.

The shape represented is the interesection of the half-spaces.

This class derives from HCP_Mutable, giving it a transformation matrix that is effectivly applied
to the halfspaces during the execution of farthest_halfspace.  Actually, this is transformation is
only approximately applied if it contains a non-uniform scale.  In that case, the returned half-space may
not be the one farthest from the given point when transformed.  However, upon iteration of the
homogeneous closest point algorithm, the planes returned by this function will lead to the correct result.
*/


#include "src/general/sets/hcp_static_polytope.h"
#include "src/general/sets/hcp_mutable.h"
#include "src/general/sets/hcp_allocating_halfspace_set.h"


/**
This class does not perform any allocation, and requires storage be available it when it performs a calculation.

(The suffix "NA" stands for "Non-Allocating.")
*/
struct HCP_Polytope_NA : public HCP_Static_Polytope, public HCP_Mutable
{
	/** Default ctor creates an blank object.   init(...) and HCP_Static_Polytope::set_shapes(...) must be called before this object is used. */
	HCP_Polytope_NA() : m_scratch_vector(nullptr) {}

	/**
	Constructor which sets the number of spatial dimensions D and sets the list of half-spaces.

	\param[in] D				the number of spatial dimensions.
	\param[in] storage			memoray available for storage, must point to storage_required(D) bytes when this object is used.
	\param[in] halfspace_count	the number of half-spaces.
	\param[in] halfspaces		an array of half-space vectors, stored as contiguous (D+1)-vectors which represent the local plane equations of the half-space bounding planes.
	*/
	HCP_Polytope_NA(uint D, void* storage, size_t halfspace_count, const real* halfspaces) : HCP_Static_Polytope(D, halfspace_count, halfspaces) { init(D, storage); }

	/** Implementation of HCP_Halfspace_Set interface. */
	real
	farthest_halfspace(real* plane, const real* point) const override
	{
		la::mul_t(m_scratch_vector, point, m_world_to_local_t, m_D + 1, m_D + 1); // The local point
		const real s = HCP_Static_Polytope::farthest_halfspace(m_scratch_vector, m_scratch_vector);	// Replace scratch vector with plane
		la::mul(plane, m_world_to_local_t, m_scratch_vector, m_D + 1, m_D + 1);	// Transform using inverse transpose
		return s / la::normal_proj_plane(plane, plane, m_D);
	}

	//// Functions required by HCP_Heap_Allocating_Halfspace_Set<typename>.  See the documentation for that class. ////
	static size_t	storage_required(uint D) { return HCP_Mutable::storage_required(D) + (D + 1) * sizeof(real); }

	void*
	init(uint D, void* storage)
	{
		m_D = D;
		storage = HCP_Mutable::init(D, storage);
		m_scratch_vector = (real*)storage;
		return m_scratch_vector + (D + 1);
	}

protected:
	real*	m_scratch_vector;
};


/**
This polytope object allocates and manages its own storage.
*/
typedef HCP_Heap_Allocating_Halfspace_Set<HCP_Polytope_NA> HCP_Polytope;


#endif // #ifndef _HCP_POLYTOPE_H_
