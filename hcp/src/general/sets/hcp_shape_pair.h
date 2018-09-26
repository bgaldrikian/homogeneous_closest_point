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

#ifndef _HCP_SHAPEPAIR_H_
#define _HCP_SHAPEPAIR_H_


/**
Implementation of HCP_Halfspace_Set, representing the union of two other half-space sets.

The shape defined is the intersection of the shapes defined by the two other half-space sets.
*/

#include "src/hcp_defs.h"
#include "src/general/sets/hcp_allocating_halfspace_set.h"
#include "math/linalg.h"

#include <limits>


/**
This class does not perform any allocation, and requires storage space be available it when it performs a calculation.

(The suffix "NA" stands for "Non-Allocating.")
*/
struct HCP_Shape_Pair_NA : public HCP_Halfspace_Set
{
	/** Default ctor creates an blank object.   init(...) and set_shapes(...) must be called before this object is used. */
	HCP_Shape_Pair_NA() : m_valid_shape_count(0), m_D(0), m_scratch_plane(nullptr) {}

	/**
	Initializing ctor.

	\param[in] D		the number of spatial dimensions.
	\param[in] storage	memory available for storage, must point to storage_required(D) bytes when this object is used.
	\param[in] shape1	a pointer to a user-defined implementation of HCP_Halfspace_Set (may be NULL).
	\param[in] shape2	a pointer to a user-defined implementation of HCP_Halfspace_Set (may be NULL).
	*/
	HCP_Shape_Pair_NA(uint D, void* storage, const HCP_Halfspace_Set* shape1, const HCP_Halfspace_Set* shape2) : m_valid_shape_count(0)
	{
		init(D, storage);
		set_shapes(shape1, shape2);
	}

	/**
	Sets the shapes in this pair.

	\param[in] shape1	a pointer to a user-defined implementation of HCP_Halfspace_Set (may be NULL).
	\param[in] shape2	a pointer to a user-defined implementation of HCP_Halfspace_Set (may be NULL).
	*/
	void
	set_shapes(const HCP_Halfspace_Set* shape1, const HCP_Halfspace_Set* shape2)
	{
		// Set shapes array (pack)
		m_valid_shape_count = 0;
		if (shape1 != 0) m_shapes[m_valid_shape_count++] = shape1;
		if (shape2 != 0) m_shapes[m_valid_shape_count++] = shape2;
	}

	/** Implementation of HCP_Halfspace_Set interface. */
	real
	farthest_halfspace(real* plane, const real* point) const override
	{
		// Choose the greatest dot product from either set
		real greatest_s = -std::numeric_limits<real>().max();
		for (uint shapeIndex = 0; shapeIndex < m_valid_shape_count; ++shapeIndex)
		{
			const real s = m_shapes[shapeIndex]->farthest_halfspace(m_scratch_plane, point);
			if (s > greatest_s)
			{
				greatest_s = s;
				la::cpy(plane, m_scratch_plane, m_D+1);
			}
		}
		return greatest_s;
	}

	//// Functions required by HCP_Heap_Allocating_Halfspace_Set<typename>.  See the documentation for that class. ////
	static size_t	storage_required(uint D) { return (D + 1) * sizeof(real); }

	void*
	init(uint D, void* storage)
	{
		m_D = D;
		m_scratch_plane = (real*)storage;
		return m_scratch_plane + (D + 1);
	}

protected:
	const HCP_Halfspace_Set*	m_shapes[2];
	uint						m_valid_shape_count;
	uint						m_D;
	real*						m_scratch_plane;
};


/**
This shape pair object allocates and manages its own storage.
*/
struct HCP_Shape_Pair : HCP_Heap_Allocating_Halfspace_Set<HCP_Shape_Pair_NA>
{
	HCP_Shape_Pair() : HCP_Heap_Allocating_Halfspace_Set<HCP_Shape_Pair_NA>() {}
	HCP_Shape_Pair(uint D, const HCP_Halfspace_Set* shape1, const HCP_Halfspace_Set* shape2) : HCP_Heap_Allocating_Halfspace_Set<HCP_Shape_Pair_NA>(D) { set_shapes(shape1, shape2); }
};


#endif // #ifndef _HCP_SHAPEPAIR_H_
