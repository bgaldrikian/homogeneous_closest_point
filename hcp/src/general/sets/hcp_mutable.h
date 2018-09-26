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

#ifndef _HCP_MUTABLE_H_
#define _HCP_MUTABLE_H_


#include "src/general/hcp_linalg.h"


/**
Class to hold a transformation matrix, used by half-space set classes to transform
data from a local space to a world space.

This class takes a local-to-world transformation matrix as input and internally stores
its inverse transpose.
*/
struct HCP_Mutable
{
	/** Internal storage required for this class to operate. */
	static size_t	storage_required(uint D) { return 2 * (D + 1)*(D + 1) * sizeof(real) + 2 * D * sizeof(uint); }

	/** ctor creates a blank object. */
	HCP_Mutable() : m_world_to_local_t(nullptr), m_inv_scratch(nullptr) {}

	/**
	Sets the number of spatial dimension and storage used by this object.

	\param[in] D		the number of spatial dimensions.
	\param[in] storage	memory available for storage, must point to storage_required(D) bytes when this object is used.

	\return a pointer to the memory just beyond the used storage (the input address plus storage_required bytes).
	*/
	void*
	init(uint D, void* storage)
	{
		m_world_to_local_t = (real*)storage;
		m_inv_scratch = (char*)(m_world_to_local_t + (D + 1)*(D + 1));
		storage = m_inv_scratch + (D + 1)*(D + 1) * sizeof(real) + 2 * D * sizeof(uint);
		la::set(m_world_to_local_t, 1, D + 1, D + 1);
		return storage;
	}

	/**
	Set the local-to-world transformation matrix.

	\param[in] D				the number of spatial dimensions.
	\param[in] local_to_world	(D+1)x(D+1) homogeneous transformation matrix, stored in column-major format.
	*/
	void
	set_local_to_world_transform(uint D, const real* local_to_world)
	{
		hcpla::invert(m_world_to_local_t, local_to_world, D + 1, m_inv_scratch);
		la::transpose(m_world_to_local_t, D + 1);
	}

	/**
	Access the stored inverse transpose of the local-to-world transform.

	\return the transpose of the world-to-local transform, stored in column-major format.
	*/
	const real*	get_world_to_local_t() const { return m_world_to_local_t; }

protected:
	real*	m_world_to_local_t;
	char*	m_inv_scratch;
};


#endif // #ifndef _HCP_MUTABLE_H_
