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

#ifndef _HCP_ALLOCATING_HALFSPACE_SET_H_
#define _HCP_ALLOCATING_HALFSPACE_SET_H_


#include "src/hcp_defs.h"


/*
Heap-allocating class for halfspace sets.

Templated class Halfspace_Set_NA must implement:

// The number of bytes required for storage (passed into the ctor or init functions).
// \param[in] D	the number of spatial dimensions.
// \return the number of bytes required for storage.
static size_t	storage_required(uint D);

// Sets the number of spatial dimension and storage used by this object.
// \param[in] D			the number of spatial dimensions.
// \param[in] storage	memory available for storage, must point to storage_required(D) bytes when this object is used.
// \return a pointer to the memory just beyond the used storage (the input address plus storage_required bytes).
void*	init(uint D, void* storage);
*/
template<typename Halfspace_Set_NA>
struct HCP_Heap_Allocating_Halfspace_Set : public Halfspace_Set_NA
{
	/** Default ctor does not allocate */
	HCP_Heap_Allocating_Halfspace_Set() : m_storage(nullptr) { Halfspace_Set_NA::init(0, nullptr); }

	/**
	Constructor which allocates enough memory for a given number of spatial dimensions.

	\param[in] D	the number of spatial dimensions.
	*/
	HCP_Heap_Allocating_Halfspace_Set(uint D) : m_storage(nullptr) { init(D); }

	/** dtor */
	virtual ~HCP_Heap_Allocating_Halfspace_Set() { delete[] m_storage; }

	/**
	Allocate enough memory for a given number of spatial dimensions.
	Uses the required functions storage_required and init for the templated class.
	*/
	void
	init(uint D)
	{
		delete[] m_storage;
		m_storage = new char[Halfspace_Set_NA::storage_required(D)];
		Halfspace_Set_NA::init(D, m_storage);
	}

private:
	char*	m_storage;
};


#endif // #ifndef _HCP_ALLOCATING_HALFSPACE_SET_H_
