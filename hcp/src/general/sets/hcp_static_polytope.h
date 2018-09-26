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

#ifndef _HCP_STATIC_POLYTOPE_H_
#define _HCP_STATIC_POLYTOPE_H_


#include "math/linalg.h"

#include <limits>


/**
Implementation of HCP_Halfspace_Set, an explicit representation of a finite number of half-spaces.

The shape represented is the interesection of the half-spaces.
*/
struct HCP_Static_Polytope : public virtual HCP_Halfspace_Set
{
	/** Default ctor, creates a blank object. */
	HCP_Static_Polytope() : m_D(0), m_halfspace_count(0), m_halfspaces(nullptr) {}

	/**
	Constructor which sets the number of spatial dimensions D, but keeps the default empty list of half-spaces.

	\param[in] D	the number of spatial dimensions.
	*/
	HCP_Static_Polytope(uint D) : m_D(D), m_halfspace_count(0), m_halfspaces(nullptr) {}

	/**
	Constructor which sets the number of spatial dimensions D and sets the list of half-spaces.

	\param[in] D				the number of spatial dimensions.
	\param[in] halfspace_count	the number of half-spaces.
	\param[in] halfspaces		an array of half-space vectors, stored as contiguous (D+1)-vectors which represent the plane equations of the half-space bounding planes.
	*/
	HCP_Static_Polytope(uint D, size_t halfspace_count, const real* halfspaces) : m_D(D) { set_halfspaces(halfspace_count, halfspaces); }

	/**
	Set the half-spaces represented by this object.

	\param[in] halfspace_count	the number of half-spaces.
	\param[in] halfspaces		an array of half-space vectors, stored as contiguous (D+1)-vectors which represent the plane equation of the half-space bounding planes.
	*/
	void
	set_halfspaces(size_t halfspace_count, const real* halfspaces)
	{
		m_halfspace_count = halfspace_count;
		m_halfspaces = halfspaces;
	}

	/** Implementation of HCP_Halfspace_Set interface. */
	real
	farthest_halfspace(real* plane, const real* point) const override
	{
		const real* halfspace = m_halfspaces;
		real greatest_s = -std::numeric_limits<real>().max();
		const real* h = m_halfspaces;
		const real* stop = h + (m_D + 1)*m_halfspace_count;
		for (; h < stop; h += (m_D + 1))
		{
			const real s = la::dot(point, h, m_D + 1);
			if (s > greatest_s)
			{
				greatest_s = s;
				halfspace = h;
			}
		}
		if (halfspace) la::cpy(plane, halfspace, m_D + 1);
		else la::zero(plane, m_D + 1);
		return greatest_s;
	}

protected:
	uint		m_D;
	size_t		m_halfspace_count;
	const real*	m_halfspaces;
};


#endif // #ifndef _HCP_STATIC_POLYTOPE_H_
