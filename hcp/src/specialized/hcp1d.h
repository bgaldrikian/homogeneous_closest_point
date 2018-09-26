// Copyright (c) 2017-2018 NVIDIA Corporation
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

#ifndef _HCP1D_H_
#define _HCP1D_H_


/** D = 1 specialized functions for the D-specific implementation of the homogeneous closest point algorithm. */

#include "src/specialized/hcpd.h"


/** Specialization of hcpd_update.  See the declaration in hcpd.h for the function description. */
template<>
inline bool
hcpd_update<1>(real p[2], real& det_basis, real S[2][2], uint& N, const real[2])
{
	if (N > 1)
	{
		if (S[0][0] * S[1][0] < 0) return false;	// If we get here with two opposing planes, we have a void simplex
		cpy<2>(S[0], S[--N]);	// Otherwise keep the newest plane only
	}
	p[0] = -S[0][0] * S[0][1];	// The nearest point is the new plane surface
	p[1] = 1;					// We always have a point (w = 1) solution
	det_basis = 1;
	return true;
}


#endif // #ifndef _HCP1D_H_
