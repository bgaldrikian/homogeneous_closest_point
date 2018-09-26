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

#ifndef _HCP2D_H_
#define _HCP2D_H_


/** D = 2 specialized functions for the D-specific implementation of the homogeneous closest point algorithm. */

#include "src/specialized/hcpd.h"


/** Specialization of hcpd_solve.  See the declaration in hcpd.h for the function description. */
template<>
inline real
hcpd_solve_feature<1>(real p[], const real q[], const real S[][2], const uint indices[], uint N)
{
	real detH = 1;

	if (N == 0)
	{
		cpy<2>(p, q);
		if (dot<2>(p, p) == 0) p[1] = 1;
	}
	else	// N == 1:
	{
		const real* S0 = S[indices[0]];
		detH = S0[0];
		if (detH != 0)
		{
			// We will only solve using q_D == 1 in this case
			p[0] = -S0[1] / detH;
			p[1] = 1;
		}
	}

	return detH;
}


/** Specialization of hcpd_search_features.  See the declaration in hcpd.h for the function description. */
template<>
inline bool
hcpd_search_features<1>(uint& columns, real p[], real& det_basis, const real S[][2], uint N, const real q[])
{
	// Initialize measure for best solution
	real p2_min = std::numeric_limits<real>().max();
	p[1] = 1;

	bool updated = false;
	bool solution_found = false;
	switch (N)
	{
	case 2:
		updated |= hcpd_test_feature<1>(p, p2_min, det_basis, columns, 0b10, N, S, q);
	case 1:
		updated |= hcpd_test_feature<1>(p, p2_min, det_basis, columns, 0b01, N, S, q);
		if (!updated && N > 1) break;
		solution_found |= updated;
	case 0:
		solution_found |= hcpd_test_feature<1>(p, p2_min, det_basis, columns, 0b00, N, S, q);
	}

	return solution_found;
}


#endif // #ifndef _HCP2D_H_
