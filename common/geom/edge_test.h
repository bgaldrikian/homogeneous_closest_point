// Copyright (c) 2014-2018 Bryan Galdrikian
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

#ifndef _EDGE_TEST_H_
#define _EDGE_TEST_H_


#include "math/linalg_ext.h"
#include "math/discrete.h"

#include <memory>


/** Tests for mutual intersection of a set of halfspaces, by testing edge intersections. **/


/*
	Planes are (D+1) reals (normal & displacement), so the planes array is expected to contain (D+1)*N elements
*/
static int
edge_test(uint D, uint N, const real* planes)
{
	if (N <= D) return true;
	real* S = (real*)alloca((D + 1)*(D + 1) * sizeof(real));
	real* S_inv = (real*)alloca((D + 1)*(D + 1) * sizeof(real));
	uint* indices = (uint*)alloca((D + 1) * sizeof(uint));
	for (Choose<uint> c(N, D+1, indices); c; ++c)
	{
		for (uint j = 0; j < D + 1; ++j) la::cpy(la::col(S, j, D + 1), la::col(planes, indices[j], D + 1), D + 1);
		const real detS = lax::inv(S_inv, S, D + 1);
		if (detS == 0) continue;
		uint pos_width_count = 0;
		const real* s_inv_Di = la::col(S_inv, D, D + 1);
		for (uint i = 0; i < D + 1; ++i) pos_width_count += (uint)(*s_inv_Di++ < 0);
		if (!pos_width_count) return 0;
	}
	return 1;
}


#endif // #ifndef _EDGE_TEST_H_
