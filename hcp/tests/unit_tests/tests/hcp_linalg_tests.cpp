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


#include "util/unit_test.h"

#include "src/general/hcp_linalg.h"

// Random distributions
#include "math/dist.h"

#include <limits>
using namespace std;


UNIT_TEST(HCP_QRP_Determinant)
{
	Random<real> rnd;

	const uint D_max = 10;
	const uint trial_count = 10000;

	real A[D_max*D_max];
	real B[D_max*D_max];
	uint P[D_max];
	uint Q[D_max];

	for (uint D = 1; D <= D_max; ++D)
	{
		const real eps = (real)0.00002;
		for (uint N = 1; N < D; ++N)
		{
			for (uint T = 0; T < trial_count; ++T)
			{
				for (uint j = 0; j < N; ++j) rnd.ball(la::col(A, j, D), D, 2);

				la::cpy(B, A, D*N);
				const real det_QR = hcpla::QRP_decompose(P, B, D, N, 0);

				hcpla::QR_create_complement_basis(B, D, N);

				la::cpy(B, A, D*N);
				const real det_LU = hcpla::LUPQ_decompose(P, Q, B, D);

				EXPECT_NEAR(det_QR, det_LU, eps*abs(det_QR));
			}
		}
	}
}
