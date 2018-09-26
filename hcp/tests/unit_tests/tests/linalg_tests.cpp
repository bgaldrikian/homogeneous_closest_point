// Copyright (c) 2018 NVIDIA Corporation
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

#include "math/linalg_ext.h"

// Random distributions
#include "math/dist.h"

#include <vector>
#include <limits>
using namespace std;


UNIT_TEST(QRP_Orthogonalize_Full)
{
	Random<real> rnd;

	const uint D_max = 10;
	const uint trial_count = 10000;

	real A[D_max*D_max];
	real QR[D_max*D_max];
	real Q[D_max*D_max];
	real QT[D_max*D_max];
	real test[D_max*D_max];
	uint P[D_max];

	for (uint D = 1; D <= D_max; ++D)
	{
		const real eps = 2 * (D+1) * numeric_limits<real>().epsilon();
		for (uint N = 1; N < D; ++N)
		{
			for (uint T = 0; T < trial_count; ++T)
			{
				for (uint j = 0; j < N; ++j)
					rnd.ball(la::col(A, j, D), D, 2);
				la::cpy(QR, A, D*D);
				const uint rankA = lax::QRP_decompose(P, QR, D, N, 0);
				lax::QR_create_orthonormal_basis(Q, D, QR, D, N);

				la::transpose(QT, Q, D, D);

				// Test column orthonormality of Q:
				la::mul(test, QT, Q, D, D, D);
				for (uint i = 0; i < D; ++i) la::elem(test, i, i, D) -= 1;
				for (uint e = 0; e < D*D; ++e)
					EXPECT(abs(test[e]) <= eps);

				// Test row orthonormality of Q:
				la::mul(test, Q, QT, D, D, D);
				for (uint i = 0; i < D; ++i) la::elem(test, i, i, D) -= 1;
				for (uint e = 0; e < D*D; ++e)
					EXPECT(abs(test[e]) <= eps);

				// Ensure that columns [rankA, D) of Q are orthogonal to columns of A
				for (uint jA = 0; jA < N; ++jA)
				{
					const real* colA = la::col(A, P[jA], D);
					for (uint jQ = rankA; jQ < D; ++jQ)
					{
						const real* colQ = la::col(Q, jQ, D);
						EXPECT(abs(la::dot(colA, colQ, D)) <= eps);
					}
				}
			}
		}
	}
}

UNIT_TEST(QR_Create_Complement_Basis)
{
	Random<real> rnd;

	const uint D_max = 10;
	const uint trial_count = 10000;

	real A[D_max*D_max];
	real QR[D_max*D_max];

	for (uint D = 1; D <= D_max; ++D)
	{
		const real eps = 2 * (D + 1) * numeric_limits<real>().epsilon();
		for (uint N = 1; N < D; ++N)
		{
			for (uint T = 0; T < trial_count; ++T)
			{
				for (uint j = 0; j < N; ++j)
					rnd.ball(la::col(A, j, D), D, 2);
				la::cpy(QR, A, D*D);
				const uint rankA = lax::QRP_decompose(0, QR, D, N, 0);
				lax::QR_create_complement_basis(QR, D, N);

				// Test column orthonormality of Q complement:
				for (uint i = rankA; i < D; ++i)
				{
					for (uint j = rankA; j < D; ++j)
					{
						real e = la::dot(la::col(QR, i, D), la::col(QR, j, D), D);
						if (i == j) e -= 1;
						EXPECT(abs(e) <= eps);
					}
				}

				// Ensure that columns of Q complement are orthogonal to columns of A
				for (uint jA = 0; jA < N; ++jA)
				{
					const real* colA = la::col(A, jA, D);
					for (uint jQ = rankA; jQ < D; ++jQ)
					{
						const real* colQ = la::col(QR, jQ, D);
						EXPECT(abs(la::dot(colA, colQ, D)) <= eps);
					}
				}
			}
		}
	}
}
