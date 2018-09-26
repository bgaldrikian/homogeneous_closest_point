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

#ifndef _HCP3D_H_
#define _HCP3D_H_


/** D = 3 specialized functions for the D-specific implementation of the homogeneous closest point algorithm. */

#include "src/specialized/hcpd.h"


/**
LUP decomposition (LU decomposition with row permutation) function optimized for 2x2 matrices.

Finds the LUP decomposition of the transpose of M.

\param[in,out]	M	on input, the matrix to decompose.  On output, the LU decomposition of M^T.
\param[out]		P	the row permutation used in the LUP decomposition.

\return the determinant of M.
*/
inline real
LUP_decompose_transpose_2(real M[2][2], uint& P)
{
	real detM = 1;

	if ((P = (uint)(sq(M[1][0]) > sq(M[0][0]))) != 0)
	{
		detM = -detM;
		swap<2>(M[0], M[1]);
	}

	detM *= M[0][0];

	if (M[0][0] != 0) M[1][1] -= (M[1][0] /= M[0][0])*M[0][1];

	detM *= M[1][1];

	return detM;
}


/**
Linear equation solver which uses an LUP decomposition, optimized for 2x2 matrices.

Solves the equation M*x = b.

\param[in,out]	x	a 2-vector.  On input, the right-hand side b.  On output, the solution x.
\param[in]		LU	the LU decomposition of M, which can be obtained from LUP_decompose_transpose_2.
\param[in]		P	the row permutation used in the decomposition, which can be obtained from LUP_decompose_transpose_2.
*/
inline void
LUP_solve_transpose_2(real x[2], const real LU[2][2], uint P)
{
	// Perform row permutation to get Pb
	swap(x[0], x[P]);
	// Forward substitute to get (L^-1)Pb
	x[1] -= LU[1][0] * x[0];
	// Back substitute to get (U^-1)(L^-1)Pb
	x[1] /= LU[1][1];
	x[0] = (x[0] - LU[0][1] * x[1]) / LU[0][0];
}


/** Specialization of hcpd_solve.  See the declaration in hcpd.h for the function description. */
template<>
inline real
hcpd_solve_feature<2>(real p[], const real q[], const real S[][3], const uint indices[], uint N)
{
	real H[2][2];
	uint P;
	real detH = 1;

	switch (N)
	{
	case 0:
		cpy<3>(p, q);
		if (dot<3>(p, p) == 0) p[2] = 1;
		break;
	case 1:
	{
		const real* S0 = S[indices[0]];

		cpy<2>(H[0], S0);
		H[1][0] = -H[0][1];
		H[1][1] = H[0][0];

		p[1] = dot<2>(q, H[1]);

		detH = LUP_decompose_transpose_2(H, P);
		if (detH == 0) return 0;

		p[0] = -q[2] * S0[2];

		LUP_solve_transpose_2(p, H, P);
		p[2] = q[2];
		if (dot<3>(p, p) == 0)
		{
			p[0] = -S0[2];
			LUP_solve_transpose_2(p, H, P);
			p[2] = 1;
		}
		break;
	}
	case 2:
	{
		const real* S0 = S[indices[0]];
		const real* S1 = S[indices[1]];

		cpy<2>(H[0], S0);
		cpy<2>(H[1], S1);

		detH = LUP_decompose_transpose_2(H, P);
		if (detH == 0) return 0;

		// We will only solve using q_D == 1 in this case
		p[0] = -S0[2];
		p[1] = -S1[2];

		LUP_solve_transpose_2(p, H, P);
		p[2] = 1;
		break;
	}
	}

	return detH;
}


/** Specialization of hcpd_search_features.  See the declaration in hcpd.h for the function description. */
template<>
inline bool
hcpd_search_features<2>(uint& columns, real p[], real& det_basis, const real S[][3], uint N, const real q[])
{
	// Initialize measure for best solution
	real p2_min = std::numeric_limits<real>().max();
	p[2] = 1;

	bool updated = false;
	bool solution_found = false;
	switch (N)
	{
	case 3:
		updated |= hcpd_test_feature<2>(p, p2_min, det_basis, columns, 0b110, N, S, q);
		updated |= hcpd_test_feature<2>(p, p2_min, det_basis, columns, 0b101, N, S, q);
	case 2:
		updated |= hcpd_test_feature<2>(p, p2_min, det_basis, columns, 0b011, N, S, q);
		if (!updated && N > 2) break;
		solution_found |= updated;
		updated = false;
	case 1:
		switch (N)
		{
		case 3: updated |= hcpd_test_feature<2>(p, p2_min, det_basis, columns, 0b100, N, S, q);
		case 2: updated |= hcpd_test_feature<2>(p, p2_min, det_basis, columns, 0b010, N, S, q);
		case 1: updated |= hcpd_test_feature<2>(p, p2_min, det_basis, columns, 0b001, N, S, q);
		}
		solution_found |= updated;
		if (!updated && N > 1) break;
	case 0:
		solution_found |= hcpd_test_feature<2>(p, p2_min, det_basis, columns, 0b000, N, S, q);
	}

	return solution_found;
}


#endif // #ifndef _HCP3D_H_
