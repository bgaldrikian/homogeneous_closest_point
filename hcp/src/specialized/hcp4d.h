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

#ifndef _HCP4D_H_
#define _HCP4D_H_


/** D = 4 specialized functions for the D-specific implementation of the homogeneous closest point algorithm. */

#include "src/specialized/hcpd.h"


/** Returns index of the extremal element in a three-element set {e0, e1, e2} based upon comparisons c_ij. The extremal index m is such that c_mn is true, or e_m == e_n, for all n. */
inline uint	ext_index(uint c_10, uint c_21, uint c_20) { return c_10 << c_21 | (c_21&c_20) << 1; }

/** Returns index (0, 1, or 2) of the minimum argument. */
inline uint	index_of_min(real x0, real x1, real x2) { return ext_index((uint)(x1 < x0), (uint)(x2 < x1), (uint)(x2 < x0)); }

/** Returns index (0, 1, or 2) of the maximum argument. */
inline uint	index_of_max(real x0, real x1, real x2) { return ext_index((uint)(x1 > x0), (uint)(x2 > x1), (uint)(x2 > x0)); }


/**
LUP decomposition (LU decomposition with row permutation) function optimized for 3x3 matrices.

Finds the LUP decomposition of the transpose of M.

\param[in,out]	M	on input, the matrix to decompose.  On output, the LU decomposition of M^T.
\param[out]		P	the row permutation used in the LUP decomposition.

\return the determinant of M.
*/
inline real
LUP_decompose_transpose_3(real M[3][4], uint P[2])	// Using 4-columns for better data alignment.  M[j][3] will be ignored.
{
	real detM = 1;

	if ((P[0] = index_of_max(sq(M[0][0]), sq(M[1][0]), sq(M[2][0]))) != 0)
	{
		detM = -detM;
		swap<4>(M[0], M[P[0]]);
	}

	detM *= M[0][0];

	if (M[0][0] != 0)
	{
		const real Lrk1 = (M[1][0] /= M[0][0]);
		const real Lrk2 = (M[2][0] /= M[0][0]);
		M[1][1] -= Lrk1*M[0][1];
		M[1][2] -= Lrk1*M[0][2];
		M[2][1] -= Lrk2*M[0][1];
		M[2][2] -= Lrk2*M[0][2];
	}

	if ((P[1] = 1 + (uint)(sq(M[2][1]) > sq(M[1][1]))) != 1)
	{
		detM = -detM;
		swap<4>(M[1], M[2]);
	}

	detM *= M[1][1];

	if (M[1][1] != 0) M[2][2] -= (M[2][1] /= M[1][1])*M[1][2];

	detM *= M[2][2];

	return detM;
}


/**
Linear equation solver which uses an LUP decomposition, optimized for 3x3 matrices.

Solves the equation M*x = b.

\param[in,out]	x	a 3-vector.  On input, the right-hand side b.  On output, the solution x.
\param[in]		LU	the LU decomposition of M, which can be obtained from LUP_decompose_transpose_3.
\param[in]		P	the row permutation used in the decomposition, which can be obtained from LUP_decompose_transpose_3.
*/
inline void
LUP_solve_transpose_3(real x[3], const real LU[3][4], const uint P[2])
{
	// Perform row permutation to get Pb
	swap(x[0], x[P[0]]);
	swap(x[1], x[P[1]]);
	// Forward substitute to get (L^-1)Pb
	x[1] -= LU[1][0] * x[0];
	x[2] -= LU[2][0] * x[0] + LU[2][1] * x[1];
	// Back substitute to get (U^-1)(L^-1)Pb
	x[2] /= LU[2][2];
	x[1] = (x[1] - LU[1][2] * x[2]) / LU[1][1];
	x[0] = (x[0] - LU[0][1] * x[1] - LU[0][2] * x[2]) / LU[0][0];
}


/** Specialization of hcpd_solve.  See the declaration in hcpd.h for the function description. */
template<>
inline real
hcpd_solve_feature<3>(real p[], const real q[], const real S[][4], const uint indices[], uint N)
{
	real H[3][4];
	uint P[2];
	real detH = 1;

	switch (N)
	{
	case 0:
		cpy<4>(p, q);
		if (dot<4>(p, p) == 0) p[3] = 1;
		break;
	case 1:
	{
		const real* S0 = S[indices[0]];
		const real recip_n2 = 1 / dot<3>(S0, S0);
		cpy<4>(p, q);
		madd<3>(p, S0, -dot<4>(p, S0)*recip_n2, p);
		madd<3>(p, S0, -dot<4>(p, S0)*recip_n2, p);	// Need an accurate way of finding perpendiculars, preferably w/o sqrt... this iteration suffices
		if (dot<4>(p, p) == 0)
		{
			const real c = -S0[3] * recip_n2;
			p[0] = c*S0[0];
			p[1] = c*S0[1];
			p[2] = c*S0[2];
			p[3] = 1;
		}
		break;
	}
	case 2:
	{
		const real* S0 = S[indices[0]];
		const real* S1 = S[indices[1]];

		const real recip_n0_sq = 1 / dot<3>(S0, S0);

		cpy<4>(H[0], S0);
		madd<3>(H[1], S0, -dot<3>(S0, S1)*recip_n0_sq, S1);
		const real n1_sq = dot<3>(H[1], H[1]);
		if (n1_sq == 0) return 0;

		const real n1_sq_over_n0_sq = n1_sq*recip_n0_sq;

		const int i = index_of_min(sq(H[0][0]) * n1_sq_over_n0_sq + sq(H[1][0]), sq(H[0][1]) * n1_sq_over_n0_sq + sq(H[1][1]), sq(H[0][2]) * n1_sq_over_n0_sq + sq(H[1][2]));
		zero<4>(H[2]);
		H[2][i] = n1_sq;
		for (int j = 0; j < 3; ++j) H[2][j] -= H[0][j] * H[0][i] * n1_sq_over_n0_sq + H[1][j] * H[1][i];
		p[2] = dot<3>(q, H[2]);

		cpy<4>(H[1], S1);

		detH = LUP_decompose_transpose_3(H, P);
		if (detH == 0) return 0;

		p[0] = -q[3] * S0[3];
		p[1] = -q[3] * S1[3];

		LUP_solve_transpose_3(p, H, P);
		p[3] = q[3];
		if (dot<4>(p, p) == 0)
		{
			p[0] = -S0[3];
			p[1] = -S1[3];
			LUP_solve_transpose_3(p, H, P);
			p[3] = 1;
		}
		break;
	}
	case 3:
	{
		const real* S0 = S[indices[0]];
		const real* S1 = S[indices[1]];
		const real* S2 = S[indices[2]];

		cpy<4>(H[0], S0);
		cpy<4>(H[1], S1);
		cpy<4>(H[2], S2);

		detH = LUP_decompose_transpose_3(H, P);
		if (detH == 0) return 0;

		// We will only solve using q_D == 1 in this case
		p[0] = -S0[3];
		p[1] = -S1[3];
		p[2] = -S2[3];

		LUP_solve_transpose_3(p, H, P);
		p[3] = 1;
		break;
	}
	}

	return detH;
}


/** Specialization of hcpd_search_solutions.  See the declaration in hcpd.h for the function description. */
template<>
inline bool
hcpd_search_features<3>(uint& columns, real p[], real& det_basis, const real S[][4], uint N, const real q[])
{
	// Initialize measure for best solution
	real p2_min = std::numeric_limits<real>().max();
	p[3] = 1;

	bool updated = false;
	bool solution_found = false;
	switch (N)
	{
	case 4:
		updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b1110, N, S, q);
		updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b1101, N, S, q);
		updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b1011, N, S, q);
	case 3:
		updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b0111, N, S, q);
		if (!updated && N > 3) break;
		solution_found |= updated;
		updated = false;
	case 2:
		switch (N)
		{
		case 4:
			updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b1100, N, S, q);
			updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b1010, N, S, q);
			updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b1001, N, S, q);
		case 3:
			updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b0110, N, S, q);
			updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b0101, N, S, q);
		case 2:
			updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b0011, N, S, q);
		}
		if (!updated && N > 2) break;
		solution_found |= updated;
		updated = false;
	case 1:
		switch (N)
		{
		case 4: updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b1000, N, S, q);
		case 3: updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b0100, N, S, q);
		case 2: updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b0010, N, S, q);
		case 1: updated |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b0001, N, S, q);
		}
		solution_found |= updated;
		if (!updated && N > 1) break;
	case 0:
		solution_found |= hcpd_test_feature<3>(p, p2_min, det_basis, columns, 0b0000, N, S, q);
	}

	return solution_found;
}


#endif // #ifndef _HCP4D_H_
