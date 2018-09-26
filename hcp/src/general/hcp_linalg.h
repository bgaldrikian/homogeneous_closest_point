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

#ifndef _HCP_LINALG_H_
#define _HCP_LINALG_H_

#include "math/linalg.h"
#include <algorithm>


/** Extended linear algebra functions used by the homogeneous closest point algorithm. */


namespace hcpla
{

/**
Create a Householder vector v out of a vector x (replaces x).
Used to create the reflection operator Q = [1 - 2*v*(v^T)/|v|^2].
Does not guarantee positive [0] component of Q*x.
Normalized such that v[0] = 1, BUT THE v[0] COMPONENT IS NOT WRITTEN.  INSTEAD,
the [0] component of Q*x is written in x[0].  x[1..D] will contain v[1..D].
The user should implement application of Q with this in mind, using v[0] = 1.
The value of stride may be set to the column size to perform row operations.
Returns the normalization factor beta = 2/|v|^2.
*/
inline real
house_in_place(real* x, uint D)
{
	real sigma = 0;
	real* xi = x + 1;
	for (uint i = 1; i < D; ++i, ++xi) sigma += (*xi)*(*xi);
	const real mu_sq = (*x)*(*x) + sigma;
	if (mu_sq == 0) return 2;	// |x| = 0
	const real Qx0 = -sqrt(mu_sq)*(((int)(*x > 0)<<1)-1);
	const real recip_v0 = 1/(*x -= Qx0);	// x[0] holds v0 now (using v[0] = x[0] + sgn(x[0])*|x|)
	xi = x;
	for (uint i = 1; i < D; ++i) *++xi *= recip_v0;
	const real v0_sq = (*x)*(*x);
	*x = Qx0;	// x[0] holds Q*x [0] component now
	return 2*v0_sq/(v0_sq + sigma);
}


/*
QR decomposition based upon Algorithm 5.4.1 (Householder QR With Column Pivoting) from Matrix Computations (4th Edition), Golub and Van Loan.

P must point to N ints if it is not NULL.
A must point to MxN reals, M >= N.

eps2 is a numerical cutoff factor for terminating the decomposition.  It should be nonnegative.

Calculates the decomposition of A: AP = QR, where Q is an MxM orthogonal matrix, R is an MxN upper-triangular matrix,
and P is an NxN permutation matrix.  Upon completion, the upper triangle of A is replaced by R, and the essential part
of the Householder vectors which encode the Q transformation are stored in the subdiagonal of A. The column permutations
due to pivoting are encoded in the array P.

Returns a value with the same magnitude as the product of the singular values of A.  If M = N, returns det(A).
*/
inline real
QRP_decompose(uint* P, real* A, uint M, uint N, real eps2)
{
	real* Arr = A;
	real d = 1;
	for (uint r = 0; r < N; ++r, Arr += M+1)
	{
		uint k = r;
		real* ar = la::col(A, r, M);
		real c2 = la::norm_2_sq(ar, M);
		const real* as = ar;
		for (uint s = r + 1; s < N; ++s)
		{
			const real cs_2 = la::norm_2_sq(++as, M);
			if (cs_2 > c2)
			{
				k = s;
				c2 = cs_2;
			}
		}
		if (c2 <= eps2) return 0;

		if (P) P[r] = k;
		if (k != r)
		{
			real* ak = la::col(A, k, M);
			for (uint i = 0; i < M; ++i) std::swap(*ar++, *ak++);
			d = -d; // Negate d to account for column swap
		}

		const real beta = house_in_place(Arr, M-r);
		const real Qx0 = *Arr;
		*Arr = 1;
		real* a = Arr+M;
		for (uint s = r+1; s < N; ++s, a += M) la::madd(a, Arr, -beta*la::dot(Arr, a, M-r), a, M-r);
		d *= -(*Arr = Qx0);	// Negate d each iteration, to account for reflection
	}

	return d;
}


/*
Build the explicit Q matrix for the complement space from a QR decomposition.
Builds in-place, after the columns of the QR decomposition.

On input, A = QR decomposition of an MxN (N <= M) matrix, with space for enough columns to form full MxM matrix
On output, A = MxM maxtrix, with first N columns untouched
Columns N through M-1 contain basis vectors spanning the space complementary to the space spanned by the decomposed matrix
*/
inline void
QR_create_complement_basis(real* A, uint M, uint N)
{
	la::zero(A + N*M, (M - N)*M);
	A += M*M - 1;
	for (uint j = M; --j >= N; A -= M + 1) *A = 1;
	real* Q = A + M;
	for (uint j = N; j--; A -= M + 1, --Q)
	{
		const real beta = 2 / (1 + la::norm_2_sq(A + 1, M - j - 1));
		real* Ac = Q;
		for (uint k = N; k < M; ++k, Ac += M)
		{
			const real beta_v_q = beta*(*Ac + la::dot(Ac + 1, A + 1, M - j - 1));
			*Ac -= beta_v_q;
			la::madd(Ac + 1, A + 1, -beta_v_q, Ac + 1, M - j - 1);
		}
	}
}


/*
Replaces a matrix M with its LUPQ decomposition.

Based upon "Outer Product LU with Complete Pivoting," from Matrix Computations (4th Edition), Golub and Van Loan.

M must point to a DxD matrix of reals stored in column-major order.
P and Q must each point to buffers of D-1 uint elements.

Upon return, M will have its LU decomposition encoded in its place, and P and Q will store row and column permutation data, respectively.

The return value of this function is the determinant of M.
*/
inline real
LUPQ_decompose(uint* P, uint* Q, real* M, uint D)
{
	real detM = 1;

	real* LUkk = M;
	for (uint k = 0; k < D - 1; ++k, LUkk += D + 1)
	{
		uint pivot_row = k;
		uint pivot_col = k;
		real abs_pivot_elem = 0;
		const real* m = LUkk;
		for (uint c = k; c < D; ++c, m += k)
		{
			for (uint r = k; r < D; ++r)
			{
				const real abs_elem = std::abs(*m++);
				if (abs_elem > abs_pivot_elem)
				{
					abs_pivot_elem = abs_elem;
					pivot_row = r;
					pivot_col = c;
				}
			}
		}

		P[k] = pivot_row;
		if (pivot_row != k)
		{
			detM = -detM;
			real* r = M + k;
			real* s = M + pivot_row;
			for (uint c = 0; c < D; ++c, r += D, s += D) std::swap(*r, *s);
		}

		Q[k] = pivot_col;
		if (pivot_col != k)
		{
			detM = -detM;
			real* c = M + k*D;
			real* s = M + pivot_col*D;
			for (uint r = 0; r < D; ++r) std::swap(*c++, *s++);
		}

		const real Ukk = *LUkk;
		detM *= Ukk;

		if (Ukk != 0)
		{
			real* e = LUkk + 1;
			for (uint r = k + 1; r < D; ++r, ++e)
			{
				real& Lrk = *e;
				Lrk /= Ukk;
				real* s = e + D;
				const real* t = LUkk + D;
				for (uint c = k + 1; c < D; ++c, s += D, t += D) *s -= Lrk*(*t);
			}
		}
	}

	detM *= *LUkk;

	return detM;
}


/*
Solve a set of linear equations using an LU decomposition of a square, full-rank matrix
along with row and column permutations.  These values can be obtained using LUP_decompose or
LUPQ_decompose (set LU to the matrix output in the M parameter of those functions).

x must point to D reals corresponding to the column vector of right-hand values of the equations.

Upon return, the D reals pointed to by x will contain the solution (column) vector.

Solve the D linear equations M*x = b using:

	hcpla::LUPQ_decompose(P, Q, M, D);
	la::cpy(x, b, D);
	hcpla::LUPQ_solve(x, M, P, Q, D);
*/
inline void
LUPQ_solve(real* x, const real* LU, const uint* P, const uint* Q, uint D)
{
	for (uint i = 0; i < D - 1; ++i) std::swap(x[i], x[P[i]]);	// Perform row permutation to get Pb
	const real* m = LU;
	for (uint r = 1; r < D; ++r) { m = LU + r; for (uint i = 0; i < r; ++i, m += D) x[r] -= (*m)*x[i]; }	// Forward substitute to get (L^-1)Pb
	for (uint r = D; r-- > 0; m -= D + 1)	// Back substitute to get (U^-1)(L^-1)Pb
	{
		const real Mrr = *m;
		const real* s = m + D;
		for (uint i = r + 1; i < D; ++i, s += D) x[r] -= (*s)*x[i];
		x[r] /= Mrr;
	}
	for (uint i = D - 1; i-- > 0;) std::swap(x[i], x[Q[i]]);	// Perform column permutation to get the solution (Q^T)(U^-1)(L^-1)Pb
}


/*
Invert the matrix using an LU decomposition of a square, full-rank matrix
along with row and column permutations.  These values can be obtained using LUPQ_decompose
(set LU to the matrix output in the M parameter of that function).

invM must point to D*D reals to hold the elements of the output column-major DxD matrix.

Upon return, the D*D reals pointed to by invM will contain the inverse of the decomposed matrix.
*/
inline void
LUPQ_invert(real* invM, const real* LU, const uint* P, const uint* Q, uint D)
{
	la::set(invM, 1, D, D);
	for (uint c = 0; c < D; ++c, invM += D) LUPQ_solve(invM, LU, P, Q, D);
}


/**
Invert a DxD matrix.

Uses LU decomposition.

invM and M must each point to D*D reals to hold the inverse and input matrix elements in column-major format.  invM and M may point to the same memory.
scratch must point to D*D*sizeof(real) + 2*(D-1)*sizeof(uint) bytes.
*/
inline void
invert(real* invM, const real* M, uint D, void* scratch)
{
	real* LU = (real*)scratch;
	uint* P = (uint*)(LU + D*D);
	uint* Q = P + (D - 1);
	la::cpy(LU, M, D*D);
	hcpla::LUPQ_decompose(P, Q, LU, D);
	hcpla::LUPQ_invert(invM, LU, P, Q, D);
}


};	// namespace hcpla


#endif // #ifndef _HCP_LINALG_H_
