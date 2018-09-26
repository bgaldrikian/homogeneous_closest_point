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

#ifndef _LINALG_EXT_H_
#define _LINALG_EXT_H_


#include "math/linalg.h"
#include "types.h"

#ifndef _MSC_VER
#include <alloca.h>
#else
#include <malloc.h>
#endif

#include <algorithm>
#include <limits>


/** Extended linear algebra functions: bidiagonalization, matrix cofactor calculation, orthogonalization, QR decomposition, LU decomposition, matrix inverse, linear equation solving, etc. */


namespace lax
{

// Create a Householder vector v out of a vector x (replaces x).
// Used to create the reflection operator Q = [1 - 2*v*(v^T)/|v|^2].
// Does not guarantee positive [0] component of Q*x.
// Normalized such that v[0] = 1, BUT THE v[0] COMPONENT IS NOT WRITTEN.  INSTEAD,
// the [0] component of Q*x is written in x[0].  x[1..D] will contain v[1..D].
// The user should implement application of Q with this in mind, using v[0] = 1.
// The value of stride may be set to the column size to perform row operations.
// Returns the normalization factor beta = 2/|v|^2.
inline real
house_in_place(real* x, uint D, size_t real_stride = 1)
{
	real sigma = 0;
	real* xi = x+real_stride;
	for (uint i = 1; i < D; ++i, xi += real_stride) sigma += (*xi)*(*xi);
	const real mu_sq = (*x)*(*x) + sigma;
	if (mu_sq == 0) return 2;	// |x| = 0
	const real Qx0 = -sqrt(mu_sq)*(((int)(*x > 0)<<1)-1);
	const real recip_v0 = 1/(*x -= Qx0);	// x[0] holds v0 now (using v[0] = x[0] + sgn(x[0])*|x|)
	xi = x;
	for (uint i = 1; i < D; ++i) *(xi += real_stride) *= recip_v0;
	const real v0_sq = (*x)*(*x);
	*x = Qx0;	// x[0] holds Q*x [0] component now
	return 2*v0_sq/(v0_sq + sigma);
}


// Householder vector, guarantees positive [0] component of Px (when |x| != 0).
// However, v cannot be stored as an "essential part" of D-1 values v[1...(D-1)] when x[1...(D-1)] = 0 and x[0] >= 0.
// This is because the vector v is zero in that case, so cannot be normalized.
inline real
house_pos(real* v, const real* x, uint D)
{
	real s = 0;
	for (uint i = 1; i < D; ++i)
	{
		v[i] = x[i];
		s += v[i]*v[i];
	}

	if (s == 0 && x[0] >= 0)
	{
		v[0] = 0;
		return 0;
	}

	const real mu = sqrt(x[0]*x[0] + s);
	v[0] = x[0] > 0 ? -s/(x[0] + mu) : x[0] - mu;
	return 2/(v[0]*v[0] + s);
}


/*
	Bidiagonalization based upon Algorithm 5.4.2 (Householder Bidiagonalization) from Matrix Computations (4th Edition), Golub and Van Loan
*/

// Bidiagonalization of a DxD matrix M, in place.  Encodes the Householder reflections in the unused elements of the result.  Returns det(M).
inline real
bidiag(real* M, uint D)
{
	if (D == 1) return *M;	// Special case, since this would have an odd number (1) of reflections and give the wrong determinant sign
	real detM = 1;
	real* Mii = M;
	real* m;
	for (uint i = 0; i < D; ++i, Mii += D+1)	// We could skip i = D-1, but this way we get an even number of reflections so that determinants preserve their signs
	{
		real b = house_in_place(Mii, D-i);
		real Qx0 = *Mii;
		*Mii = 1;
		m = Mii+D;
		const uint i1 = i+1;
		for (uint j = i1; j < D; ++j, m += D) la::madd(m, Mii, -b*la::dot(Mii, m, D-i), m, D-i);
		*Mii = Qx0;
		if (i < D-2)
		{
			real* M_i1 = Mii+D;
			b = house_in_place(M_i1, D-i1, D);
			Qx0 = *M_i1;
			*M_i1 = 1;
			for (uint j = 1; j < D-i; ++j)
			{
				real bva = 0;
				m = M_i1;
				for (uint k = i1; k < D; ++k, m += D) bva += (*m)*m[j];
				bva *= b;
				m = M_i1;
				for (uint k = i1; k < D; ++k, m += D) m[j] -= bva*(*m);
			}
			*M_i1 = Qx0;
		}
		detM *= *Mii;
	}
	return detM;
}


// Apply Householder reflections encoded in bidiagonalized matrix A (from bidiag) to M in reverse order, in place
inline void
inv_bidiag_transform(real* M, const real* A, uint D)
{
	if (D == 1) return;
	const real* Aii = A + D*D-1;
	const real* a;
	real* M0i = M + (D-1)*D;
	for (uint i = D; i--; Aii -= D+1, M0i -= D)
	{
		const uint i1 = i+1;
		if (i < D-2)
		{
			const real* Ai1 = Aii+D;
			real b = 1;
			a = Ai1+D;
			for (uint j = i1+1; j < D; ++j, a += D) b += (*a)*(*a);
			b = 2/b;
			for (uint j = 0; j < D; ++j)
			{
				a = Ai1;
				real* m = M0i + D;
				real bva = m[j];	// Implicit v[0] = 1
				for (uint k = i1+1; k < D; ++k) bva += *(a += D)*(m += D)[j];
				bva *= b;
				a = Ai1;
				m = M0i + D;
				m[j] -= bva;	// Implicit v[0] = 1
				for (uint k = i1+1; k < D; ++k) (m += D)[j] -= *(a += D)*bva;
			}
		}
		real b = 1;
		a = Aii+1;
		for (uint j = i1; j < D; ++j, ++a) b += (*a)*(*a);
		b = 2/b;
		real * m = M;
		a = Aii+1;
		for (uint j = 0; j < D; ++j, m += D)
		{
			const real bva = b*(m[i] + la::dot(a, m+i1, D-i1));	// Implicit v[0] = 1
			m[i] -= bva;	// Implicit v[0] = 1
			la::madd(m+i1, a, -bva, m+i1, D-i1);	
		}
	}
}


/*
	Calculate the cofactors of an upper bidiagonal matrix

	B and cofB must point to DxD real (column-major) matrices
	
	B is the input bidiagonal matrix.  Only the diagonal and first superdiagonal are read.
*/
inline void
cof_bidiag(real* cofB, const real* B, uint D)
{
	la::zero(cofB, D*D);
	real* c = cofB;
	for (uint j = 0; j < D; c += ++j)
	{
		for (uint i = j; i < D; ++i)	// cofactor matrix of an upper-bidagonal matrix is lower-triangular		
		{
			real minor = 1;
			const real* b = B;
			uint m;
			for (m = 0; m < j; ++m, b += D+1) minor *= *b;
			for (b += D; m < i; ++m, b += D+1) minor *= *b;
			for (++b, ++m; m < D; ++m, b += D+1) minor *= *b;
			*c++ = (1-(int)(((i^j)&1)<<1))*minor;
		}
	}
}


/*
	Cofactor function which uses bidiagonalization.
	
	M and cofM must point to DxD real (column-major) matrices
	scratch must point to D*D*sizeof(real)

	Returns det(M).
*/
inline real
cof_w(real* cofM, const real* M, uint D, void* scratch)
{
	// Bidiagonalize M
	real* B = (real*)scratch;
	la::cpy(B, M, D*D);
	const real detM = bidiag(B, D);

	// Calculate cof(M) by first finding the cofactor matrix of B, an upper-bidiagonal matrix
	cof_bidiag(cofM, B, D);

	// Inverse transform to obtain cof(M)
	inv_bidiag_transform(cofM, B, D);

	return detM;
}

// This version creates its own scratch space
inline real
cof(real* cofM, const real* M, uint D)
{
	return cof_w(cofM, M, D, alloca(D*D*sizeof(real)));
}


/*
	Determinant of a square DxD matrix, using the cof function.
*/
inline real
det(const real* M, uint D)
{
	return cof((real*)alloca(D*D*sizeof(real)), M, D);
}


/*
	perp - acccurate perpendicular to D-1 D-dimensional vectors.  D must be at least 2.

	W must point to D*D reals.  The first (D-1)*D reals are the input D-1 D-vectors.  The last D reals is the output D-vector.
	For convenience, returns a pointer to the output vector.

	scratch must point to 2*D*D*sizeof(real)
*/
inline real*
perp_w(real* W, uint D, void* scratch)
{
	real* p = la::col(W, D-1, D);	// alias for last vector in W
	la::zero(p, D);
	real* C = (real*)scratch;
	scratch = C + D*D;
	cof_w(C, W, D, scratch);
	const real sign = (real)(((int)(D&1)<<1)-1);
	la::mul(p, la::col(C, D-1, D), sign, D);
	const real W2 = la::norm_2_sq(p, D);
	if (W2 != 0)	// Iterative improvement to (W[0] ... W[D-2])(p) = (0)
	{
		cof_w(C, W, D, scratch);
		real* q = (real*)scratch;
		la::zero(q, D);
		real* w = W;
		real* c = C;
		for (uint i = 0; i < D-1; ++i, w += D, c += D) la::madd(q, c, sign*la::dot(w, p, D), q, D);
		la::madd(p, q, -1/W2, p, D);
	}
	return p;
}

// This version creates its own scratch space
inline real*
perp(real* W, uint D)
{
	return perp_w(W, D, alloca(2 * D*D * sizeof(real)));
}


/*
	Cofactor matrix, high-precision.  Uses perp to calculate each column.

	M and cofM must point to DxD real (column-major) matrices
	scratch must point to 3*D*D*sizeof(real) bytes
*/
inline void
cofx_w(real* cofM, const real* M, uint D, void* scratch)
{
	real* W = (real*)scratch;
	scratch = W + D*D;
	real* p = la::col(W, D-1, D);	// alias for last vector in W
	la::cpy(W, la::col(M, 1, D), D*(D-1));
	perp_w(W, D, scratch);
	la::cpy(cofM, p, D);
	real* w = W;
	for (uint minor_col_to_replace = 0; minor_col_to_replace < D-1; ++minor_col_to_replace, w += D, M += D)
	{
		la::cpy(w, M, D);
		perp_w(W, D, scratch);
		cofM += D;
		la::mul(cofM, p, (real)(((int)(minor_col_to_replace&1)<<1)-1), D);
	}
}

// This version creates its own scratch space
inline void
cofx(real* cofM, const real* M, uint D)
{
	cofx_w(cofM, M, D, alloca(3 * D * D * sizeof(real)));
}


/*
	Orthogonalize MxN matrix Q to first order.

	scratch must point to (M+N)*N*sizeof(real) bytes
*/
inline void
orthogonalize(real* Q, uint M, uint N, void* scratch)
{
	real* minusHalfQTQ = (real*)scratch;
	real* minusHalfQQTQ = minusHalfQTQ + N*N;
	real* e = minusHalfQTQ;
	real* Qj = Q;
	for (uint j = 0; j < N; ++j, Qj += M)	// Create -(1/2)*Q^T*Q
	{
		e += j;
		*e++ = -(real)0.5*la::dot(Qj, Qj, M);
		real* t = e+N-1;
		real* Qi = Qj+M;
		for (uint i = j+1; i < N; ++i, ++e, t += N, Qi += M) *e = *t = -(real)0.5*la::dot(Qi, Qj, M);
	}
	la::mul(minusHalfQQTQ, Q, minusHalfQTQ, M, N, N);	// Create -(1/2)*Q*(Q^T*Q)
	la::madd(Q, Q, (real)1.5, minusHalfQQTQ, M*N);	// Q = (3/2)*Q - (1/2)*Q*(Q^T*Q)
}


/*
	Newton iteration to improve NxN matrix inverse.  Requires original matrix M and inverse invM to be improved.

	scratch must point to 2*N*N*sizeof(real) bytes
*/
inline void
improve_inverse(real* invM, const real* M, uint N, void* scratch)
{
	real* invMxM = (real*)scratch;
	real* invMxMxinvM = invMxM + N*N;
	la::mul(invMxM, invM, M, N, N, N);
	la::mul(invMxMxinvM, invMxM, invM, N, N, N);
	la::msub(invM, invM, 2, invMxMxinvM, N*N);
}


/*
	R component of QR decomposition based upon Algorithm 5.2.1 (Householder QR) from Matrix Computations (4th Edition), Golub and Van Loan

	A must point to MxN reals, M >= N
	scratch must point to M*sizeof(real) bytes

	Calculates the R component of the decomposition A = QR, where Q is an MxM orthogonal matrix and R is an MxN upper-triangular matrix.
	The diagonal elements R(i,i) will be positive semidefinite (definite if A has full rank).  Upon completion, A is replaced by R.
*/
inline void
upper_triangularize_w(real* A, uint M, uint N, void* scratch)
{
	real* v = (real*)scratch;
	for (uint j = 0; j < N; ++j, A += M+1)
	{
		const real beta = house_pos(v+j, A, M-j);
		real* a = A;
		for (uint k = j; k < N; ++k, a += M) la::madd(a, v+j, -beta*la::dot(v+j, a, M-j), a, M-j);
	}
}


/*
	QR decomposition based upon Algorithm 5.1.1 (Householder QR) from Matrix Computations (4th Edition), Golub and Van Loan

	A must point to MxN reals, M >= N

	Calculates decomposition of A: A = QR, where Q is an MxM orthogonal matrix and R is an MxN upper-triangular matrix.
	Upon completion, the upper triangle of A is replaced by R, and the essential part of the Householder vectors which
	encode the Q transformation are stored in the subdiagonal of A.
*/
inline void
QR_decompose(real* A, uint M, uint N)
{
	for (uint j = 0; j < N; ++j, A += M+1)
	{
		const real beta = house_in_place(A, M-j);
		const real Qx0 = *A;
		*A = 1;
		real* a = A+M;
		for (uint k = j+1; k < N; ++k, a += M) la::madd(a, A, -beta*la::dot(A, a, M-j), a, M-j);
		*A = Qx0;
	}
}


/*
	QR decomposition based upon Algorithm 5.4.1 (Householder QR With Column Pivoting) from Matrix Computations (4th Edition), Golub and Van Loan

	P must point to N ints if it is not NULL.
	A must point to MxN reals, M >= N.

	eps2 is a numerical cutoff factor for terminating the decomposition.  It should be nonnegative.

	Calculates the decomposition of A: AP = QR, where Q is an MxM orthogonal matrix, R is an MxN upper-triangular matrix,
	and P is an NxN permutation matrix.  Upon completion, the upper triangle of A is replaced by R, and the essential part
	of the Householder vectors which encode the Q transformation are stored in the subdiagonal of A. The column permutations
	due to pivoting are encoded in the array P.

	Returns the rank of A, equal to the number of columns in the decomposition.
*/
inline uint
QRP_decompose(uint* P, real* A, uint M, uint N, real eps2)
{
	real* Arr = A;
	uint r = 0;
	for (; r < N; ++r, Arr += M+1)
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
		if (c2 <= eps2) break;

		if (P) P[r] = k;
		real* ak = la::col(A, k, M);
		for (uint i = 0; i < M; ++i) std::swap(*ar++, *ak++);

		const real beta = house_in_place(Arr, M-r);
		const real Qx0 = *Arr;
		*Arr = 1;
		real* a = Arr+M;
		for (uint s = r+1; s < N; ++s, a += M) la::madd(a, Arr, -beta*la::dot(Arr, a, M-r), a, M-r);
		*Arr = Qx0;
	}

	return r;
}


/*
	Solves (A^T)*x = b for x, where A is an MxN matrix, M >= N, and x and b are N-vectors.
	QR and P are the decomposition of A given by QRP_decompose.
	x must point to M reals, and should have the first N elements initialized to b.
	x will hold the solution upon return.
*/
inline void
QRP_solve_transpose(real* x, const real* QR, uint M, uint N, const uint* P)
{
	// Pad x with 0
	la::zero(x+N, M-N);

	// Apply P^T
	for (uint i = 0; i < N; ++i) std::swap(x[i], x[P[i]]);

	// Solve (R^T)*x = b for x (using x for b), by row-wise forward substitution
	x[0] /= *QR++;
	for (uint i = 1; i < N; ++i)
	{
		QR += M-i;
		for (uint j = 0; j < i; ++j) x[i] -= (*QR++)*x[j];
		x[i] /= *QR++;
	}

	// Now calculate Q*x (replacing x)
	for (uint j = N; j--; QR -= M+1)
	{
		const real beta_v_y = 2*(x[j] + la::dot(QR, x+j+1, M-j-1))/(1 + la::norm_2_sq(QR, M-j-1));
		x[j] -= beta_v_y;
		la::madd(x+j+1, QR, -beta_v_y, x+j+1, M-j-1);
	}
}


/*
	b = P*(QR^T)*x

	Using the QR decomposition from QR_decompose or QRP^T decomposition from QRP_decompose, apply (QR)^T or P*(QR)^T to the vector x.
	QR must point to an MxN matrix, M >= N, and P must either point to an array of N ints or be NULL (representing the identity).
	x must point to an M-vector, and the result b is written to the first N elements of x.
*/
inline void
QRP_apply_transpose(real* x, const real* QR, uint M, uint N, const uint* P)
{
	// Calculate (Q^T)*x
	QR -= M;
	for (uint j = 0; j < N; ++j)
	{
		QR += M+1;
		const real beta_v_x = 2*(x[j] + la::dot(QR, x+j+1, M-j-1))/(1 + la::norm_2_sq(QR, M-j-1));
		x[j] -= beta_v_x;
		la::madd(x+j+1, QR, -beta_v_x, x+j+1, M-j-1);
	}

	// Apply R^T to obtain (QR)^T*x
	for (uint i = N; i--; QR -= M-i)
	{
		x[i] *= *(--QR);
		for (uint j = i; j--;) x[i] += *(--QR)*x[j];
	}

	// Apply P to obtain P*(QR)^T*x
	for (uint i = N; i--;) std::swap(x[i], x[P[i]]);
}


/*
Solves (R^T)*R*x = b for x, where R is an upper triangluar MxN matrix, M >= N, and x and b are N-vectors.

x must be initialized to the value of b.  Upon return it will be the solution to the above equation.
*/
inline void
solve_symmetrized(real* x, const real* R, uint M, uint N)
{
	// Solve (R^T)*y = b for y (using x for y), by row-wise forward substitution
	x[0] /= *R;
	for (uint i = 1; i < N; ++i)
	{
		R += M + 1 - i;
		for (uint j = 0; j < i; ++j) x[i] -= (*R++)*x[j];
		x[i] /= *R;
	}

	// Solve R*x = y for x, by column-wise backsubstitution
	for (uint j = N; --j;)
	{
		x[j] /= *R;
		for (uint i = j; i--;) x[i] -= x[j] * (*--R);
		R -= M + 1 - j;
	}
	x[0] /= *R;
}


/*
Replaces a matrix M with its LUPQ decomposition.

Based upon "Outer Product LU with Complete Pivoting," from Matrix Computations (4th Edition), Golub and Van Loan.

M must point to a DxD matrix of reals stored in column-major order.
P and Q must each point to buffers of D-1 uint elements.

Upon return, M will have its LU decomposition encoded in its place, and P and Q will store row and column permutation data, respectively.

The return value of this function is the determinant of M.
*/
enum PivotingStrategy { CompletePivoting, RookPivoting, PartialPivoting };

template<int Pivoting = CompletePivoting>
inline real
LUPQ_decompose(uint* P, uint* Q, real* M, uint D)	// If Pivoting == PartialPivoting, Q is ignored
{
	real detM = 1;

	real* LUkk = M;
	for (uint k = 0; k < D - 1; ++k, LUkk += D + 1)
	{
		uint pivot_row;
		uint pivot_col;
		real abs_pivot_elem = 0;

		switch (Pivoting)
		{
		default:
		case CompletePivoting:
		{
			pivot_row = k;
			pivot_col = k;
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
			break;
		}
		case RookPivoting:
		{
			pivot_row = k;
			pivot_col = k;
			bool scan = true;
			for (bool scan_rows = true; scan; scan_rows = !scan_rows)
			{
				scan = false;
				if (scan_rows)
				{
					const real* m = LUkk + (pivot_col - k)*D;
					for (uint r = k; r < D; ++r)
					{
						const real abs_elem = std::abs(*m++);
						if (abs_elem > abs_pivot_elem)
						{
							abs_pivot_elem = abs_elem;
							pivot_row = r;
							scan = true;
						}
					}
				}
				else
				{
					const real* m = LUkk + (pivot_row - k);
					for (uint c = k; c < D; ++c, m += D)
					{
						const real abs_elem = std::abs(*m);
						if (abs_elem > abs_pivot_elem)
						{
							abs_pivot_elem = abs_elem;
							pivot_col = c;
							scan = true;
						}
					}
				}
			}
			break;
		}
		case PartialPivoting:
		{
			(void)pivot_col;
			pivot_row = k;
			const real* m = LUkk;
			for (uint r = k; r < D; ++r)
			{
				const real abs_elem = std::abs(*m++);
				if (abs_elem > abs_pivot_elem)
				{
					abs_pivot_elem = abs_elem;
					pivot_row = r;
				}
			}
			break;
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

		if (Pivoting != PartialPivoting)
		{
			Q[k] = pivot_col;
			if (pivot_col != k)
			{
				detM = -detM;
				real* c = M + k*D;
				real* s = M + pivot_col*D;
				for (uint r = 0; r < D; ++r) std::swap(*c++, *s++);
			}
		}
		else
		{
			(void)Q;
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
Replaces a matrix M with its LUP decomposition.

Based upon "Outer Product LU with Partial Pivoting," from Matrix Computations (4th Edition), Golub and Van Loan.

M must point to a DxD matrix of reals stored in column-major order.
P must point to buffers of D-1 uint elements.

Upon return, M will have its LU decomposition encoded in its place, and P will store row permutation data, respectively.

The return value of this function is the determinant of M.
*/
inline real
LUP_decompose(uint* P, real* M, uint D)
{
	return LUPQ_decompose<PartialPivoting>(P, 0, M, D);
}


/*
	Solve a set of linear equations using an LU decomposition of a square, full-rank matrix
	along with row and column permutations.  These values can be obtained using LUP_decompose or
	LUPQ_decompose (set LU to the matrix output in the M parameter of those functions).

	x must point to D reals corresponding to the column vector of right-hand values of the equations.

	Upon return, the D reals pointed to by x will contain the solution (column) vector.

	Solve the D linear equations M*x = b using:

		lax::LUPQ_decompose(P, Q, M, D);
		la::cpy(x, b, D);
		lax::LUPQ_solve(x, M, P, Q, D);

	Iterative improvement requires one keep a copy of the orginal matrix:

		la::cpy(LU, M, D*D);
		lax::LUPQ_decompose(P, Q, LU, D);
		la::cpy(x, b, D);
		lax::LUPQ_solve(x, LU, P, Q, D);
		la::mul(d, M, x, D, D);
		la::sub(d, d, b, D);
		lax::LUPQ_solve(d, LU, P, Q, D);
		la::sub(x, x, d, D);
*/
inline void
LU_solve(real* x, const real* LU, uint D)
{
	const real* m = LU;
	for (uint r = 1; r < D; ++r) { m = LU + r; for (uint i = 0; i < r; ++i, m += D) x[r] -= (*m)*x[i]; }	// Forward substitute to get (L^-1)Pb
	for (uint r = D; r-- > 0; m -= D + 1)	// Back substitute to get (U^-1)(L^-1)Pb
	{
		const real Mrr = *m;
		const real* s = m + D;
		for (uint i = r + 1; i < D; ++i, s += D) x[r] -= (*s)*x[i];
		x[r] /= Mrr;
	}
}


inline void
LUP_solve(real* x, const real* LU, const uint* P, uint D)
{
	for (uint i = 0; i < D - 1; ++i) std::swap(x[i], x[P[i]]);	// Perform row permutation to get Pb
	LU_solve(x, LU, D);
}


inline void
LUPQ_solve(real* x, const real* LU, const uint* P, const uint* Q, uint D)
{
	LUP_solve(x, LU, P, D);
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


/*
	Calculate the inverse of a square, full-rank matrix.

	M must point to D*D reals holding the elements of the input matrix, stored in column-major format.
	invM must point to a buffer of D*D reals to hold the output (inverse) matrix, stored in column-
	major format.

	scratch must point to D*D*sizeof(real) + 2*(D-1)*sizeof(uint) bytes.

	Upon return, if M was non-singular then the invM buffer will be filled with the elements of the
	inverse of M.  Otherwise, the invM buffer is not modified.

	Note: invM and M may point to the same buffer.

	The function returns the determinant of M.
*/
inline real
inv_w(real* invM, const real* M, uint D, void* scratch)
{
	real* LU = (real*)scratch;
	uint* Pi_r = (uint*)(LU + D*D);
	uint* Pi_c = Pi_r + (D - 1);
	la::cpy(LU, M, D * D);
	const real detM = LUPQ_decompose(Pi_r, Pi_c, LU, D);
	if (detM != 0)
	{
		la::zero(invM, D*D);
		for (uint c = 0; c < D; ++c, invM += D)
		{
			invM[c] = 1;
			LUPQ_solve(invM, LU, Pi_r, Pi_c, D);
		}
	}
	return detM;
}

// This version creates its own scratch space
inline real
inv(real* invM, const real* M, uint D)
{
	return inv_w(invM, M, D, alloca(D * D * sizeof(real) + 2 * (D - 1) * sizeof(uint)));
}


/*
Orthogonalize the columns of a matrix using Modified Graham-Schmidt.

\param[out]		P				a buffer of N uint elements, which will hold the column permutation used in the orthogonalization
\param[in,out]	A				the matrix to orthogonalize
\param[in]		column_stride	the number of real entries per column of A
\param[in]		M				the number of rows of A to operate on (must be equal to or less than column_stride)
\param[in]		N				the number of columns of A to operate on
\param[in]		eps2			a cutoff factor, column normalization will terminate when square norms fall below this

\return the number of columns in the orthogonalization, equal to the column rank of A
*/
inline uint
MGS_orthogonalize(uint* P, real* A, uint column_stride, uint M, uint N, real eps2)
{
	for (uint j = 0; j < N; ++j) P[j] = j;
	real* Aj = A;
	const uint rank_max = M < N ? M : N;
	for (uint j = 0; j < rank_max; ++j, Aj += column_stride)
	{
		real col_norm2 = la::norm_2_sq(Aj, M);
		real* Ak = Aj;
		real* A_col = Aj;
		uint col = j;
		for (uint k = j + 1; k < N; ++k)
		{
			A_col += column_stride;
			const real norm2 = la::norm_2_sq(A_col, M);
			if (norm2 > col_norm2)
			{
				Ak = A_col;
				col_norm2 = norm2;
				col = k;
			}
		}
		if (col_norm2 <= eps2) return j;
		for (uint i = 0; i < M; ++i) std::swap(Ak[i], Aj[i]);
		std::swap(P[col], P[j]);
		la::mul(Aj, Aj, 1 / sqrt(col_norm2), M);
		A_col = Aj;
		for (uint k = j + 1; k < N; ++k)
		{
			A_col += column_stride;
			la::madd(A_col, Aj, -la::dot(A_col, Aj, M), A_col, M);
		}
	}
	return rank_max;
}


/*
Complete a basis of column vectors.

The input matrix *MUST* consist of orthonormal vectors packed into the initial columns.

\param[in,out]	A				the matrix upon which to operate, which must hold M columns
\param[in]		column_stride	the number of real entries per column of A
\param[in]		M				the number of rows of A to operate on (must be equal to or less than column_stride)
\param[in]		N				the number of orthonormal column vectors initially in A, must be less than or equal to M

Note: if M == N then this function does nothing as the basis will already be complete.

Upon return, columns N to M-1 will contain normalized basis vectors that span the
subspace complementary to the subspace spanned by the vectors in columns 0 to N-1.
*/
inline void
complete_basis(real* A, uint column_stride, uint M, uint N)
{
	if (N >= M) return;	// Nothing to do if N == M, and N should never be greater than M.
	// Store the current row square norm in the last column of A.  This is needed only until we fill in the last column with a basis vector.
	real* row_norm2 = la::col(A, M - 1, column_stride);
	uint i_min = 0;
	real row_norm2_min = std::numeric_limits<real>().max();
	for (uint i = 0; i < M; ++i)
	{
		real* aij = A + i;
		real norm2 = 0;
		for (uint j = 0; j < N; ++j, aij += column_stride) norm2 += (*aij)*(*aij);
		row_norm2[i] = norm2;
		if (norm2 < row_norm2_min)
		{
			i_min = i;
			row_norm2_min = norm2;
		}
	}
	// Now build the complement subspace while updating row_norm2
	real* Aj = la::col(A, N, M);
	for (uint j = N; j < M; ++j, Aj += column_stride)
	{
		la::zero(Aj, M);
		Aj[i_min] = 1;
		real* Ak = A;
		for (uint k = 0; k < j; ++k, Ak += column_stride) la::madd(Aj, Ak, -Ak[i_min], Aj, M);
		row_norm2_min = std::numeric_limits<real>().max();
		real col_norm2 = 0;
		for (uint i = 0; i < M; ++i)
		{
			const real aij2 = Aj[i] * Aj[i];
			col_norm2 += aij2;
			if (j < M - 1 && (row_norm2[i] += aij2) < row_norm2_min)	// Do not write to row_norm2 for the last column (j = M - 1)
			{
				i_min = i;
				row_norm2_min = row_norm2[i];
			}
		}
		la::mul(Aj, Aj, 1 / sqrt(col_norm2), M);
	}
}


/*
Create a complete orthonormal basis split into two subspaces, one which spans the space spanned by
the columns of the original input matrix, and one which spans the complementary space.

\param[out]		P				a buffer of N uint elements, which will hold the column permutation used in the orthogonalization
\param[in,out]	A				the matrix upon which to operate, which must hold M columns
\param[in]		column_stride	the number of real entries per column of A
\param[in]		M				the number of rows of A to operate on (must be equal to or less than column_stride)
\param[in]		N				the number of orthonormal column vectors initially in A, must be less than or equal to M
\param[in]		eps2			a cutoff factor, column normalization will terminate when square norms fall below this

This function returns the column rank R of the original matrix A.  Upon return, the first R
columns of A will hold orthonormal vectors that span the space spanned by the columns of the
original matrix A.  The last M-R columns will hold orthonormal vectors that span the complementary
space.

Upon return P will hold R values, where R is the rank of the original matrix A.

\return the column rank of the original matrix A.
*/
inline uint
MGS_orthogonalize_full(uint* P, real* A, uint column_stride, uint M, uint N, real eps2)
{
	const uint R = MGS_orthogonalize(P, A, column_stride, M, N, eps2);
	complete_basis(A, column_stride, M, R);
	return R;
}


/*
Build the explicit Q matrix from a QR decomposition

Q = output MxL maxtrix of basis column vectors (1 <= L <= M)
QR = input decomposition of an MxN matrix
*/
inline void
QR_create_orthonormal_basis(real* Q, uint L, const real* QR, uint M, uint N)
{
	la::set(Q, 1, M, L);
	QR += M*(N - 1) + N;
	for (uint j = N; j--; QR -= M + 1)
	{
		const real beta = 2 / (1 + la::norm_2_sq(QR, M - j - 1));
		for (uint k = j; k < L; ++k)
		{
			real* Qc = &la::elem(Q, j, k, M);
			const real beta_v_q = beta*(*Qc + la::dot(Qc + 1, QR, M - j - 1));
			*Qc -= beta_v_q;
			la::madd(Qc + 1, QR, -beta_v_q, Qc + 1, M - j - 1);
		}
	}
}


/*
Build the explicit Q matrix for the complement space from a QR decomposition
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

};	// namespace lax


#endif // #ifndef _LINALG_EXT_H_
