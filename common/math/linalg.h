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

#ifndef _LINALG_H_
#define _LINALG_H_


#include "types.h"

#include <string.h>
#include <cmath>


/** Basic linear algebra functions. */


namespace la
{

/** Zeroes the elements of a vector (or matrix, if D = M*N). */
inline	void		zero(real* r, uint D)												{ memset(r, 0, D*sizeof(real)); }

/** Sets the elements of a vector to the given value. */
inline	void		set(real* r, real c, uint D)										{ for (uint i = 0; i < D; ++i) r[i] = c; }

/** Copies a vector (or matrix, if D = M*N). */
inline	void		cpy(real* r, const real* v, uint D)									{ memmove(r, v, D*sizeof(real)); }

/** Negates a vector (or matrix, if D = M*N). */
inline	void		neg(real* r, const real* v, uint D)									{ for (uint i = 0; i < D; ++i) r[i] = -v[i]; }

/** Multiplies a vector (or matrix, if D = M*N) by a scalar. */
inline	void		mul(real* r, const real* v, real c, uint D)							{ for (uint i = 0; i < D; ++i) r[i] = c*v[i]; }

/** Adds two vectors (or matrices, if D = M*N). */
inline	void		add(real* r, const real* u, const real* v, uint D)					{ for (uint i = 0; i < D; ++i) r[i] = u[i] + v[i]; }

/** Sub - subtracts two vectors (or matrices, if D = M*N). */
inline	void		sub(real* r, const real* u, const real* v, uint D)					{ for (uint i = 0; i < D; ++i) r[i] = u[i] - v[i]; }

/** Multiply-add.  r = u*c + v, for D-vectors (or MxN matrices, if M*N = D) r, u, and v. */
inline	void		madd(real* r, const real* u, real c, const real* v, uint D)			{ for (uint i = 0; i < D; ++i) r[i] = c*u[i] + v[i]; }

/** Multiply-subtract.  r = u*c - v, for D-vectors (or MxN matrices, if M*N = D) r, u, and v. */
inline	void		msub(real* r, const real* u, real c, const real* v, uint D)			{ for (uint i = 0; i < D; ++i) r[i] = c*u[i] - v[i]; }

/** \return the inner product of two D-vectors. */
inline	real		dot(const real* u, const real* v, uint D)							{ real r = 0; for (uint i = 0; i < D; ++i) r += u[i]*v[i]; return r; }

/** \return the square of the vector 2-norm. */
inline	real		norm_2_sq(const real* u, uint D)									{ return dot(u, u, D); }

/** \return the vector 2-norm. */
inline	real		norm_2(const real* u, uint D)										{ return sqrt(norm_2_sq(u, D)); }

/** \return the square of the 2-norm distance between two vectors. */
inline	real		diff_norm_2_sq(const real* u, const real* v, uint D)				{ real r = 0; for (uint i = 0; i < D; ++i) { const real d = u[i]-v[i]; r += d*d; } return r; }

/** \returns the square of the vector 2-norm of the sum of two vectors. */
inline	real		sum_norm_2_sq(const real* u, const real* v, uint D)					{ real r = 0; for (uint i = 0; i < D; ++i) { const real s = u[i]+v[i]; r += s*s; } return r; }

/**
Creates the normalized vector.  r and v may point to the same vector.
\return the input vector's length.
*/
inline	real
normal(real* r, const real* v, uint D)
{
	const real l2 = norm_2_sq(v, D);
	if (l2 == 0) return 0;
	const real recipL = 1/sqrt(l2);
	mul(r, v, recipL, D);
	return recipL*l2;
}

/**
Creates the normalized homogeneous point ([D]-element = 1).  Input and output vectors are (D+1)-dimensional.  r and v may point to the same vector.
\return the input vector's [D]-element
*/
inline	real
normal_proj_point(real* r, const real* v, uint D)
{
	const real w = v[D];
	if (w == 0) return 0;
	mul(r, v, 1/w, D);
	r[D] = 1;
	return w;
}

/**
Creates the normalized homogeneous plane (|D-vector| = 1).  Input and output vectors are (D+1)-dimensional.  r and v may point to the same vector.
\return the vector 2-norm of input vector's first D elements.
*/
inline	real
normal_proj_plane(real* r, const real* v, uint D)
{
	const real n2 = norm_2_sq(v, D);
	if (n2 == 0) return 0;
	const real recipN = 1/sqrt(n2);
	mul(r, v, recipN, D+1);
	return recipN*n2;
}

/** Sets the main diagonal elements of an MxN matrix to the given value, and off-diagonal elements to zero. */
inline	void		set(real* m, real c, uint M, uint N)								{ zero(m, M*N); for (uint i = 0; i < N; ++i, m += M) m[i] = c; }

/** \return a reference the (i,j) element of an M-rowed matrix.  Mutable and const versions. */
inline	real&		elem(real* m, uint i, uint j, uint M)								{ return m[i + M*j]; }
inline	const real&	elem(const real* m, uint i, uint j, uint M)							{ return m[i + M*j]; }

/** \return the head of the N-vector column of an M-rowed matrix m, indexed by colN.  Mutable and const versions. */
inline	real*		col(real* m, uint colN, uint M)										{ return m+M*colN; }
inline	const real*	col(const real* m, uint colN, uint M)								{ return m+M*colN; }

/** Sets a column of an M-rowed matrix. */
inline	void		set_col(real* m, uint colN, const real* v, uint M)					{ cpy(m + colN*M, v, M); }

/** Gets a column of an M-rowed matrix. */
inline	void		get_col(real* r, const real* m, uint colN, uint M)					{ cpy(r, m + colN*M, M); }

/** Sets a row of an MxN matrix. */
inline	void		set_row(real* m, uint rowN, const real* v, uint M, uint N)			{ m += rowN; for (uint i = 0; i < N; ++i, m += M) *m = *v++; }

/** Cets a row of an MxN matrix. */
inline	void		get_row(real* r, const real* m, uint rowN, uint M, uint N)			{ m += rowN; for (uint i = 0; i < N; ++i, m += M) *r++ = *m; }

/** Transpose an NxN matrix in-place. */
inline	void		transpose(real* m, uint N)											{ for (uint j = 0; j < N-1; ++j) { m += j+1; real* t = m+N-1; for (uint i = j+1; i < N; ++i, ++m, t += N) { const real x = *m; *m = *t; *t = x; } } }

/** Transpose an MxN matrix.  m and n cannot point to the same matrix. */
inline	void		transpose(real* r, const real* m, uint M, uint N)					{ for (uint i = 0; i < N; ++i, m += M) set_row(r, i, m, N, M); }

/** Multiply r = m*v, where m is MxN, v is an N-vector, and r is an M-vector.  r and v cannot point to the same vector. */
inline	void		mul(real* r, const real* m, const real* v, uint M, uint N)			{ zero(r, M); for (uint j = 0; j < N; ++j) for (uint i = 0; i < M; ++i) r[i] += (*m++)*v[j]; }

/** Multiply r^T = (v^T)*m, where m is MxN, v is an M-vector, and r is an N-vector.  r and v cannot point to the same vector. */
inline	void		mul_t(real* r, const real* v, const real* m, uint M, uint N)		{ for (uint j = 0; j < N; ++j) { r[j] = 0; for (uint i = 0; i < M; ++i) r[j] += (*m++)*v[i]; } }

/** Multiply r = m*n, where m is MxL, n is LxN , and r is MxN.  r can not point to the same matrix as m or n. */
inline	void		mul(real* r, const real* m, const real* n, uint M, uint L, uint N)	{ for (uint j = 0; j < N; ++j, r += M, n += L) mul(r, m, n, M, L); }

};	// namespace la


#endif // #ifndef _LINALG_H_
