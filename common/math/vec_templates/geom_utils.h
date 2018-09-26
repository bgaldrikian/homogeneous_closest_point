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

#ifndef _GEOM_UTILS_H_
#define _GEOM_UTILS_H_


/* Geometry utilities */

#include "math/vec_templates/mat.h"
#include "math/linalg_ext.h"


// Determinant
template<typename F, int D>
inline F
det(const Mat<F,D,D>& m)
{
	Mat<F,D,D> LU = m;
	uint P[D], Q[D];
	return lax::LUPQ_decompose(P, Q, &LU(0,0), D);
}


// Inverse
template<typename F, int D>
inline Mat<F,D,D>
inv(const Mat<F,D,D>& m)
{
	Mat<F,D,D> LU = m;
	uint P[D], Q[D];
	lax::LUPQ_decompose(P, Q, &LU(0, 0), D);
	Mat<F,D,D> invM;
	lax::LUPQ_invert(&invM(0,0), &LU(0,0), P, Q, D);
	return invM;
}


/* Accurate perpendicular to D-1 D-dimensional vectors, in a function struct for partial specialization */
template<typename F, int D>
struct Perp
{
	Vec<F,D>
	calculate(const Vec<F,D>* v)	// Requires D-1 vectors
	{
		Mat<F,D,D> W;
		for (int i = 0; i < D-1; ++i) W(i) = v[i];
		W(D-1) = (F)0;
		Mat<F,D,D> C;
		lax::cof(&C(0,0), &W(0,0), D);
		const F sign = (F)(((D&1)<<1)-1);
		W(D-1) = sign*C(D-1);
		const F W2 = W(D-1).length_squared();
		if (W2 != (F)0)	// Iterative improvement to (v[0] ... v[D-2])(p) = (0)
		{
			lax::cof(&C(0,0), &W(0,0), D);
			Vec<F,D> q((F)0);
			for (int i = 0; i < D-1; ++i) q += (sign*(W(i)|W(D-1)))*C(i);
			W(D-1) -= q/W2;
		}
		return W(D-1);
	}
};

/* Accurate perpendicular to D-1 D-dimensional vectors, specialized for D = 1 */
template<typename F>	struct Perp<F,1>	{ Vec<F,1>	calculate(const Vec<F,1>*)		{ return Vec<F,1>((F)0); } };

/* Accurate perpendicular to D-1 D-dimensional vectors, specialized for D = 2 */
template<typename F>	struct Perp<F,2>	{ Vec<F,2>	calculate(const Vec<F,2>* v)	{ return ~(*v); } };


/* Accurate perpendicular to D-1 D-dimensional vectors */
template<typename F, int D>	Vec<F,D>	perp(const Vec<F,D>* v)	{ return Perp<F,D>().calculate(v); }	// Requires D-1 vectors


// Cofactor matrix
template<typename F, int D>
Mat<F,D,D>
cof(const Mat<F,D,D>& m)
{
	Vec<F,D> w[D-1];
	for (int minor_col = 0; minor_col < D-1; ++minor_col) w[minor_col] = m(minor_col+1);
	Mat<F,D,D> c;
	c(0) = perp<F,D>(w);
	for (int minor_col_to_replace = 0; minor_col_to_replace < D-1; ++minor_col_to_replace)
	{
		w[minor_col_to_replace] = m(minor_col_to_replace);
		c(minor_col_to_replace+1) = (F)(((minor_col_to_replace&1)<<1)-1)*perp<F,D>(w);
	}
	return c;
}

// Specialization of cof for D = 1
template<typename F>	Mat<F,1,1>	cof(const Mat<F,1,1>& m)	{ return Mat<F,1,1>((F)1); }


// Generalization of cross product in D+1 dimensions
template<typename F, int D>
Vec<F,D+1>
cross_D(const Vec<F,D+1> v[D])
{
	Mat<F,D,D> m;
	// Initialize minor for first row
	for (int minor_row = 0; minor_row < D; ++minor_row)
	{
		for (int minor_col = 0; minor_col < D; ++minor_col)
		{
			m(minor_row, minor_col) = v[minor_col](minor_row+1);
		}
	}
	Vec<F,D+1> p;
	p(0) = det(m);
	for (int minor_row_to_replace = 0; minor_row_to_replace < D; ++minor_row_to_replace)
	{
		for (int minor_col = 0; minor_col < D; ++minor_col)
		{
			m(minor_row_to_replace, minor_col) = v[minor_col](minor_row_to_replace);
		}
		p(minor_row_to_replace+1) = (F)(((minor_row_to_replace&1)<<1)-1)*det(m);
	}
	return p;
}


// Accurate perpendicular to D (D+1)-Dimensional vectors
template<typename F, int D>
bool
perp_D_safe(Vec<F,D+1>& perp, const Vec<F,D+1> v[D])
{
	perp = cross_D<F,D>(v);
	const F s2 = perp.length_squared();
	if (s2 == (F)0) return false;
	for (int iter = 0; iter < 2; ++iter)
	{
		// Iterative improvement to (v[0] ... v[D-1])(perp) = (0)
		Vec<F,D+1> w[D];
		for (int i = 1; i < D; ++i) w[i] = v[i];
		Vec<F,D+1> q((F)0);
		for (int i = 0; i < D; ++i)
		{
			w[i] = perp;
			q += (v[i]|perp)*cross_D<F,D>(w);
			w[i] = v[i];
		}
		perp += q/s2;
	}
	return true;
}


// Cofactor of a square matrix element
template<typename F, int D>
inline F
cof(const Mat<F,D,D>& m, int rowN, int colN)
{
	Mat<F,D-1,D-1> n;
	for (int i = 0, k = 0; i < D; ++i) if (i != rowN) { for (int j = 0, l = 0; j < D; ++j) if (j != colN) { n(k,l) = m(i,j); ++l; } ++k; }
	return (F)((int)1-(int)(((rowN+colN)&1)<<1))*det(n);
}

// Cofactor of a square matrix element - specialization for 1x1 matrix
template<typename F>	inline F	cof(const Mat<F,1,1>&, int, int) { return (F)1; }


template<typename F>
inline bool
intersect_planes(Vec<F,3>& pos, Vec<F,3>& dir, const Vec<F,4>& plane0, const Vec<F,4>& plane1)
{
	// Set up 3x3 set of linear equations
	const Vec<F,3> cross = Vec<F,3>(plane0.plane_n())^Vec<F,3>(plane1.plane_n());
	const F detM = cross.length_squared();
	if (detM == (F)0) return false;

	Mat<F,3,3> M;
	M.set_row(0, plane0.plane_n());
	M.set_row(1, plane1.plane_n());
	M.set_row(2, cross);

	// Create inverse
	Mat<F,3,3> invM;
	invM.set_col(0, Vec<F,3>(plane1.plane_n())^cross);
	invM.set_col(1, cross^Vec<F,3>(plane0.plane_n()));
	invM.set_col(2, cross);
	invM /= detM;

	// rhs for position solution
	Vec<F,3> b;
	b(0) = -plane0.plane_d();
	b(1) = -plane1.plane_d();
	b(2) = (F)0;

	// Solve for pos
	Vec<F,3> x = invM*b;
	x -= invM*(M*x-b);	// Iterate for accuracy
	pos = x;

	// rhs for direction solution
	b(0) = (F)0;
	b(1) = (F)0;
	b(2) = cross.length();

	// Solve for dir
	x = invM*b;
	x -= invM*(M*x-b);	// Iterate for accuracy
	dir = x;

	return true;
}


#endif // #ifndef _GEOM_UTILS_H_
