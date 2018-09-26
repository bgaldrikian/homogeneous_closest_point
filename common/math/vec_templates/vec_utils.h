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

#ifndef _VEC_UTILS_H_
#define _VEC_UTILS_H_


/* Vector math utilities */

#include "math/vec_templates/mat.h"


/* Vector builder */
template<typename F, int D>
inline Vec<F,D+1>
operator << (const Vec<F,D>& a, F b)
{
	return Vec<F,D+1>(a, b);
}


/* Project a point onto a plane */
template<typename F, int D>
inline Vec<F,D>
project(const Vec<F,D>& point, const Vec<F,D+1>& plane)
{
	return point - (Vec<F,D+1>(point, (F)1)|plane)*plane.plane_n();
}


/* Outer product */
template<typename F, int D>
inline Mat<F,D,D>
operator * (const Vec<F,D>& a, const Vec<F,D>& b) { Mat<F,D,D> r; for (int i = 0; i < D; ++i) for (int j = 0; j < D; ++j) r(i,j) = a(i)*b(j); return r; }


// Hodge star in 3D
template<typename F>
inline Mat<F,3,3>
operator * (const Vec<F,3>& v)
{
	Mat<F,3,3> starV;
	starV(0,0) = (F)0;
	starV(0,1) = -v(2);
	starV(0,2) = v(1);
	starV(1,0) = v(2);
	starV(1,1) = (F)0;
	starV(1,2) = -v(0);
	starV(2,0) = -v(1);
	starV(2,1) = v(0);
	starV(2,2) = (F)0;
	return starV;
}

// Rotation in 2D through angle theta
template<typename F>
inline Mat<F,2,2>
rotationMatrix(F theta)
{
	const F c = cos(theta);
	const F s = sin(theta);
	Mat<F,2,2> m;
	m(0,0) = c;
	m(0,1) = -s;
	m(1,0) = s;
	m(1,1) = c;
	return m;
}

// Rotation in 3D about axis (normalized), through angle theta
template<typename F>
inline Mat<F,3,3>
rotationMatrix(F theta, const Vec<F,3>& axis)
{
	const F c = cos(theta);
	return Mat<F,3,3>(c) + ((F)1-c)*(axis*axis) + sin(theta)*(*axis);
}

// Build a homogeneous transform with axis-aligned scaling.  The DxD block will be rotation*scale.
template<typename F, int D>
inline Mat<F,D+1,D+1>
homogeneousTransform(const Mat<F,D,D>& rotation, const Vec<F,D>& scale, const Vec<F,D>& position)
{
	Mat<F,D+1,D+1> tm;
	for (int i = 0; i < D; ++i)
	{
		for (int j = 0; j < D; ++j)
		{
			tm(i,j) = rotation(i,j)*scale(j);
		}
		tm(i,D) = position(i);
	}
	for (int j = 0; j < D; ++j)
	{
		tm(D,j) = (F)0;
	}
	tm(D,D) = (F)1;
	return tm;
}


#endif // #ifndef _VEC_UTILS_H_
