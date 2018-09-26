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

#ifndef _MATH_EXT_H_
#define _MATH_EXT_H_


/*** Extensions to the basic vector operations ***/

#include "math/vec_templates/vec.h"


/* Dimension-specific wedge products */

template<typename F>
inline Vec<F,2>
operator ~ (const Vec<F,2>& a)
{
	Vec<F,2> w;
	w(0) = a(1);
	w(1) = -a(0);
	return w;
}

template<typename F>
inline F
operator ^ (const Vec<F,2>& a, const Vec<F,2>& b)
{
	return a(0)*b(1) - a(1)*b(0);
}

template<typename F>
inline Vec<F,3>
operator ^ (const Vec<F,3>& a, const Vec<F,3>& b)
{
	Vec<F,3> c;
	c(0) = a(1)*b(2) - a(2)*b(1);
	c(1) = a(2)*b(0) - a(0)*b(2);
	c(2) = a(0)*b(1) - a(1)*b(0);
	return c;
}

template<typename F>
inline Vec<F,6>
operator ^ (const Vec<F,4>& a, const Vec<F,4>& b)
{
	Vec<F,6> w;
	w(0) = a(1)*b(0) - a(0)*b(1);
	w(1) = a(2)*b(0) - a(0)*b(2);
	w(2) = a(3)*b(0) - a(0)*b(3);
	w(3) = a(2)*b(1) - a(1)*b(2);
	w(4) = a(3)*b(1) - a(1)*b(3);
	w(5) = a(3)*b(2) - a(2)*b(3);
	return w;
}

template<typename F>
inline Vec<F,4>
operator ^ (const Vec<F,6>& a, const Vec<F,4>& b)
{
	Vec<F,4> w;
	w(0) = - a(5)*b(1) + a(4)*b(2) - a(3)*b(3);
	w(1) =   a(5)*b(0) - a(2)*b(2) + a(1)*b(3);
	w(2) = - a(4)*b(0) + a(2)*b(1) - a(0)*b(3);
	w(3) =   a(3)*b(0) - a(1)*b(1) + a(0)*b(2);
	return w;
}

template<typename F>
inline Vec<F,4>
operator ^ (const Vec<F,4>& a, const Vec<F,6>& b)
{
	return b^a;
}


// Accurate perpendicular to two 3-dimensional vectors
template<typename F>
inline Vec<F,3>	
perp(const Vec<F,3>& v0, const Vec<F,3>& v1)
{
	Vec<F,3> p = v0^v1;	// Cross-product gives perpendicular
	const F s2 = p.length_squared();
	if (s2 != (F)0)	p += ((v0|p)*(p^v1) + (v1|p)*(v0^p))/s2;	// Iterative improvement to (v0 v1)(p) = (0)
	return p;
}

// Accurate perpendicular to three 4-dimensional vectors
template<typename F>
inline Vec<F,4>
perp(const Vec<F,4>& v0, const Vec<F,4>& v1, const Vec<F,4>& v2)
{
	const Vec<F,6> v0v1 = v0^v1;
	Vec<F,4> p = v0v1^v2;	// Double-wedge gives a perpendicular vector
	const F s2 = p.length_squared();
	if (s2 != (F)0)	// Iterative improvement to (v0 v1 v2)|(p) = (0)
	{
		const Vec<F,6> v1v2 = v1^v2;
		const Vec<F,6> v2v0 = v2^v0;
		const F recip_s2 = (F)1/s2;
		p += ((v0|p)*(p^v1v2) + (v1|p)*(p^v2v0) + (v2|p)*(p^v0v1))*recip_s2;	
		p += ((v0|p)*(p^v1v2) + (v1|p)*(p^v2v0) + (v2|p)*(p^v0v1))*recip_s2;
	}
	return p;
}


#endif // #ifndef _MATH_EXT_H_
