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

#ifndef _VEC_H_
#define _VEC_H_


/*** Templated simple vector class ***/

#include <cmath>


template<typename F, int D>
class Vec
{
public:
				Vec()										{}
				Vec(F v)									{ for (int i = 0; i < D; ++i) e[i] = v; }
				Vec(const Vec<F,D+1>& v)					{ for (int i = 0; i < D; ++i) e[i] = v(i); }
				Vec(const Vec<F,D>& v, F s)					{ for (int i = 0; i < D-1; ++i) e[i] = v.e[i]; e[D-1] = s; }

				operator Vec<F,D+1>()				const	{ Vec<F,D+1> r; for (int i = 0; i < D; ++i) r(i) = e[i]; r(D) = (F)0; return r; }

	F&			operator ()	(int index)						{ return e[index]; }
	const F&	operator ()	(int index)				const	{ return e[index]; }

	Vec<F,D>	operator -	()						const	{ Vec<F,D> r; for (int i = 0; i < D; ++i) r.e[i] = -e[i]; return r; }

	Vec<F,D>	operator *	(F s)					const	{ Vec<F,D> r; for (int i = 0; i < D; ++i) r.e[i] = e[i]*s; return r; }
	Vec<F,D>	operator /	(F s)					const	{ return (*this)*((F)1/s); }
	Vec<F,D>&	operator *=	(F s)							{ for (int i = 0; i < D; ++i) e[i] *= s; return *this; }
	Vec<F,D>&	operator /=	(F s)							{ *this *= ((F)1/s); return *this; }

	Vec<F,D>	operator +	(const Vec<F,D>& v)		const	{ Vec<F,D> r; for (int i = 0; i < D; ++i) r.e[i] = e[i]+v.e[i]; return r; }
	Vec<F,D>	operator -	(const Vec<F,D>& v)		const	{ Vec<F,D> r; for (int i = 0; i < D; ++i) r.e[i] = e[i]-v.e[i]; return r; }
	Vec<F,D>&	operator +=	(const Vec<F,D>& v)				{ for (int i = 0; i < D; ++i) e[i] += v.e[i]; return *this; }
	Vec<F,D>&	operator -=	(const Vec<F,D>& v)				{ for (int i = 0; i < D; ++i) e[i] -= v.e[i]; return *this; }

	F			operator |	(const Vec<F,D>& v)		const	{ F r = (F)0; for (int i = 0; i < D; ++i) r += e[i]*v.e[i]; return r; }

	F			length_squared()					const	{ return (*this)|(*this); }
	F			length()							const	{ return sqrt(length_squared()); }

	Vec<F,D>	normal()							const	{ Vec<F,D> r = *this; r.vector_normalize(); return r; }

	Vec<F,D>	plane_n()							const	{ Vec<F,D> r; for (int i = 0; i < D-1; ++i) r.e[i] = e[i]; r.e[D-1] = (F)0; return r; }
	F			plane_d()							const	{ return e[D-1]; }

	F			vector_normalize()
				{
					const F l2 = length_squared();
					const F recipL = l2 > (F)0 ? (F)1/sqrt(l2) : (F)0;
					*this *= recipL;
					return recipL*l2;
				}

	F			point_normalize()
				{
					const F w = e[D-1];
					if (w != (F)0) *this /= w;
					return w;
				}

	F			plane_normalize()
				{
					const F l2 = plane_n().length_squared();
					const F recipL = l2 > (F)0 ? (F)1/sqrt(l2) : (F)0;
					*this *= recipL;
					return recipL*l2;
				}

private:
	F	e[D];
};

// Left scalar multiply
template<typename F, int D>
inline Vec<F,D>
operator * (F s, const Vec<F,D>& m)
{
	return m*s;
}


#endif // #ifndef _VEC_H_
