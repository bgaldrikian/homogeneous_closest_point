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

#ifndef _MAT_H_
#define _MAT_H_


/*** Templated simple matrix class ***/

#include "math/vec_templates/vec.h"


template<typename F, int R, int C>
class Mat
{
public:
					Mat()										{}
					Mat(F v)									{ for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) (*this)(i,j) = (i == j ? v : (F)0); }

	F&				operator ()	(int rowN, int colN)			{ return e[colN][rowN]; }
	const F&		operator ()	(int rowN, int colN)	const	{ return e[colN][rowN]; }

	Vec<F,R>&		operator ()	(int colN)						{ return (Vec<F,R>&)e[colN]; }
	const Vec<F,R>&	operator ()	(int colN)				const	{ return (const Vec<F,R>&)e[colN]; }

	void			set_col(int colN, const Vec<F,R>& col)		{ for (int i = 0; i < R; ++i) (*this)(i,colN) = col(i); }
	Vec<F,R>		get_col(int colN)					const	{ Vec<F,R> col; for (int i = 0; i < R; ++i) col(i) = (*this)(i,colN); return col; }
	void			set_row(int rowN, const Vec<F,C>& row)		{ for (int i = 0; i < C; ++i) (*this)(rowN,i) = row(i); }
	Vec<F,C>		get_row(int rowN)					const	{ Vec<F,C> row; for (int i = 0; i < C; ++i) row(i) = (*this)(rowN,i); return row; }

	Mat<F,R,C>		operator -	()						const	{ Mat<F,R,C> r; for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) r(i,j) = -(*this)(i,j); return r; }

	Mat<F,R,C>		operator *	(F s)					const	{ Mat<F,R,C> r; for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) r(i,j) = (*this)(i,j)*s; return r; }
	Mat<F,R,C>		operator /	(F s)					const	{ return (*this)*((F)1/s); }
	Mat<F,R,C>&		operator *=	(F s)							{ for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) (*this)(i,j) *= s; return *this; }
	Mat<F,R,C>&		operator /=	(F s)							{ *this *= ((F)1/s); return *this; }

	Mat<F,R,C>		operator +	(const Mat<F,R,C>& m)	const	{ Mat<F,R,C> r; for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) r(i,j) = (*this)(i,j)+m(i,j); return r; }
	Mat<F,R,C>		operator -	(const Mat<F,R,C>& m)	const	{ Mat<F,R,C> r; for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) r(i,j) = (*this)(i,j)-m(i,j); return r; }
	Mat<F,R,C>&		operator +=	(const Mat<F,R,C>& m)			{ for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) (*this)(i,j) += m(i,j); return *this; }
	Mat<F,R,C>&		operator -=	(const Mat<F,R,C>& m)			{ for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) (*this)(i,j) -= m(i,j); return *this; }

	Vec<F,R>		operator *	(const Vec<F,C>& v)		const	{ Vec<F,R> r((F)0); for (int j = 0; j < R; ++j) for (int i = 0; i < C; ++i) r(j) += (*this)(j,i)*v(i); return r; }

	template<int S>
	Mat<F,R,S>		operator *	(const Mat<F,C,S>& m)	const	{ Mat<F,R,S> r((F)0); for (int j = 0; j < S; ++j) for (int i = 0; i < R; ++i) for (int k = 0; k < C; ++k) r(i,j) += (*this)(i,k)*m(k,j); return r; }

	Mat<F,C,R>		T()									const	{ Mat<F,C,R> t; for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) t(j,i) = (*this)(i,j); return t; }

private:
	F	e[C][R];
};

// Left scalar multiply
template<typename F, int R, int C>
inline Mat<F,R,C>
operator * (F s, const Mat<F,R,C>& m)
{
	return m*s;
}

// Left row-vector multiply
template<typename F, int R, int C>
inline Vec<F,C>
operator * (const Vec<F,R>& v, const Mat<F,R,C>& m)
{
	Vec<F,C> r((F)0); for (int j = 0; j < C; ++j) for (int i = 0; i < R; ++i) r(j) += v(i)*m(i,j); return r;
}


#endif // #ifndef _MAT_H_
