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

#ifndef _HCPD_MATH_H_
#define _HCPD_MATH_H_


#include "types.h"


/** Math and vector operations templated for use by D-specific specializations of the homogeneous closest point algorithm, for D = 2, 3, and 4. */

/** Scalar functions. */
inline real sq(real x) { return x*x; }
inline void swap(real& x, real& y) { const real t = x; x = y; y = t; }

/** Swap D-vectors x and y. */
template<uint D> inline void swap(real* x, real* y) { struct V { real e[D]; }; const V t = *(V*)x; *(V*)x = *(V*)y; *(V*)y = t; }

/** Divide a D-vector by a scalar: r = u/c.  Doing this instead of multiplying by 1/c in order to increase accuracy. */
template<uint D> inline void div(real* r, const real* u, real c);	// Only use specialized versions
template<> inline void div<2>(real* r, const real* u, real c) { r[0] = u[0] / c; r[1] = u[1] / c; }
template<> inline void div<3>(real* r, const real* u, real c) { r[0] = u[0] / c; r[1] = u[1] / c; r[2] = u[2] / c; }
template<> inline void div<4>(real* r, const real* u, real c) { r[0] = u[0] / c; r[1] = u[1] / c; r[2] = u[2] / c;  r[3] = u[3] / c; }
template<> inline void div<5>(real* r, const real* u, real c) { r[0] = u[0] / c; r[1] = u[1] / c; r[2] = u[2] / c;  r[3] = u[3] / c;  r[4] = u[4] / c; }

/** zero a D-vector v. */
template<uint D> inline void zero(real* v);	// Only use specialized versions
template<> inline void zero<2>(real* v) { struct V { real x, y; inline V() : x(0), y(0) {} }; *(V*)v = V(); }
template<> inline void zero<3>(real* v) { struct V { real x, y, z; inline V() : x(0), y(0), z(0) {} }; *(V*)v = V(); }
template<> inline void zero<4>(real* v) { struct V { real x, y, z, u; inline V() : x(0), y(0), z(0), u(0) {} }; *(V*)v = V(); }
template<> inline void zero<5>(real* v) { struct V { real x, y, z, u, v; inline V() : x(0), y(0), z(0), u(0), v(0) {} }; *(V*)v = V(); }

/** Copy a D-vector w into v. */
template<uint D> inline void cpy(real* v, const real* w) { struct V { real e[D]; }; *(V*)v = *(V*)w; }

/** Calculate the square distance (2-norm) between two D-vectors v and w. */
template<uint D> inline real diff_norm_2_sq(const real* v, const real* w);	// Only use specialized versions
template<> inline real diff_norm_2_sq<1>(const real* v, const real* w) { return sq(v[0] - w[0]); }
template<> inline real diff_norm_2_sq<2>(const real* v, const real* w) { return sq(v[0] - w[0]) + sq(v[1] - w[1]); }
template<> inline real diff_norm_2_sq<3>(const real* v, const real* w) { return sq(v[0] - w[0]) + sq(v[1] - w[1]) + sq(v[2] - w[2]); }

/** Multiply a D-vector by a scalar and add another D-vector: r = u*c + v. */
template<uint D> inline void madd(real* r, const real* u, real c, const real* v);	// Only use specialized versions
template<> inline void madd<1>(real* r, const real* u, real c, const real* v) { r[0] = u[0] * c + v[0]; }
template<> inline void madd<2>(real* r, const real* u, real c, const real* v) { r[0] = u[0] * c + v[0]; r[1] = u[1] * c + v[1]; }
template<> inline void madd<3>(real* r, const real* u, real c, const real* v) { r[0] = u[0] * c + v[0]; r[1] = u[1] * c + v[1]; r[2] = u[2] * c + v[2]; }
template<> inline void madd<4>(real* r, const real* u, real c, const real* v) { r[0] = u[0] * c + v[0]; r[1] = u[1] * c + v[1]; r[2] = u[2] * c + v[2]; r[3] = u[3] * c + v[3]; }

/** Dot product of D-vectors v and w. */
template<uint D> inline real dot(const real* v, const real* w);	// Only use specialized versions
template<> inline real dot<2>(const real* v, const real* w) { return v[0] * w[0] + v[1] * w[1]; }
template<> inline real dot<3>(const real* v, const real* w) { return v[0] * w[0] + v[1] * w[1] + v[2] * w[2]; }
template<> inline real dot<4>(const real* v, const real* w) { return v[0] * w[0] + v[1] * w[1] + v[2] * w[2] + v[3] * w[3]; }
template<> inline real dot<5>(const real* v, const real* w) { return v[0] * w[0] + v[1] * w[1] + v[2] * w[2] + v[3] * w[3] + v[4] * w[4]; }


#endif // #ifndef _HCPD_MATH_H_
