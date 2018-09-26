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

#ifndef _DIST_H_
#define _DIST_H_

#include <random>
#include "math/linalg.h"


/** Various random distributions for a templated floating-point type F. */
template<typename F>
class Random
{
public:
			Random(uint32_t s = std::mt19937::default_seed) : m_generator(s) {}

	void	seed(uint32_t s)								{ m_generator.seed(s); }

	void	discard(uint64_t skip)							{ m_generator.discard(skip); }

	F		real()											{ return m_uniform(m_generator); }

	/** P(x) = 1/(maximum-minimum), minimum <= x < maximum. */
	F		uniform(F minimum, F maximum)					{ return (maximum-minimum)*real() + minimum; }

	/** P(x) = e^-x, 0 <= x < infinity. */
	F		exponential()									{ return m_exponential(m_generator); }

	/** P(x) = (1/2pi)e^(-(1/2)x^2), -infinity < x < infinity. */
	F		gaussian()										{ return m_gaussian(m_generator); }

	/** D-dimensional gaussian distribution. */
	void	gaussian(F* v, uint32_t D)						{ for (size_t i = 0; i < D; ++i) v[i] = gaussian(); }

	/** Point uniformly distributed on the surface of a D-ball of radius R. */
	void	sphere(F* v, uint32_t D, F R = (F)1)			{ do gaussian(v, D); while (la::normal(v, v, D) == 0); la::mul(v, v, R, D); }

	/** Point uniformly distributed within the volume of a D-ball of radius R. */
	void	ball(F* v, uint32_t D, F R = (F)1)				{ sphere(v, D, R*pow(real(), (F)1/(F)D)); }

	/** Point uniformly distributed within the volume of a D-spherical shell with inner radius Ri and outer radius Ro. */
	void	spherical_shell(F* v, uint32_t D, F Ri, F Ro)	{ F R; while ((R = Ro*pow(real(), (F)1 / (F)D)) < Ri) {} sphere(v, D, R); }

private:
	std::mt19937						m_generator;
	std::uniform_real_distribution<F>	m_uniform;
	std::exponential_distribution<F>	m_exponential;
	std::normal_distribution<F>			m_gaussian;
};


#endif // #ifndef _DIST_H_
