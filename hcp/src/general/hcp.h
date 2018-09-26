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

#ifndef _HCP_H_
#define _HCP_H_


/** Implementation of the homogeneous closest point algorithm. */

#include "src/hcp_defs.h"
#include <stddef.h>


/**
Execute the homogeneous closest point algorithm.  This version exposes internal state (thus the 'i' suffix).

\param[in,out]	S			must point to (D+1)*(D+1) reals.  Can be initialized for "warm starts," when N_init > 0 (see below).  It will be filled with the working simplex upon return.
\param[out]		N			the updated number of columns ((D+1)-vectors) in S upon return.
\param[out]		p			must point to D+1 reals.  It will be filled with the solution upon return if the function returns 1.
\param[in]		D			the number of spatial dimensions.  Must be positive.
\param[in]		halfspaces	the user-defined implementation of HCP_Halfspace_Set.
\param[in]		scratch		user-provided scratch space.  It must point to hcp_solve_i_scratch_required(D) bytes.
\param[in]		q			the the objective vector.  It must point to D+1 reals representing a normalized projective point (q[D] == 1) or non-zero direction (q[D] == 0, |q| != 0).
\param[out]		rho			(optional) if not NULL, *rho will be set to the absolute value of the determinant of the linear system solved for the returned solution.  0 < rho <= 1.  Default is NULL.
							NOTE: a crude approximation for the condition number is (1+sqrt(1-rho))/rho.
\param[in]		N_init		(optional) the number of columns ((D+1)-vectors) in the initial simplex S.  Set to zero for a cold start.
\param[in]		tol			(optional) the requested solution accuracy, relative to the size of the output vector.
\return 1, 0, or -1.  The meaning of each value is given below.

Return values:
	 1: the halfspace set has a non-empty intersection
	 0: the halfspace set has an empty intersection
	-1: the function failed to find a result
*/
int hcp_solve_i(real* S, uint& N, real* p, uint D, const HCP_Halfspace_Set& halfspaces, void* scratch, const real* q, real* rho = 0, uint N_init = 0, real tol = hcp_default_tolerance());

/* Scratch space required for hcp_solve_i for a given D > 0. */
inline size_t hcp_solve_i_scratch_required(unsigned D) { return (2*D*D+D+1)*sizeof(real) + (6*D-4)*sizeof(uint); }


/**
Execute the homogeneous closest point algorithm.  This version hides its internal state within scratch memory.
It also does not require a normalized objective vector q (see the requirements for q in hcp_solve_i).

\param[in]		D			the number of spatial dimensions.  Must be positive.
\param[in]		halfspaces	the user-defined implementation of HCP_Halfspace_Set.
\param[in]		scratch		user-provided scratch space.  It must point to hcp_solve_scratch_required(D) bytes.
\param[in,out]	q			(optional) is a pointer the objective vector on input, and (if the function returns 1) the solution vector on output, if it is not NULL.
							If NULL, the objective vector will be the origin (0, ..., 0, 1).  Otherwise, it must point to D+1 reals representing a projective point (q[D] != 0) or direction (q[D] == 0).
\param[out]		rho			(optional) if not NULL, *rho will be set to the absolute value of the determinant of the linear system solved for the returned solution.  0 < rho <= 1.  Default is NULL.
							NOTE: a crude approximation for the condition number is (1+sqrt(1-rho))/rho.
\param[in]		tol			(optional) the requested solution accuracy, relative to the size of the output vector.

Return values:
	 1: the halfspace set has a non-empty intersection
	 0: the halfspace set has an empty intersection
	-1: the function failed to find a result
*/
int hcp_solve(uint D, const HCP_Halfspace_Set& halfspaces, void* scratch, real* q = 0, real* rho = 0, real tol = hcp_default_tolerance());


/** Scratch space required for hcp_solve for a given D > 0. */
inline size_t hcp_solve_scratch_required(uint D) { return hcp_solve_i_scratch_required(D) + (D*D+4*D+3)*sizeof(real); }


/**
This class calls hcp_solve with storage that it manages.
*/
class HCPA
{
public:
	/**
	ctor

	\param[in]	D	the number of spatial dimensions.  Must be positive.
	*/
	HCPA(uint D) : m_D(D), m_storage(new char[hcp_solve_scratch_required(D)]) {}
	~HCPA() { delete [] m_storage; }

	/**
	Execute the homogeneous closest point algorithm (hcp_solve).

	\param[in]		halfspaces	the user-defined implementation of HCP_Halfspace_Set.
	\param[in,out]	q			(optional) is a pointer the objective vector on input, and (if the function returns 1) the solution vector on output, if it is not NULL.
								If NULL, the objective vector will be the origin (0, ..., 0, 1).  Otherwise, it must point to D+1 reals representing a projective point (q[D] != 0) or direction (q[D] == 0).
	\param[out]		rho			(optional) if not NULL, *rho will be set to the absolute value of the determinant of the linear system solved for the returned solution.  0 < rho <= 1.  Default is NULL.
								NOTE: a crude approximation for the condition number is (1+sqrt(1-rho))/rho.
	\param[in]		tol			(optional) the requested solution accuracy, relative to the size of the output vector.

	Return values:
		 1: the halfspace set has a non-empty intersection
		 0: the halfspace set has an empty intersection
		-1: the function failed to find a result
	*/
	int		solve(const HCP_Halfspace_Set& halfspaces, real* q = 0, real* rho = 0, real tol = hcp_default_tolerance()) { return hcp_solve(m_D, halfspaces, m_storage, q, rho, tol); }

	/**
	\return the number of spatial dimensions with which this object is set (and allocated) to perform calculations.
	*/
	uint	get_D() const { return m_D; }

private:
	uint	m_D;		//!< The number of spatial dimensions
	char*	m_storage;	//!< A pointer to the allocated scratch space
};


#endif // #ifndef _HCP_H_
