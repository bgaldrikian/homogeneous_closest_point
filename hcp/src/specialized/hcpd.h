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

#ifndef _HCPD_H_
#define _HCPD_H_


#include "src/hcp_defs.h"
#include "src/specialized/hcpd_math.h"


/*****************************************************************************************
NOTE: There is no need to include this file directly.

Instead, include one or more of the specializations hcp1d.h, hcpd2.h, hcp3d.h, or hcp4d.h.

A specialization is needed in order for the code to compile.
******************************************************************************************/


/** Templated HCP functions for dimension-specific implementations of the homogeneous closest point algorithm. */

/**
Support function for hcpd_test_feature<D>.

This function is not defined for general D.  Specializations are defined in hcp1d.h, hcp2d.h, hcp3d.h, and hcp4d.h.

An implementation if this function must find the nearest vector (p) to an input vector (q) which solves:

	dot<D+1>(p, S[indices[i]]) = 0

for all i = 0, ... , N-1.

\param[out]	p		the solution, a (D+1)-dimensional vector.
\param[in]	q		the target, a (D+1)-dimensional vector.
\param[in]	S		an array of (D+1)-dimensional vectors, of at least the size such that each value in indices[0], ..., indices[N-1] is a valid element.
					In this algorithm we may assume that S stores a (D+1)x(D+1) array.
\param[in]	indices	stores the indices elements of the S array (equivalently, columns of the S matrix) used.
\param[in]	N		the size of the indices array.

\return the square of the determinant of the linear system solved to find p.  If 0 is returned, a solution is not found.
*/
template<uint D>
inline real
hcpd_solve_feature(real p[], const real q[], const real S[][D + 1], const uint indices[], uint N);


/*
Support function for hcpd_search_features<D>.

Given N_set plane equations stored in columns 0 through (N_set-1) of S, and a subset of columns of S defined by the bits in columns_test,
a solution p to the equations

	dot<D+1>(p, S[c]) = 0

for all test column indices c, is found.  If that solution also satisfies

	dot<D+1>(p, S[c]) <= 0

for all indices 0 <= c < N_set which are _not_ marked by bits in columns_test, then the solution p is valid.

Valid solutions are compared with the current solution by values stored in the measure array.  If a valid solution has better measure parameters,
then it is closer to the current solution and the current solution is copied to the output vector p along with update measure parameters,
the square determinant of the basis used to solve for p, and the columns used to calculate p.
*/
template<uint D>
inline bool
hcpd_test_feature(real p[D + 1], real& p2_min, real& det_basis, uint& columns, uint columns_test, uint N_set, const real S[D + 1][D + 1], const real q[D + 1])
{
	real p_test[D + 1];	// test solution

	// Create index array
	uint indices[D];
	uint N = 0;
	for (uint j = 0; j < N_set; ++j) if (((columns_test >> j) & 1) != 0) indices[N++] = j;

	// Solve for closest point on the feature defined by indices
	const real det_basis_test = hcpd_solve_feature<D>(p_test, q, S, indices, N);
	if (det_basis_test == 0) return false;

	// Don't accept a bounded solution if we have an unbounded solution already
	if (p_test[D] > p[D]) return false;

	// Make sure solution is valid (not outside the other half-spaces in S)
	for (uint j = 0; j < N_set; ++j) if ((((columns_test >> j) & 1) == 0) && dot<D + 1>(p_test, S[j]) > 0) return false;

	// If the solution isn't better than the current one, don't use it
	const real p2 = q[D] != 0 ? diff_norm_2_sq<D>(p_test, q) : -dot<D + 1>(p_test, q);
	if (p[D] == p_test[D] && p2 >= p2_min) return false;

	// Use this solution
	cpy<D + 1>(p, p_test);
	p2_min = p2;
	det_basis = det_basis_test;
	columns = columns_test;

	return true;
}


/*
Support function for hcpd_update<D>.

This function is not defined for general D.  Specializations are defined in hcp1d.h, hcp2d.h, hcp3d.h, and hcp4d.h.

An implementation if this function must find the best solution p for every subset of the columns
of S, an N x (D+1) matrix with columns holding plane equations, such that 

	dot<D+1>(p, S[c]) = 0

for all test column indices c, and

	dot<D+1>(p, S[c]) <= 0

for all other columns c of S.

\param[out]	columns		the columns used to find the best solution, represented as bits in an unsigned integer.
\param[out]	p			the best solution, a (D+1)-dimensional vector.
\param[out] det_basis	the determinant of the linear system solved to find p.
\param[in]	S			an array of N (D+1)-dimensional vectors representing plane equations.
\param[in]	N			the number of columns in S.
\param[in]	q			the target, a (D+1)-dimensional vector.

\return true iff a solution is found.
*/
template<uint D>
inline bool
hcpd_search_features(uint& columns, real p[], real& det_basis, const real S[][D + 1], uint N, const real q[]);


/*
Support function for hcpd_solve_i<D>.

Finds the closest point represented by a homogeneous (D+1)-vector p to the target (D+1)-vector q, which satisfies

	dot<D+1>(p, S[c]) <= 0

for all columns c of S, an N x (D+1) matrix with columns holding (D+1)-dimensional vectors representing plane equations.

On output, only the columns of S which held strict equality in the above equation are kept.  If no solution is found, the function returns false.

\param[out]		p			the solution, a (D+1)-dimensional vector.
\param[out]		det_basis	the determinant of the linear system solved to find p.
\param[in,out]	S			an array of N (D+1)-dimensional vectors representing plane equations.  On output, only the columns of S which lead to strict equality in the above equation are kept.
\param[in,out]	N			the number of columns in S.  On output, the number of columns of S kept.
\param[in]		q			the target, a (D+1)-dimensional vector.

\return true iff a solution is found.
*/
template<uint D>
inline bool
hcpd_update(real p[], real& det_basis, real S[][D + 1], uint& N, const real q[])
{
	// Handle N == 0 case
	if (N == 0)
	{
		cpy<D + 1>(p, q);
		return true;
	}

	// Reduce problem into new plane
	real v[D];		// Householder vector
	real R[D][D];	// Reduced simplex matrix
	real qR[D];		// Reduced objective

	// Save off index of the new plane, which will serve as the original number of working planes (old planes)
	const uint N_old = N - 1;

	// Restrict solution to new plane, since it must lie there if it exists.  Using a Householder reflection which maps the new plane normal to the last coordinate axis:
	const real* h_new = S[N_old];
	const real Vnd = h_new[D - 1];
	const real sgnVnd = Vnd >= 0 ? (real)1 : -(real)1;
	cpy<D>(v, h_new);
	v[D - 1] += sgnVnd;
	const real negBeta = -1 / (1 + Vnd*sgnVnd);

	// Project all plane equations and the objective into the new plane
	const real sgnVnd_d_new = sgnVnd*h_new[D];
	for (uint n = 0; n < N_old; ++n)
	{
		madd<D>(R[n], v, negBeta*dot<D>(v, S[n]), S[n]);	// Perform the D-dimensional transformation of the normal vector
		R[n][D - 1] = S[n][D] + sgnVnd_d_new*R[n][D - 1];	// Replace Rn[D-1] with the correct plane displacement for the reduced plane in the space of the new plane (note, these won't be normalized, which is OK)
	}
	madd<D - 1>(qR, v, negBeta*dot<D>(v, q), q);
	qR[D - 1] = q[D];

	uint columns;		// Bits representing columns of S used in feature calculations
	const bool solution_found = hcpd_search_features<D-1>(columns, p, det_basis, R, N_old, qR);

	// Translate p to original space
	p[D] = p[D - 1];
	p[D - 1] *= sgnVnd_d_new;
	madd<D>(p, v, negBeta*dot<D>(v, p), p);

	// If no combination of planes produces a feature, we have a void simplex
	if (!solution_found) return false;

	// ... otherwise update the simplex
	const uint N_set = N;
	N = 0;
	for (uint j = 0; j < N_set; ++j, columns >>= 1) if (columns & 1) cpy<D + 1>(S[N++], S[j]);
	cpy<D + 1>(S[N++], h_new);	// Include last column

	return true;
}


/**
Execute the homogeneous closest point algorithm (templated D).  This version exposes internal state (thus the 'i' suffix).

\param[in,out]	S			will be filled with the working simplex upon return.  Can be initialized for "warm starts," when N_init > 0 (see below)
\param[out]		N			the updated number of columns ((D+1)-vectors) in S upon return.
\param[out]		p			will be filled with the solution upon return if the function returns 1.
\param[in]		halfspaces	the user-defined implementation of HCP_Halfspace_Set.
\param[in]		q			the objective vector.  It must contain 5 reals representing a normalized projective point (q[4] == 1) or direction (q[4] == 0, |q| == 1).
\param[out]		rho			(optional) if not NULL, *rho will be set to the absolute value of the determinant of the linear system solved for the returned solution.  0 < rho <= 1.  Default is NULL.
							NOTE: a crude approximation for the condition number is (1+sqrt(1-rho))/rho.
\param[in]		N_init		(optional) the number of columns ((D+1)-vectors) in the initial simplex S.  Set to zero for a cold start.
\param[in]		tol			the requested solution accuracy, relative to the size of the output vector.
\return 1, 0, or -1.  The meaning of each value is given below.

Return values:
	 1: the halfspace set has a non-empty intersection
	 0: the halfspace set has an empty intersection
	-1: the function failed to find a result
*/
template<uint D>
inline int
hcpd_solve_i(real S[D + 1][D + 1], uint& N, real p[D + 1], const HCP_Halfspace_Set& halfspaces, const real q[D + 1], real* rho = nullptr, uint N_init = 0, real tol = hcp_default_tolerance())
{
	// Constants
	const uint max_iteration_count = 20 * D;	// Maximum allowed iterations of main loop.  If exceeded, error value (-1) is returned
	const real tol_sq = tol*tol;	// Using tolerance squared

	// Initialize
	N = N_init;

	// Iterate until a stopping condition is met or the maximum number of iterations is reached
	real det_basis = 1;	// Used to approximate solution robustness
	int result = -1;
	for (uint i = 0; i < max_iteration_count; ++i)
	{
		if (!hcpd_update<D>(p, det_basis, S, N, q)) { result = 0; break; }	// Void found
		const real delta = halfspaces.farthest_halfspace(S[N], p);
		if (delta <= 0 || delta*delta <= tol_sq*dot<D + 1>(p, p)) { result = 1; break; }	// Intersection found
		++N;	// Add half-space to S
	}

	// If robustness is given, fill it with the absolute value of the determinant of the basis matrix
	if (rho) *rho = det_basis > 0 ? det_basis : -det_basis;

	return result;
}


/*
Execute the homogeneous closest point algorithm (templated D).  This version hides its internal state.
It also does not require a normalized objective vector q (see the requirements for q in hcpd_solve_i<D>).

\param[in]		halfspaces	user-defined implementation of HCP_Halfspace_Set.
\param[in,out]	q			a pointer the objective vector on input, and (if the function returns 1) the solution vector on output, if it is not NULL.
							If NULL, the objective vector will be the origin (0, ..., 0, 1).  Otherwise, it must point to D+1 reals representing a projective point (q[D] != 0) or direction (q[D] == 0).
\param[out]		rho			(optional) if not NULL, *rho will be set to the absolute value of the determinant of the linear system solved for the returned solution.  0 < rho <= 1.  Default is NULL.
							NOTE: a crude approximation for the condition number is (1+sqrt(1-rho))/rho.
\param[in]		tol			the requested solution accuracy, relative to the size of the output vector.
\return 1, 0, or -1.		The meaning of each value is given below.

Return values:
	 1: the halfspace set has a non-empty intersection
	 0: the halfspace set has an empty intersection
	-1: the function failed to find a result
*/
template<uint D>
inline int
hcpd_solve(const HCP_Halfspace_Set& halfspaces, real* q = nullptr, real* rho = nullptr, real tol = hcp_default_tolerance())
{
	real S[D + 1][D + 1];
	real p[D + 1];
	real objective[D + 1];

	// Objective = q if it is not NULL, otherwise it is the origin represented in homogeneous coordinates
	if (q)
	{
		if (q[D] != 0) div<D + 1>(objective, q, q[D]);
		else cpy<D + 1>(objective, q);
	}
	else
	{
		zero<D + 1>(objective);
		objective[D] = 1;
	}

	uint N;
	const int result = hcpd_solve_i<D>(S, N, p, halfspaces, objective, rho, 0, tol);

	// If q is given, fill it with the solution
	if (q) cpy<D + 1>(q, p);

	return result;
}


#endif // #ifndef _HCPD_H_
