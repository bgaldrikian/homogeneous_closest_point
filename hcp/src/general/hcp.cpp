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


/*** Implementation of the homogeneous closest point algorithm ***/

#include "src/general/hcp.h"
#include "src/general/hcp_linalg.h"
#include "math/discrete.h"

#include <limits>


/*
Support function for hcp_test_feature.

S must contain N plane equations, stored as (D+1)-vectors.  Each vector is of the form (n^T d)^T, where n is a D-vector.  Points x in the corresponding plane satisfy (n^T)x + d = 0.
q must contain a D-vector.

The result p satisfies (S^T)p = 0, where the right-hand side is an N-dimensional 0 vector, and be the nearest vector to q.

q_D is the w-component of q.
If q_D != 0, the interpretation is that p is the nearest point to q that lies in the intersection of all planes stored in S.
If q_D == 0, the interpretation is that p is the vector closest to q that is perpendicular to the normals of all planes stored in S.

scratch must be D*D*sizeof(real) + 2*D*sizeof(uint) bytes.

returns determinant of linear system used for solution.
*/
inline real
hcp_solve_feature(real* p, const real* q, const real* S, const uint* columns, uint D, uint N, void* scratch)
{
	// Set up local storage from scratch
	real* H = (real*)scratch;
	uint* P = (uint*)(H + D * D);
	uint* Q = P + D;	// Pi and Qi really only need to be (D-1) in size, but we do this to keep the D = 0 case from being problematic

	// Handle N = 0 case
	if (N == 0)
	{
		la::cpy(p, q, D + 1);
		if (la::norm_2_sq(p, D + 1) == 0) p[D] = 1;
		return 1;
	}

	// Copy normals into H
	for (uint j = 0; j < N; ++j) la::set_col(H, j, la::col(S, columns[j], D + 1), D);

	// Create complement basis for H
	real detH;
	if ((detH = hcpla::QRP_decompose(P, H, D, N, 0)) == 0) return 0;
	hcpla::QR_create_complement_basis(H, D, N);

	if ((p[D] = q[D]) == 0)
	{
		// In this case we simply project, p = (CC^T)q, where C is the complement basis for H
		la::zero(p, D);
		const real* h = la::col(H, N, D);
		for (uint j = N; j < D; ++j, h += D) la::madd(p, h, la::dot(q, h, D), p, D);
		if (la::norm_2_sq(p, D) != 0) return detH;	// If p != 0, return solution.  Otherwise solve for a bounded (p[D] = 1) solution below
		p[D] = 1;
	}

	// Set up rhs for complement part here, before transposing
	for (uint j = N; j < D; ++j) p[j] = la::dot(q, la::col(H, j, D), D);

	// Replace first N columns of H with normals, and transpose for linear equation setup
	for (uint j = 0; j < N; ++j) la::set_col(H, j, la::col(S, columns[j], D + 1), D);
	la::transpose(H, D);

	// LU decompose
	detH = hcpla::LUPQ_decompose(P, Q, H, D);
	if (detH == 0) return 0;

	// Set up rhs and solve
	for (uint j = 0; j < N; ++j) p[j] = -p[D] * la::elem(S, D, columns[j], D + 1);
	hcpla::LUPQ_solve(p, H, P, Q, D);

	return detH;
}


/*
Support function for hcp_search_features.

Finds the closest point on the given combination of planes to the target q, and checks to see if it meets the constraints given by the remaining plane equations.

If the solution meets the constraints, it is checked against the current solution to see if it is more optimal, and if so, replaces the current solution.

scratch must be (D*D+D+1)*sizeof(real) + 2*D*sizeof(uint) bytes
*/
inline bool
hcp_test_feature(real* p, real& p2_min, real& det_basis, uint* columns, uint& N_new, const uint* columns_test, uint N_test, uint N_old, const real* S, const real* q, uint D, void* scratch)
{
	// Set up local storage from scratch
	real* p_test = (real*)scratch;	// test solution
	scratch = p_test + (D + 1);

	// Solve for closest point on the feature defined by column indices
	const real det_basis_test = hcp_solve_feature(p_test, q, S, columns_test, D, N_test, scratch);
	if (det_basis_test == 0) return false;

	// Make sure solution is valid (not outside the other half-spaces in S)
	for (Complement<uint> j(columns_test, N_test); j < N_old; ++j) if (la::dot(p_test, la::col(S, j, D + 1), D + 1) > 0) return false;

	// Don't accept a bounded solution if we have an unbounded solution already
	if (p_test[D] > p[D]) return true;

	// If the solution isn't better than the current one, don't use it
	const real p2 = q[D] != 0 ? la::diff_norm_2_sq(p_test, q, D) : -la::dot(p_test, q, D);
	if (p[D] == p_test[D] && p2 >= p2_min) return true;

	// Use this solution
	la::cpy(p, p_test, D + 1);
	p2_min = p2;
	det_basis = det_basis_test;
	memcpy(columns, columns_test, N_test * sizeof(uint));
	N_new = N_test;

	return true;
}


/**
Searches for solutions based off of combinations of plane equations from the columns of S.

S must be a (D+1)xN matrix with column vectors which are plane equations.
q is a (D+1)-dimensional homogeneous target vector.
scratch must be (D*D+D+1)*sizeof(real) + (5*D+2)*sizeof(uint) bytes

Upon return, the best solution is written to p, a user-supplied array of D+1 reals, and the determinant of the
linear system used to find this solution is written to det_basis.  The indices of the columns of S which are used
to find p are written to 'columns', and N is replaced by the number of columns written to 'columns'.

If a solution is found, this function returns true, otherwise it returns false.
*/
inline bool
hcp_search_features(uint* columns, real* p, real& det_basis, const real* S, uint& N, const real* q, uint D, void* scratch)
{
	// Set up local storage from scratch
	uint* indices = (uint*)scratch;	// Temporary index list
	uint* column_set = indices + D;
	uint* columns_test = column_set + (D + 1);
	scratch = columns_test + (D + 1);

	// The original number of working planes (old planes)
	const uint column_count = N - 1;

	// Initialize column indices
	for (uint i = 0; i < column_count; ++i) column_set[i] = i;

	// Initialize measure for best solution
	real p2_min = std::numeric_limits<real>().max();
	p[D] = 1;	// Initialize projective component for solution type comparison

	// Bit used to mark columns in column index set.  This assumes there will never be more than 2^31 columns
	constexpr uint kept_bit = 1 << (std::numeric_limits<uint>().digits - 1);

	// Start with the largest number of halfspaces in the subset, and decrease
	uint column_set_size = column_count;
	uint N_test = column_count >= D ? D : column_count;
	do
	{
		for (Choose<uint> c(column_set_size, N_test, indices); c; ++c)
		{
			for (uint i = 0; i < N_test; ++i) columns_test[i] = column_set[indices[i]] & ~kept_bit;
			if (hcp_test_feature(p, p2_min, det_basis, columns, N, columns_test, N_test, column_count, S, q, D, scratch)) for (uint i = 0; i < N_test; ++i) column_set[indices[i]] |= kept_bit;
		}
		const uint old_column_set_size = column_set_size;
		column_set_size = 0;
		for (uint j = 0; j < old_column_set_size; ++j) if (column_set[j] & kept_bit) column_set[column_set_size++] = column_set[j] & ~kept_bit;
		if (old_column_set_size == N_test && column_set_size == 0) column_set_size = old_column_set_size;
		if (N_test > column_set_size) N_test = column_set_size;
	} while (N_test--);

	// If p2_min is not updated, we have not found a solution
	return p2_min < std::numeric_limits<real>().max();
}


/*
Simplex update function.

Given a simplex matrix S with N columns of (D+1)-dimensional plane vectors, and a target homogeneous (D+1)-vector q,

If Poly(S) is not empty then:

* p is set to the closest (D+1)-vector (in the first D coordinates) to the q in the set Poly(S).  
* S is updated to contain only the columns with plane vectors for planes that contain the point p
* N is updated to the number of columns in the updated matrix S
* det_basis is set to the determinant of the basis used to perform the calculation which returns p
* the function returns true

If Poly(S) is empty, then
* p is unchanged
* the function returns false

scratch must point to (2*D*D+D+1)*sizeof(real) + (6*D-4)*sizeof(uint) bytes
*/
inline bool
hcp_update(real* p, real& det_basis, real* S, uint& N, const real* q, uint D, void* scratch)
{
	// Handle N == 0 case
	if (N == 0)
	{
		la::cpy(p, q, D + 1);
		return true;
	}

	// Set up local storage from scratch
	uint* columns = (uint*)scratch;			// Column indices used to update simplex matrix columns
	real* v = (real*)(columns + (D - 1));	// Householder vector
	real* R = v + D;						// Reduced simplex matrix
	real* qR = R + D*D;						// Reduced objective
	scratch = qR + D;

	// Save off index of the new plane, which will serve as the original number of working planes (old planes)
	const uint N_old = N - 1;

	// Restrict solution to new plane, since it must lie there if it exists.  Using a Householder reflection which maps the new plane normal to the last coordinate axis:
	const real* h_new = la::col(S, N_old, D + 1);
	const real Vnd = h_new[D - 1];
	const real sgnVnd = Vnd >= 0 ? (real)1 : -(real)1;
	la::cpy(v, h_new, D);
	v[D - 1] += sgnVnd;
	const real negBeta = -1 / (1 + Vnd*sgnVnd);

	// Project all plane equations and the objective into the new plane
	real* Rn = R;
	const real* Sn = S;
	const real sgnVnd_d_new = sgnVnd*h_new[D];
	for (uint n = 0; n < N_old; ++n, Rn += D, Sn += D+1)
	{
		la::madd(Rn, v, negBeta*la::dot(v, Sn, D), Sn, D);	// Calculate the [D-1] component of the vector, since it is used to calculate the plane displacement
		Rn[D - 1] = Sn[D] + sgnVnd_d_new*Rn[D - 1];	// Replace Rn[D-1] with the correct plane displacement for the reduced plane in the space of the new plane (note, these won't be normalized, which is OK)
	}
	la::madd(qR, v, negBeta*la::dot(v, q, D), q, D - 1);
	qR[D - 1] = q[D];

	// Search for solutions in reduced space
	const bool solution_found = hcp_search_features(columns, p, det_basis, R, N, qR, D - 1, scratch);

	// Translate p to original space
	p[D] = p[D - 1];
	p[D - 1] *= sgnVnd_d_new;
	la::madd(p, v, negBeta*la::dot(v, p, D), p, D);

	// If a solution was found, update the simplex
	if (solution_found)
	{
		for (uint j = 0; j < N; ++j) la::cpy(la::col(S, j, D + 1), la::col(S, columns[j], D + 1), D + 1);
		la::cpy(la::col(S, N++, D + 1), h_new, D + 1);	// Include last column
	}

	return solution_found;
}


/*
scratch must be (2*D*D+D+1)*sizeof(real) + (6*D-4)*sizeof(uint) bytes.
*/
int
hcp_solve_i(real* S, uint& N, real* p, uint D, const HCP_Halfspace_Set& halfspaces, void* scratch, const real* q, real* rho, uint N_init, real tol)
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
		if (!hcp_update(p, det_basis, S, N, q, D, scratch)) { result = 0; break; }	// Void found
		const real delta = halfspaces.farthest_halfspace(la::col(S, N, D + 1), p);
		if (delta <= 0 || delta*delta <= tol_sq*la::norm_2_sq(p, D + 1)) { result = 1; break; }	// Intersection found
		++N;	// Add half-space to S
	}

	// If robustness is given, fill it with the absolute value of the determinant of the basis matrix
	if (rho) *rho = det_basis > 0 ? det_basis : -det_basis;

	return result;
}


/*
scratch must be (3*D*D+5*D+4)*sizeof(real) + (6*D-4)*sizeof(uint) bytes
*/
int
hcp_solve(uint D, const HCP_Halfspace_Set& halfspaces, void* scratch, real* q, real* rho, real tol)
{
	// Set up local storage from scratch
	real* S = (real*)scratch;
	real* p = S + (D + 1)*(D + 1);
	real* objective = p + (D + 1);
	scratch = objective + (D + 1);

	// Objective = q if it is not NULL, otherwise it is the origin represented in homogeneous coordinates
	if (q) { if (q[D] != 0) la::normal_proj_point(objective, q, D); else la::cpy(objective, q, D + 1); }
	else { la::zero(objective, D); objective[D] = 1; }

	uint N;
	const int result = hcp_solve_i(S, N, p, D, halfspaces, scratch, objective, rho, 0, tol);

	// If q is given, fill it with the solution
	if (q) la::cpy(q, p, D + 1);

	return result;
}
