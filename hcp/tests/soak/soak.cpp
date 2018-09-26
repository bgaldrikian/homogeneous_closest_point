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


// Sample utilities
#include "util/sample_utils.h"

// Defines the homogeneous closest point algorithm
#include "src/general/hcp.h"

// Defines HCP_Static_Polytope, an implementation of a set of half-spaces for use with the hcp_solve function.
#include "src/general/sets/hcp_static_polytope.h"

// Defines HCP_Shape_Pair, an implementation of a pair of HCP_HalfSpaceSets for use with the hcp_solve function.
#include "src/general/sets/hcp_shape_pair.h"

// Specialized and optimized versions of HCP
#include "src/specialized/hcp1d.h"
#include "src/specialized/hcp2d.h"
#include "src/specialized/hcp3d.h"
#include "src/specialized/hcp4d.h"
#include "src/specialized/sets/hcpd_static_polytope.h"
#include "src/specialized/sets/hcpd_shape_pair.h"

// Seidel's LP algorithm to compare with HCP results
#include "geom/seidel.h"

// Edge test to use as a tie-breaker in case of differing results
#include "geom/edge_test.h"

// Known effective average radii of various polytopes
#include "geom/radii.h"

// Random distributions
#include "math/dist.h"

// ToI for backup tests
#include "src/general/hcp_toi.h"

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;


static bool marginalVoid(int HCP_result, uint D, uint N, const real* S, const real eps = numeric_limits<real>().epsilon())
{
	if (HCP_result < 0) return false;
	if (HCP_result == 0 && N != D + 1) return false;
	if (HCP_result == 1 && N != D) return false;
	Buffer<real> S_shrink((D + 1)*(D + 1));
	la::cpy(S_shrink, S, (D + 1)*(D + 1));
	for (uint j = 0; j < D + 1; ++j) la::elem(S_shrink, D, j, D + 1) += eps;
	Buffer<real> invS((D + 1)*(D + 1));
	const real detS = lax::inv(invS, S_shrink, D + 1);
	for (uint j = 0; j < D; ++j) if (la::elem(invS, j, D, D + 1) <= 0) return false;
	return true;
}


#define ResizingBuffer(_typename, _buffername, _size)					\
static std::vector<_typename> _buffername##Vec;							\
if (_buffername##Vec.size() < _size) _buffername##Vec.resize(_size);	\
_typename* _buffername = _buffername##Vec.data()


struct toi_result : public pair<int, real>
{
	toi_result() : hcp_result(first), toi(second) {}
	int&	hcp_result;
	real&	toi;
};


static toi_result
calc_toi(uint D, const HCP_Static_Polytope& poly1, const HCP_Static_Polytope& poly2)
{
	toi_result result;
	HCPA_TOI TOI(D);
	ResizingBuffer(real, v1, D);
	ResizingBuffer(real, v2, D);
	la::zero(v1, D);
	la::zero(v2, D);
	v2[0] = 1;
	result.hcp_result = TOI.toi(result.toi, poly1, v1, poly2, v2, true);
	return result;
}


static void
check_results
(
	ostream& out,
	bool& report,
	uint64_t& errorCount,
	uint D,
	uint N,
	uint S_N,
	const real* S,
	uint plane_count,
	const real* planes,
	const HCP_Halfspace_Set& halfspace_set,
	const HCP_Static_Polytope& poly1,
	const HCP_Static_Polytope& poly2,
	const real* objective,
	bool closest_point_test,
	real rho,
	int result_test,
	const char* name_test,
	const real* solution_test,
	int result_compare,
	const char* name_compare,
	const real* solution_compare
)
{
	ResizingBuffer(real, testPlane, D + 1);
	ResizingBuffer(real, testDiff, D);
	ResizingBuffer(char, hcp_solve_scratch, hcp_solve_scratch_required(D));

	int marginal_result_compare = result_compare;
	if (result_test >= 0 && (result_compare >= 0 && result_compare != result_test))
	{
		const bool check_void_simplex = result_test == 0 && S_N == D + 1;	// We'll simply test the plane subset which the test declared to be void in this case
		marginal_result_compare = -1;
		bool use_full_edge_result = D <= 4 && N <= 1000;	// Otherwise this takes too long
		if (check_void_simplex)
		{
			marginal_result_compare = marginalVoid(result_test, D, S_N, S) ? 0 : 1;
			use_full_edge_result = marginal_result_compare != result_test;
		}
		if (use_full_edge_result) marginal_result_compare = edge_test(D, plane_count, planes);
		if (marginal_result_compare != result_test)
		{
			// Use ToI to confirm result
			const toi_result result = calc_toi(D, poly1, poly2);
			if (result.hcp_result != 1 || (result_test == 0 && result.toi < hcp_default_tolerance()) || (result_test == 1 && result.toi > -hcp_default_tolerance()))
			{
				bool within_tolerance = false;
				if (result_test == 1) within_tolerance = halfspace_set.farthest_halfspace(testPlane, solution_test) <= hcp_default_tolerance()*la::norm_2(solution_test, D + 1);
				if (!within_tolerance)
				{
					out << "\nError, " << name_test << " result (" << result_test << ") does not match " << name_compare << " result (" << result_compare << ")";
					if (check_void_simplex || use_full_edge_result) out << " or edge result(" << marginal_result_compare << ")";
					out << ": ";
					out.flush();
					report = true;
					++errorCount;
				}
			}
		}
	}

	if (!closest_point_test)
	{
		if (result_test == 1)
		{
			const real HCP_delta = halfspace_set.farthest_halfspace(testPlane, solution_test);
			if (HCP_delta > (hcp_default_tolerance() + numeric_limits<real>().epsilon())*la::norm_2(solution_test, D + 1))	// numeric_limits<real>().epsilon() term to account for solution normalization roundoff error
			{
				out << "\nError, " << name_test << " solution is not valid.  Max delta = " << HCP_delta << ": ";
				out.flush();
				report = true;
				++errorCount;
			}
			else
			if (result_compare == 1 && halfspace_set.farthest_halfspace(testPlane, solution_compare) <= hcp_default_tolerance()*la::norm_2(solution_compare, D + 1))
			{
				if (solution_test[D] != solution_compare[D])
				{
					if (solution_compare[D] == 0 && la::dot(solution_compare, objective, D) > hcp_default_tolerance()*la::norm_2(solution_compare, D + 1)*la::norm_2(objective, D + 1))
					{
						out << "\nError, " << name_test << " solution is " << (solution_test[D] == 0 ? "un" : "") << "bounded, " << name_compare << " solution is " << (solution_compare[D] == 0 ? "un" : "") << "bounded: ";
						out.flush();
						report = true;
						++errorCount;
					}
				}
				else
				{
					la::sub(testDiff, solution_compare, solution_test, D);
					const real delta = la::dot(testDiff, objective, D) / la::norm_2(objective, D);
					const real tol = 2 * hcp_default_tolerance()*la::norm_2(solution_test, D + 1);
					const real kappa = (1 + sqrt(1 - rho*rho)) / rho;
					if (delta > tol*kappa)
					{
						const toi_result result = calc_toi(D, poly1, poly2);
						if ((result.hcp_result != 1 || std::abs(result.toi) > tol) && !marginalVoid(result_test, D, S_N, S))
						{
							out << "\nError, " << (solution_test[D] == 0 ? "un" : "") << "bounded solution is less optimal than " << name_compare << " solution by " << delta << " (" << delta / tol << " tol): ";
							out.flush();
							report = true;
							++errorCount;
						}
					}
				}
			}
		}
	}
}


static int
soak(uint D, bool closest_point_test, bool unbounded, uint32_t N, uint32_t trialGroupCount, uint32_t trialGroupSize, float overrideOffset = 0, const char* outfile = nullptr, const uint64_t jumpTrialIndex = ~(uint64_t)0, uint verbosity = 1)
{
	ofstream* fout = outfile == nullptr ? nullptr : new ofstream(outfile);			
	ostream& out = fout != nullptr ? *fout : cout;

	Random<real> rnd;

	const uint32_t N1 = N / 2;
	const uint32_t N2 = N - N1;

	const real R = overrideOffset <= 0 ? (real)find_polytope_radius(D, N1) : overrideOffset;
	if (R == 0)
	{
		cout << "\nsoak: polytope radius for (D, N_B) = (" << D << ", " << N1 << ") unknown.\n";
		cout << "Please provide an an override radius.\n\n";
		return -1;
	}

	Buffer<real> planes(N*(D+1));
	HCP_Static_Polytope poly1(D, N1, planes);
	HCP_Static_Polytope poly2(D, N2, planes + N1*(D+1));

	HCP_Shape_Pair shapePair(D, &poly1, &poly2);

	Buffer<real> S((D + 1)*(D + 1));
	Buffer<char> HCP_scratch(hcp_solve_i_scratch_required(D));
	Buffer<real> discrim(D + 1);
	Buffer<real> HCP_objective(D + 1);
	Buffer<real> HCP_solution(D + 1);

	Buffer<real> HCP_solution_old(D + 1);

	SeidelLP<real> SLP;

	Buffer<real> SLP_objective(D + 1);
	Buffer<real> SLP_solution(D + 1);

	Buffer<real> S_HCPD((D + 1)*(D + 1));
	HCPD_Static_Polytope<1> poly1_1(N1, planes);
	HCPD_Static_Polytope<1> poly1_2(N2, planes + N1 * 2);
	HCPD_Static_Polytope<2> poly2_1(N1, planes);
	HCPD_Static_Polytope<2> poly2_2(N2, planes + N1 * 3);
	HCPD_Static_Polytope<3> poly3_1(N1, planes);
	HCPD_Static_Polytope<3> poly3_2(N2, planes + N1 * 4);
	HCPD_Static_Polytope<4> poly4_1(N1, planes);
	HCPD_Static_Polytope<4> poly4_2(N2, planes + N1 * 5);
	HCPD_Shape_Pair<1> shapePair1(&poly1_1, &poly1_2);
	HCPD_Shape_Pair<2> shapePair2(&poly2_1, &poly2_2);
	HCPD_Shape_Pair<3> shapePair3(&poly3_1, &poly3_2);
	HCPD_Shape_Pair<4> shapePair4(&poly4_1, &poly4_2);
	real HCPD_solution[5];

	uint64_t failureCount = 0;
	uint64_t errorCount = 0;

	uint64_t intersectCount = 0;
	uint64_t resultCount = 0;

	out << "\nTesting...\n\n";
	out.flush();

	const bool performSLP = D <= 10;	// Only test against SLP for small D

	if (performSLP) SLP.init(D, N);

	const bool singleTrial = jumpTrialIndex != ~(uint64_t)0;

	const bool displayProgressBar = !singleTrial && verbosity >= 1;

	for (uint32_t groupCount = 0; groupCount < (singleTrial ? 1 : trialGroupCount); ++groupCount)
	{
		if (displayProgressBar)
		{
			out << "\r[................................]";	// group progress bar
			out.flush();
		}

		unsigned lastBarLen = ~0;
		for (uint32_t trial = 0; trial < trialGroupSize; ++trial)
		{
			const uint64_t trialIndex = singleTrial ? jumpTrialIndex : ((uint64_t)trial + (uint64_t)trialGroupSize*(uint64_t)groupCount);

			rnd.seed((uint32_t)(trialIndex & 0xFFFFFFFF));
			rnd.discard(trialIndex >> 32);

			if (!singleTrial)
			{
				if (displayProgressBar)
				{
					const unsigned barLen = (32 * trial) / trialGroupSize;
					if (barLen != lastBarLen)
					{
						char bar[32];
						for (unsigned i = 0; i < barLen; ++i) bar[i] = '#';
						bar[barLen] = '\0';
						out << "\r[" << bar;
						out.flush();
						lastBarLen = barLen;
					}
				}
			}
			else
			{
				out << "Trial index " << trialIndex;
				out.flush();
			}

			la::zero(discrim, D + 1);
			if (unbounded)	rnd.sphere(discrim, D);	// For unbounded tests, require that the plane normals have a positive dot product with a random vector, ensuring an unbounded intersection

			// Create planes
			real* plane = planes;
			for (uint32_t i = 0; i < N; ++i, plane += D + 1)
			{
				rnd.sphere(plane, D);
				plane[D] = -1;
				if (la::dot(plane, discrim, D) < 0) la::neg(plane, plane, D);
			}

			// Create offset transform, apply it to each polytope
			if (!unbounded)	// Only transform for bounded case
			{
				plane = &planes[N1*(D + 1)];
				for (uint32_t i = 0; i < N2; ++i, plane += (D + 1)) plane[D] += plane[0] * 2 * R;
			}

			// Create objective point
			rnd.sphere(HCP_objective, D, pow((real)10, rnd.uniform(-(real)1, (real)1)));
			HCP_objective[D] = closest_point_test ? (real)1 : (real)0;

			// Run HCP
			la::cpy(HCP_solution, HCP_objective, D + 1);
			uint S_N;
			real rho;
			const int HCP_result = hcp_solve_i(S, S_N, HCP_solution, D, shapePair, HCP_scratch, HCP_objective, &rho);
			if (HCP_solution[D] != 0) la::normal_proj_point(HCP_solution, HCP_solution, D);
			else la::normal(HCP_solution, HCP_solution, D);

			// Run HCPD
			bool ran_HCPD = false;
			uint S_N_HCPD;
			real rho_HCPD;
			real* S_HCPD_buffer = (real*)S_HCPD;
			int HCPD_result = -1;
			if (D < 5)
			{
				la::cpy(HCPD_solution, HCP_objective, D + 1);
				ran_HCPD = true;	// Assume true
				switch (D)
				{
				case 1: HCPD_result = hcpd_solve_i<1>(*static_cast<real(*)[][2]>((void*)S_HCPD_buffer), S_N_HCPD, HCPD_solution, shapePair1, HCP_objective, &rho_HCPD);	break;
				case 2: HCPD_result = hcpd_solve_i<2>(*static_cast<real(*)[][3]>((void*)S_HCPD_buffer), S_N_HCPD, HCPD_solution, shapePair2, HCP_objective, &rho_HCPD);	break;
				case 3: HCPD_result = hcpd_solve_i<3>(*static_cast<real(*)[][4]>((void*)S_HCPD_buffer), S_N_HCPD, HCPD_solution, shapePair3, HCP_objective, &rho_HCPD);	break;
				case 4: HCPD_result = hcpd_solve_i<4>(*static_cast<real(*)[][5]>((void*)S_HCPD_buffer), S_N_HCPD, HCPD_solution, shapePair4, HCP_objective, &rho_HCPD);	break;
				default: ran_HCPD = false;	break;
				}
				if (HCPD_solution[D] != 0) la::normal_proj_point(HCPD_solution, HCPD_solution, D);
				else la::normal(HCPD_solution, HCPD_solution, D);
			}

			// Run SLP
			int SLP_result = -1;
			if (performSLP)
			{
				SLP.set_half_spaces(&planes[0], N, D);
				la::zero(SLP_objective, D + 1);
				if (closest_point_test) SLP_objective[0] = 1.0f;
				else la::neg(SLP_objective, HCP_objective, D);
				SLP.set_objective_fn_coefficients(SLP_objective);
				switch (SLP.solve(SLP_solution))
				{
				case 0: SLP_result = 0; break;
				case 1:
				case 2: SLP_result = 1; break;
				default: SLP_result = -1;
				}
				if (abs(SLP_solution[D]) > hcp_default_tolerance()) la::normal_proj_point(SLP_solution, SLP_solution, D);
				else
				{
					la::normal(SLP_solution, SLP_solution, D);
					SLP_solution[D] = 0;
				}
			}

			bool report = false;

			if (HCP_result < 0)
			{
				out << "\nHCP failure: ";
				out.flush();
				report = true;
				++failureCount;
			}

			// Compare HCP with SLP
			if (performSLP)
			{
				check_results(out, report, errorCount, D, N, S_N, S, N, planes, shapePair, poly1, poly2, HCP_objective, closest_point_test, rho, HCP_result, "HCP", HCP_solution, SLP_result, "SLP", SLP_solution);
			}

			// Comapre HCPD with HCP
			if (ran_HCPD)
			{
				if (HCPD_result < 0)
				{
					out << "\nHCPD failure: ";
					out.flush();
					report = true;
					++failureCount;
				}

				check_results(out, report, errorCount, D, N, S_N_HCPD, S_HCPD, N, planes, shapePair, poly1, poly2, HCP_objective, closest_point_test, rho_HCPD, HCPD_result, "HCPD", HCPD_solution, HCP_result, "HCP", HCP_solution);
			}

			if (report)
			{
				out << "trial index " << trialIndex << "\n";
				out << "Continuing...\n";
				if (displayProgressBar) out << "\r[................................]";	// redraw progress bar
				out.flush();
				lastBarLen = ~0;
			}

			intersectCount += (uint64_t)(HCP_result == 1);
			resultCount += (uint64_t)(HCP_result >= 0);

			if (singleTrial) break;
		}

		if (singleTrial) break;

		if (displayProgressBar)
		{
			out << "\r[################################] ";
			out.flush();
		}
	
		const float intersectFraction = resultCount > 0 ? (float)intersectCount / resultCount : 0.0f;

		if (displayProgressBar)
		{
			char output[1000];
			sprintf(output, "%dx%d trials:  %5.1f%% intersection.", groupCount + 1, trialGroupSize, 100.0f*intersectFraction);
			out << output;
			out.flush();
		}
	}

	out << "\n\nTest(s) completed with " << failureCount << " failure(s) and " << errorCount << " error(s).";

	const float intersectFraction = resultCount > 0 ? (float)intersectCount / resultCount : 0.0f;
	char output[1000];
	sprintf(output, "  %5.1f%% intersection.\n", 100.0f*intersectFraction);
	out << output;
	out.flush();

	delete fout;

	return 0;
}


int
main(int argc, char** argv)
{
	int argn = 2;
	if (argc < argn || argv[1][0] == '?')
	{
		cout << "\nPerforms soak tests of the HCP algorithm.\n";
		cout << "\nUsage:\n\n";
		cout << "soak dim [w W] [u U] [n N] [r R] [t trialIndex] [m M] [g G] [v V] [outfile]\n\n";
		cout << "\tdim = the dimensionality of the test.\n";
		cout << "\tW = 0 (linear optimization) or 1 (closest point test).  Default is 1.\n";
		cout << "\tU = 0 or 1, determines if all intersections will be unbounded.  Default is 0.\n";
		cout << "\tN = total number of halfspaces in tests.  Default is 100.\n";
		cout << "\tR = override radius to use.  If 0, uses pre-calculated radii.  Default is 0.\n";
		cout << "\ttrialIndex = the single trial number to run, if given.\n";
		cout << "\tM = number of groups of trials to perform.  Default is 1000.\n";
		cout << "\tG = number of trials per group.  Default is 1000000.\n";
		cout << "\tV = verbosity.  0 = result only, 1 = progress output.  Default is 1.\n";
		cout << "\toutfile = an optional output filename.  If none is given, the standard output is used.\n";
		return -1;
	}

	const int dim = atoi(argv[1]);
	if (dim < 1)
	{
		cout << "\nsoak: first (dim) argument must be positive.\n";
		return -1;
	}

	bool S = false;
	char* outfile = nullptr;
	bool closest_point_test = true;
	bool unbounded = false;
	uint32_t N = 100;
	float R = 0;
	uint32_t trialGroupSize = 1000000;
	uint32_t trialGroupCount = 1000;
	uint64_t trialIndex = ~(uint64_t)0;
	uint verbosity = 1;
	for (; argn < argc; ++argn)
	{
		if (argv[argn][0] == 'w' && argv[argn][1] == '\0')
		{
			if (++argn < argc) closest_point_test = (0 != atoi(argv[argn]));
		}
		else
		if (argv[argn][0] == 'u' && argv[argn][1] == '\0')
		{
			if (++argn < argc) unbounded = (0 != atoi(argv[argn]));
		}
		else
		if (argv[argn][0] == 'n' && argv[argn][1] == '\0')
		{
			if (++argn < argc) N = (uint32_t)atoi(argv[argn]);
		}
		else
		if (argv[argn][0] == 'r' && argv[argn][1] == '\0')
		{
			if (++argn < argc) R = (float)atof(argv[argn]);
		}
		else
		if (argv[argn][0] == 't' && argv[argn][1] == '\0')
		{
			if (++argn < argc) trialIndex = (uint64_t)strtoll(argv[argn], nullptr, 0);
		}
		else
		if (argv[argn][0] == 'm' && argv[argn][1] == '\0')
		{
			if (++argn < argc) trialGroupCount = (uint32_t)atoi(argv[argn]);
		}
		else
		if (argv[argn][0] == 'g' && argv[argn][1] == '\0')
		{
			if (++argn < argc) trialGroupSize = (uint32_t)atoi(argv[argn]);
		}
		else
		if (argv[argn][0] == 'v' && argv[argn][1] == '\0')
		{
			if (++argn < argc) verbosity = (uint)atoi(argv[argn]);
		}
		else
		if (argn < argc)
		{
			outfile = argv[argn];
		}
		else
		{
			cout << "Unknown argument \"" << argv[argn] << "\" ignored.\n";
		}
	}

	return soak((uint)dim, closest_point_test, unbounded, N, trialGroupCount, trialGroupSize, R, outfile, trialIndex, verbosity);
}
