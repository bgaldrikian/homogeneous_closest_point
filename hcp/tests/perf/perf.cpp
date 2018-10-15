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

#include <vector>
#include <iostream>
#include <cstdlib>
#include <memory.h>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <numeric>
#include <stdarg.h>
using namespace std;


// Sample utilities
#include "util/sample_utils.h"

// Defines the homogeneous closest point algorithm
#include "src/general/hcp.h"

// Defines HCP_Static_Polytope, an implementation of a set of half-spaces for use with the hcp_solve function.
#include "src/general/sets/hcp_static_polytope.h"

// Random distributions
#include "math/dist.h"

// Seidel's LP algorithm to compare with HCP results
#include "geom/seidel.h"

// Known effective average radii of various polytopes
#include "geom/radii.h"

// D-specialized and optimized versions of HCP
#include "src/specialized/hcp1d.h"
#include "src/specialized/hcp2d.h"
#include "src/specialized/hcp3d.h"
#include "src/specialized/hcp4d.h"
#include "src/specialized/sets/hcpd_static_polytope.h"
#include "src/specialized/sets/hcpd_shape_pair.h"


ArrayWithString(uint, autoD, 2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 100);

ArrayWithString(uint, autoN, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000);


class MessageFilter
{
public:
	MessageFilter() : m_verbosity(0) {}
	MessageFilter(int verbosity) : m_verbosity(verbosity) {}

	void	set_verbosity(int verbosity) { m_verbosity = verbosity; }
	int		get_verbosity() const { return m_verbosity; }

	void	print(int msg_level, const char* format, ...) const
	{
		if (msg_level > m_verbosity) return;
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		fflush(stdout);
	}

private:
	int	m_verbosity;
};


// Output data
template<typename F>
struct Stats
{
	Stats() : mean(0), median(0), rms_mean(0), mad_mean(0), rms_median(0), mad_median(0) {}

	F	mean;
	F	median;
	F	rms_mean;
	F	mad_mean;
	F	rms_median;
	F	mad_median;
};

typedef Stats<float>	StatsF;


struct StatsOutput
{
	enum Stat
	{
		kNone,
		kMean,
		kMedian,

		kStatCount
	};

	enum Dev
	{
		kRMS,
		kMAD,

		kDevCount
	};

	StatsOutput() : stat(kNone) {}

	Stat	stat;
	Dev		dev;
};


template<typename F>
static Stats<F>
calculate_stats(const vector<F>& values)
{
	Stats<F> stats;

	if (values.size())
	{
		// Calculate mean and std_dev
		const F mean = accumulate(values.begin(), values.end(), (F)0) / values.size();
		stats.mean = mean;

		// Calculate median
		vector<F> s = values;
		sort(s.begin(), s.end());
		const F median = (s[(s.size() - 1) / 2] + s[s.size() / 2]) / 2;
		stats.median = median;

		// Calculate rms deviation from mean
		vector<F> v = values;
		transform(v.begin(), v.end(), v.begin(), [mean](F x) { return x - mean; });
		stats.rms_mean = sqrt(inner_product(v.begin(), v.end(), v.begin(), (F)0) / max(v.size() - 1, (size_t)1));	// With Bessel's correction

		// Calculate rms deviation from median
		v = values;
		transform(v.begin(), v.end(), v.begin(), [median](F x) { return x - median; });
		stats.rms_median = sqrt(inner_product(v.begin(), v.end(), v.begin(), (F)0) / max(v.size() - 1, (size_t)1));	// With Bessel's correction

		// Calculate median absolute deviation from mean
		transform(s.begin(), s.end(), v.begin(), [mean](F x) { return abs(x - mean); });
		sort(v.begin(), v.end());
		stats.mad_mean = (v[(v.size() - 1) / 2] + v[v.size() / 2]) / 2;

		// Calculate median absolute deviation from median
		transform(s.begin(), s.end(), v.begin(), [median](F x) { return abs(x - median); });
		sort(v.begin(), v.end());
		stats.mad_median = (v[(v.size() - 1) / 2] + v[v.size() / 2]) / 2;
	}

	return stats;
}


struct ObjectiveType
{
	enum
	{
		kDirection,
		kPoint,

		kCount
	};
};


const char* ObjectiveTypeNames[] =
{
	"LO",
	"CP"
};


typedef tuple<string, string, uint, uint>	SKey;

struct TestSetup
{
	uint			D;
	uint			N;
	vector<uint>	hAxisVals;
	vector<string>	seriesNames;
};


struct TestParams
{
	TestParams() {}
	TestParams(uint in_objType, uint in_D, uint in_N, uint in_g, uint in_m) : objType(in_objType), D(in_D), N(in_N), g(in_g), m(in_m) {}

	uint objType;
	uint D;
	uint N;
	uint g;
	uint m;
};


class TestGroup
{
public:
	TestGroup(const MessageFilter* msg = nullptr) : m_D(0), m_N(0), m_M(0), m_msg(msg) {}
	virtual ~TestGroup() {}

	virtual const char*	name() const = 0;
	virtual const char*	arr_name() const = 0;

	/**
		The maximum number of spatial dimensions this test should run.  Default is numeric_limits<uint>().max().
	*/
	virtual uint		max_D() { return numeric_limits<uint>().max(); }

	/**
		The default number of arrangements per group to run for this test.
	*/
	virtual uint		default_M(uint /*D*/, uint /*N*/) { return 100; }

	/**
		The default number of groups to run for this test.
	*/
	virtual uint		default_G(uint /*D*/, uint /*N*/) { return 10000; }

	/**
		Sets internal storage
	*/
	void				set_params(uint D, uint N, uint M)
						{
							m_D = m_N = m_M = 0;
							if (N & 1) { if (m_msg) m_msg->print(0, "Error: N must be even.\n"); return; }
							if (D > max_D()) return;
							m_D = D;
							m_N = N;
							m_M = M;
							if (!v_alloc()) { m_D = m_N = m_M = 0; }
							m_retvals.resize(m_M);
						}

	/**
		planes must point to M*N*(D+1) reals representing N planes for each arrangement (expected to be tangent to a unit sphere with arbitrary center)
	*/
	virtual void		set_planes(const real* planes) = 0;

	/**
		Determine if an objective type is supported.  Default assumption is that all are supported.  Override if this is not true.
	*/
	virtual bool		objective_type_supported(uint objType) { return objType <= ObjectiveType::kCount; }

	/**
		objective must point to M*(D+1) reals representing a homogeneous objective vector for each arrangement
	*/
	virtual void		set_objectives(const real* objectives) = 0;

	/**
		If arr_num >= M(), averages the time to run all M() arrangements and returns time in microseconds.
		If arr_num < M(), runs only arrangement arr_num (return value is arbitrary for this case).
	*/
	virtual float		run() = 0;

	bool				valid() const { return m_D > 0; }

	uint				N() const { return m_N; }
	uint				M() const { return m_M; }
	uint				D() const { return m_D; }

	const int*			retvals() const { return m_retvals.data(); }

	string				pack_test_id(const TestParams& params) const
						{
							return
								string("t") + name() + ";" +
								string("a") + arr_name() + ";" +
								string("w") + to_string(params.objType) + ";" +
								string("D") + to_string(params.D) + ";" +
								string("N") + to_string(params.N) + ";" +
								string("g") + to_string(params.g) + ";" +
								string("m") + to_string(params.m);
						}

	bool				unpack_test_id(TestParams& params, const char* test_id) const
						{
							string testName;
							string arrName;
							params.objType = params.g = params.m = ~0;
							params.D = params.N = 0;
							istringstream iss(test_id);
							for (string token; std::getline(iss, token, ';');)
							{
								const char* param = token.c_str();
								const char p_tok = *param++;
								switch (p_tok)
								{
								case 't': testName = string(param);				break;
								case 'a': arrName = string(param);				break;
								case 'w': params.objType = (uint)atoi(param);	break;
								case 'D': params.D = (uint)atoi(param);			break;
								case 'N': params.N = (uint)atoi(param);			break;
								case 'g': params.g = (uint)atoi(param);			break;
								case 'm': params.m = (uint)atoi(param);			break;
								default: if (m_msg) m_msg->print(1, "Unknown parameter token %c ignored.\n", p_tok);
								}
							}
							return testName == name() && arrName == arr_name() && params.objType < 2 && params.D && params.N && ~params.g && ~params.m;
						}

protected:
	virtual bool		v_alloc() = 0;

	Clock					m_timer;
	uint					m_N;
	uint					m_M;
	uint					m_D;
	vector<int>				m_retvals;
	const MessageFilter*	m_msg;
};


template<typename AG>	// AG = Arrangement Generator
class HCP_TestGroup : public TestGroup
{
public:
	HCP_TestGroup(const MessageFilter* msg = nullptr) : TestGroup(msg) {}

	virtual const char*	arr_name() const override { return AG().name(); }

	virtual uint
	default_M(uint D, uint /*N*/)
	{
		if (D <= 3) return 100;
		if (D <= 7) return 10;
		return 1;
	}

	virtual uint
	default_G(uint D, uint /*N*/)
	{
		if (D <= 10) return 10000;
		if (D <= 30) return 1000;
		return 100;
	}

	virtual void
	set_planes(const real* planes) override
	{
		AG().run(m_planes, planes, m_D, m_N, m_M, m_msg);
	}

	virtual void
	set_objectives(const real* objectives) override
	{
		m_objectives = objectives;
	}

protected:
	virtual bool
	v_alloc() override
	{
		m_planes.resize(m_M*m_N*(m_D + 1));
		return true;
	}

	const real*		m_objectives;
	MemBuf<real,16>	m_planes;
};


template<typename AG>	// AG = Arrangement Generator
class HCPGeneral_TestGroup : public HCP_TestGroup<AG>
{
public:
	HCPGeneral_TestGroup(const MessageFilter* msg = nullptr) : HCP_TestGroup<AG>(msg) {}

	virtual const char*	name() const override { return "HCP"; }

	virtual float
	run() override
	{
		if (!this->valid()) return 0.0f;
		char* scratch = m_scratch;
		real* S = m_S;
		real* p = m_p;
		const HCP_Static_Polytope* polytope = m_polytopes.data();
		uint S_N;
		const real* objective = this->m_objectives;
		int* retval = this->m_retvals.data();
		const long long startTime = this->m_timer.ticks();
		for (uint m = this->m_M; m--; objective += this->m_D + 1) *retval++ = hcp_solve_i(S, S_N, p, this->m_D, *polytope++, scratch, objective);
		const long long endTime = this->m_timer.ticks();
		return 1000000.0f*(float)(endTime - startTime) / (float)(this->m_timer.ticks_per_second()*this->m_M);	// Time in microseconds
	}

private:
	virtual bool
	v_alloc() override
	{
		HCP_TestGroup<AG>::v_alloc();
		m_S.resize((this->m_D + 1)*(this->m_D + 1));
		m_p.resize(this->m_D + 1);
		m_scratch.resize(hcp_solve_i_scratch_required(this->m_D));
		m_polytopes.resize(this->m_M);
		real* arr = this->m_planes;
		for (size_t m = 0; m < this->m_M; ++m, arr += this->m_N*(this->m_D + 1)) m_polytopes[m] = HCP_Static_Polytope(this->m_D, this->m_N, arr);
		return true;
	}

	MemBuf<real,16>				m_S;
	MemBuf<real,16>				m_p;
	MemBuf<char,16>				m_scratch;
	vector<HCP_Static_Polytope>	m_polytopes;
};


template<uint DIM, typename AG>	// AG = Arrangement Generator
class HCPD_TestGroup : public HCP_TestGroup<AG>
{
public:
	HCPD_TestGroup(const MessageFilter* msg = nullptr) : HCP_TestGroup<AG>(msg) {}

	virtual const char*	name() const override { return s_name.c_str(); }

	virtual float
	run() override
	{
		if (!this->valid()) return 0.0f;
		const HCPD_Static_Polytope<DIM>* polytope = m_polytopes.data();
		real S[DIM + 1][DIM + 1];
		uint S_N;
		real p[DIM + 1];
		const real* objective = this->m_objectives;
		int* retval = this->m_retvals.data();
		const long long startTime = this->m_timer.ticks();
		for (uint m = this->m_M; m--; objective += this->m_D + 1) *retval++ = hcpd_solve_i<DIM>(S, S_N, p, *polytope++, objective);
		const long long endTime = this->m_timer.ticks();
		return 1000000.0f*(float)(endTime - startTime) / (float)(this->m_timer.ticks_per_second()*this->m_M);	// Time in microseconds
	}

private:
	virtual bool
	v_alloc() override
	{
		if (this->m_D != DIM) { if (this->m_msg) this->m_msg->print(0, "Error: D must equal templated value %D.\n", DIM); return false; }
		HCP_TestGroup<AG>::v_alloc();
		m_polytopes.resize(this->m_M);
		real* arr = this->m_planes;
		for (size_t m = 0; m < this->m_M; ++m, arr += this->m_N*(this->m_D + 1)) m_polytopes[m] = HCPD_Static_Polytope<DIM>(this->m_N, arr);
		return true;
	}

	vector< HCPD_Static_Polytope<DIM> >	m_polytopes;

	static const string					s_name;
};

template<uint DIM, typename AG>	// AG = Arrangement Generator
const string	HCPD_TestGroup<DIM, AG>::s_name = string("HCP") + to_string(DIM);


template<typename AG>	// AG = Arrangement Generator
class Seidel_TestGroup : public TestGroup
{
public:
	Seidel_TestGroup(const MessageFilter* msg = nullptr) : TestGroup(msg) {}

	virtual const char*	name() const override { return "Seid"; }
	virtual const char*	arr_name() const override { return AG().name(); }

	virtual uint		max_D() { return 30; }

	virtual uint
	default_M(uint D, uint /*N*/)
	{
		if (D <= 3) return 100;
		if (D <= 7) return 10;
		return 1;
	}

	virtual uint
	default_G(uint D, uint /*N*/)
	{
		if (D <= 5) return 10000;
		if (D <= 10) return 1000;
		if (D <= 20) return 100;
		return 10;
	}

	virtual void
	set_planes(const real* planes) override
	{
		AG().run(m_planes, planes, m_D, m_N, m_M, m_msg);
		real* arr = m_planes;
		for (SeidelLP<real>& slp : m_SeidelLPs)
		{
			slp.set_half_spaces(arr, m_N, m_D);
			arr += m_N*(m_D + 1);
		}
	}

	virtual bool		objective_type_supported(uint objType) { return objType == ObjectiveType::kDirection; }

	virtual void
	set_objectives(const real* objectives) override
	{
		for (SeidelLP<real>& slp : m_SeidelLPs)
		{
			slp.set_objective_fn_coefficients(objectives);
			objectives += m_D + 1;
		}
	}

	virtual float
	run() override
	{
		if (!this->valid()) return 0.0f;
		int* retval = m_retvals.data();
		const long long startTime = this->m_timer.ticks();
		SeidelLP<real>* slp = m_SeidelLPs.data();
		for (uint m = this->m_M; m--; ++slp) *retval++ = slp->is_feasible();
		const long long endTime = this->m_timer.ticks();
		return 1000000.0f*(float)(endTime - startTime) / (float)(this->m_timer.ticks_per_second()*this->m_M);	// Time in microseconds
	}

private:
	virtual bool
	v_alloc() override
	{
		m_planes.resize(m_M*m_N*(m_D + 1));

		m_SeidelLPs.resize(m_M);
		for (SeidelLP<real>& slp : m_SeidelLPs) slp.init(m_D, m_N);

		return true;
	}

	MemBuf<real,16>				m_planes;
	vector< SeidelLP<real> >	m_SeidelLPs;
};


struct PolyArrangementType
{
	enum
	{
		kUnbounded,
		kCentered,
		kMarginal,
		kSeparate,

		kPolyArrangementCount
	};
};


static const char* ArrangementTypeName(uint type)
{
	switch (type)
	{
	case PolyArrangementType::kUnbounded:	return "UNB";
	case PolyArrangementType::kCentered:	return "OVL";
	case PolyArrangementType::kMarginal:	return "MRG";
	case PolyArrangementType::kSeparate:	return "SEP";
	}

	return "unknown";
}


template<uint ArrangementType>
struct Poly_AG
{
	static const char*	name() { return ArrangementTypeName(ArrangementType); }

	uint				type() const { return ArrangementType; }

	bool
	run(real* out, const real* in, uint D, uint N, uint M, const MessageFilter* msg)
	{
		if (D == 0) return false;

		real Delta;
		switch (ArrangementType)
		{
		case PolyArrangementType::kUnbounded:
		case PolyArrangementType::kCentered:	Delta = 0.0f;	break;
		case PolyArrangementType::kMarginal:	Delta = 1.0f;	break;
		case PolyArrangementType::kSeparate:	Delta = 10.0f;	break;
		default: return false;
		}

		const real R = (real)find_polytope_radius(D, N / 2);
		if (Delta != 0 && R <= 0)	// Allow R == 0 if Delta == 0
		{
			if (msg) msg->print(1, "[Poly_Pair_AG] D = %d, N_B = %d: R* unknown.\n", D, N / 2);
			return false;
		}

		la::cpy(out, in, M*N*(D + 1));
		const real offset = static_cast<real>(R*Delta);

		if (offset != 0)
		{
			real* plane = out;
			const uint N_B = N / 2;
			for (size_t m = 0; m < M; ++m)
			{
				for (size_t n = 0; n < N_B; ++n, plane += D + 1) plane[D] += offset*plane[0];
				for (size_t n = 0; n < N_B; ++n, plane += D + 1) plane[D] -= offset*plane[0];
			}
		}

		if (ArrangementType == PolyArrangementType::kUnbounded)
		{
			real* plane = out;
			for (size_t i = 0; i < M*N; ++i, plane += D + 1) if (plane[0] >= 0) la::neg(plane, plane, D);
		}

		return true;
	}
};


static string create_data_filename(const string& file_prefix, const TestGroup& test, size_t objN, uint N, uint D)
{
	return file_prefix + test.arr_name() + "_" + test.name() + "_" + ObjectiveTypeNames[objN] + "_N" + to_string(N) + "_D" + to_string(D);
}


static FILE* open_data_file_w(const string& data_dir, const string& filename, const MessageFilter& msg)
{
	const string& path = data_dir + "/" + filename;
	FILE* f = fopen(path.c_str(), "a");
	if (f == nullptr) f = fopen(path.c_str(), "w");	// Try "w" mode if "a" does not work
	if (f == nullptr) msg.print(0, "\nCould not open data file for write.  Be sure directory %s exists.\n", data_dir);
	return f;
}


static FILE* open_data_file_r(const string& data_dir, const string& filename, const MessageFilter& msg)
{
	const string& path = data_dir + "/" + filename;
	FILE* f = fopen(path.c_str(), "r");
	if (f == nullptr) msg.print(1, "\nCould not open data file %s for read.\n", path.c_str());
	return f;
}


static void
dual_printf(FILE* fp, const char* format, ...)
{
	va_list args;
	va_start(args, format);
	vprintf(format, args);
	fflush(stdout);
	if (fp) vfprintf(fp, format, args);
	va_end(args);
}


static void
output_results(const TestSetup& setup, const map<SKey,StatsF>& stats, const StatsOutput& stats_out, const string& data_dir, const MessageFilter& msg)
{
	for (uint arrN = 0; arrN < PolyArrangementType::kPolyArrangementCount; ++arrN)
	{
		FILE* f = nullptr;
		const string filename = data_dir + "/" + (setup.D ? "N" : "D") + "_" + ArrangementTypeName(arrN);
		if ((f = fopen(filename.c_str(), "w")) == nullptr)
		{
			msg.print(0, "Could not create stats file %s.\n", filename);
			f = nullptr;
		}
		dual_printf(f, "# ");
		if (setup.D) dual_printf(f, "D = %d, N", setup.D);
		else dual_printf(f, "N = %d, D", setup.N);
		dual_printf(f, " sweep, arrangement type = %s.  Times in microseconds.\n", ArrangementTypeName(arrN));
		dual_printf(f, "%-8s", setup.D ? "N" : "D");
		for (const string& series_name : setup.seriesNames) dual_printf(f, "%-11s%-11s", series_name.c_str(), (string("d(") + series_name + ")").c_str());
		dual_printf(f, "\n");
		for (uint val : setup.hAxisVals)
		{
			uint D, N;
			if (setup.D)
			{
				D = setup.D;
				N = val;
			}
			else
			{
				N = setup.N;
				D = val;
			}
			dual_printf(f, "%-8d", val);
			for (const string& series_name : setup.seriesNames)
			{
				map<SKey, StatsF>::const_iterator it = stats.find(SKey(string(ArrangementTypeName(arrN)), series_name, N, D));
				if (it != stats.end())
				{
					float stat = 0;
					float dev = 0;
					switch (stats_out.stat)
					{
					case StatsOutput::kMean:
						stat = it->second.mean;
						switch (stats_out.dev)
						{
						case StatsOutput::kRMS: dev = it->second.rms_mean;	break;
						case StatsOutput::kMAD: dev = it->second.mad_mean;	break;
						}
						break;
					case StatsOutput::kMedian:
						stat = it->second.median;
						switch (stats_out.dev)
						{
						case StatsOutput::kRMS: dev = it->second.rms_median;	break;
						case StatsOutput::kMAD: dev = it->second.mad_median;	break;
						}
						break;
					}
					if (stat > 0) dual_printf(f, "%-11.5g%-11.5g", stat, dev);
					else dual_printf(f, "                      ");
				}
			}
			dual_printf(f, "\n");
		}
		if (f != nullptr) fclose(f);
		printf("\n");
	}
}


/**
	D = # of spatial dimensions
	N = # (total) of half-spaces per arrangement
	M = # of arrangements per measurement group
	G = # of measurement groups
*/
static bool
run(map<SKey,StatsF>& stats, const TestSetup& setup, const string& data_dir, bool compile, const vector<TestGroup*>& tests, uint D, uint N, uint M_in, uint G_in, uint start_group, TestParams* params, bool flush_data, const string& file_prefix, const MessageFilter& msg)
{
	for (TestGroup* test : tests)
	{
		const uint M = M_in ? M_in : test->default_M(D, N);
		test->set_params(D, N, M);
	}

	const real R = (real)find_polytope_radius(D, N / 2);
	const real point_Ri = 2 * R;
	const real point_Ro = 10 * R;

	uint M_max = 0;
	for (TestGroup* test : tests) M_max = max(M_max, test->M());

	MemBuf<real,16> poly_planes(M_max*N*(D + 1));
	MemBuf<real,16> point_objectives(M_max*(D + 1));
	MemBuf<real,16> direction_objectives(M_max*(D + 1));

	const real* objective_lists[ObjectiveType::kCount] = { direction_objectives, point_objectives };

	// [tests.size(), 2, G]
	vector< vector< vector<float> > > t_meas(tests.size());
	vector<uint> G(tests.size());
	uint G_max = 0;
	for (size_t testN = 0; testN < tests.size(); ++testN)
	{
		vector< vector<float> >& o = t_meas[testN];
		o.resize(ObjectiveType::kCount);
		G[testN] = G_in ? G_in : tests[testN]->default_G(D, N);
		for (vector<float>& g : o) g.resize(G[testN]);
		G_max = max(G_max, G[testN]);
	}

	if (!compile)
	{
		vector<uint> objTypes;
		if (!params)
		{
			objTypes.resize(ObjectiveType::kCount);
			iota(objTypes.begin(), objTypes.end(), 0);
		}
		else objTypes = { params->objType };

		Random<real> rnd;

		/*
			log_2(D) : 4 bits (D <= 65535)
			log_2(N) : 5 bits
			g		 : 14 bits (G <= 16384)
			m		 : 7 bits (M <= 128)
		*/
		const uint32_t seedDN = ilogb(D) << 26 | ilogb(N) << 21;

		msg.print(2, "groups (%%):");
		int last_print = -1;
		for (uint g = start_group; g < G_max; ++g)
		{
			const int print_n = (100 * g) / G_max;
			if (print_n != last_print)
			{
				msg.print(2, " %d", print_n);
				last_print = print_n;
			}
			real* poly_plane = poly_planes;
			real* point = point_objectives;
			real* direction = direction_objectives;
			for (uint m = 0; m < M_max; ++m, point += D + 1, direction += D + 1)
			{
				const uint32_t seed = seedDN | (!params ? (g << 7 | m) : (params->g << 7 | params->m));
				rnd.seed(seed);
				rnd.sphere(point, D, rnd.uniform(point_Ri, point_Ro)); point[D] = 1;
				rnd.sphere(direction, D); direction[D] = 0;
				for (uint n = 0; n < N; ++n, poly_plane += D + 1) { rnd.sphere(poly_plane, D); poly_plane[D] = -1; }
			}

			for (uint testN = 0; testN < tests.size(); ++testN)
			{
				TestGroup& test = *tests[testN];
				if (!test.valid()) continue;
				if (g >= G[testN]) break;
				test.set_planes(poly_planes);
				for (uint objType : objTypes)
				{
					if (!test.objective_type_supported(objType)) continue;
					test.set_objectives(objective_lists[objType]);
					const float t = test.run();
					t_meas[testN][objType][g] = t;
					if (!params)
					{
						for (uint m = 0; m < test.M(); ++m)
						{
							if (test.retvals()[m] < 0)
							{
								const string id = test.pack_test_id(TestParams(objType, D, N, g, m));
								msg.print(0, "\n[%s_%s] obj = %s, D = %d, N = %d, group = %d, arr = %d: error.  Test ID: %s\n", test.name(), test.arr_name(), ObjectiveTypeNames[objType], D, N, g, m, id.c_str());
							}
						}
					}
					if (flush_data)
					{
						const string filename = create_data_filename(file_prefix, test, objType, N, D);
						FILE* f = open_data_file_w(data_dir, filename, msg);
						if (f == nullptr) return false;
						fprintf(f, "%.9g\n", t);
						fclose(f);
					}
				}
			}
		}
		msg.print(2, " %d", 100);
	}

	if (compile || !flush_data)	// Only write to data file here if we weren't doing it after each measurement, or use the same loops to read data.
	{
		for (size_t testN = 0; testN < tests.size(); ++testN)
		{
			TestGroup* test = tests[testN];
			for (size_t objN = 0; objN < ObjectiveType::kCount; ++objN)
			{
				const string series_name = string(test->name()) + "/" + ObjectiveTypeNames[objN];
				if (find(setup.seriesNames.begin(), setup.seriesNames.end(), series_name) != setup.seriesNames.end())
				{
					const string filename = create_data_filename(file_prefix, *test, objN, N, D);
					if (!compile)
					{
						FILE* f = open_data_file_w(data_dir, filename, msg);
						if (f == nullptr) return false;
						for (float t : t_meas[testN][objN]) fprintf(f, "%.9g\n", t);
						fclose(f);
					}
					else
					{
						vector<float>& t = t_meas[testN][objN];
						fill(t.begin(), t.end(), 0.0f);
						FILE* f = open_data_file_r(data_dir, filename, msg);
						if (f != nullptr)
						{
							t.clear();	// Get actual size from file
							float read_val;
							while (1 == fscanf(f, "%f", &read_val)) t.push_back(read_val);
							fclose(f);
						}
						stats[SKey(string(test->arr_name()), series_name, N, D)] = calculate_stats(t);
					}
				}
			}
		}
	}

	return true;
}


static bool
setup_and_run(const string& run_type, const string& data_dir, const StatsOutput& stats_out, uint M_in, uint G_in, uint start_group, bool exclusivelyHCP, bool flush_data, const MessageFilter& msg)
{
	char sweep_type = '\0';
	if (run_type == "N") sweep_type = 'N';
	else
	if (run_type == "D") sweep_type = 'D';

	// Create tests
	vector<TestGroup*> tests;
	tests.push_back(new HCPGeneral_TestGroup< Poly_AG<PolyArrangementType::kCentered> >(&msg));
	tests.push_back(new HCPGeneral_TestGroup< Poly_AG<PolyArrangementType::kMarginal> >(&msg));
	tests.push_back(new HCPGeneral_TestGroup< Poly_AG<PolyArrangementType::kSeparate> >(&msg));
	tests.push_back(new HCPGeneral_TestGroup< Poly_AG<PolyArrangementType::kUnbounded> >(&msg));
	if (sweep_type != 'D')	// Load all tests if sweeping N or running from a test ID.  For the latter case, we will search for the correct test.
	{
		tests.push_back(new HCPD_TestGroup< 3, Poly_AG<PolyArrangementType::kCentered> >(&msg));
		tests.push_back(new HCPD_TestGroup< 3, Poly_AG<PolyArrangementType::kMarginal> >(&msg));
		tests.push_back(new HCPD_TestGroup< 3, Poly_AG<PolyArrangementType::kSeparate> >(&msg));
		tests.push_back(new HCPD_TestGroup< 3, Poly_AG<PolyArrangementType::kUnbounded> >(&msg));
	}
	if (!exclusivelyHCP)
	{
		tests.push_back(new Seidel_TestGroup< Poly_AG<PolyArrangementType::kCentered> >(&msg));
		tests.push_back(new Seidel_TestGroup< Poly_AG<PolyArrangementType::kMarginal> >(&msg));
		tests.push_back(new Seidel_TestGroup< Poly_AG<PolyArrangementType::kSeparate> >(&msg));
		tests.push_back(new Seidel_TestGroup< Poly_AG<PolyArrangementType::kUnbounded> >(&msg));
	}

	TestSetup setup;

	for (TestGroup* test : tests)
		for (uint objType = 0; objType < ObjectiveType::kCount; ++objType)
			if (test->objective_type_supported(objType))
			{
				const string series_name = string(test->name()) + "/" + ObjectiveTypeNames[objType];
				if (find(setup.seriesNames.begin(), setup.seriesNames.end(), series_name) == setup.seriesNames.end()) setup.seriesNames.push_back(series_name);
			}

	setup.D = setup.N = 0;

	map<SKey, StatsF> stats;

	const bool compile = stats_out.stat != StatsOutput::kNone;

	switch (sweep_type)
	{
	case 'N':
		{
		setup.D = 3;
		setup.hAxisVals.assign(autoN, autoN + sizeof(autoN) / sizeof(autoN[0]));
		msg.print(1, "D = %d, sweep ", setup.D);
		msg.print(2, "\n");
		msg.print(1, "N = ");
		for (uint N : setup.hAxisVals)
		{
			msg.print(1, "%d ", N);
			if (!run(stats, setup, data_dir, compile, tests, setup.D, N, M_in, G_in, start_group, nullptr, flush_data, "t_NSWEEP_", msg)) return false;
			if (N != setup.hAxisVals.back()) msg.print(2, "\nN = ");
		}
		msg.print(1, "\n");
		}
		break;
	case 'D':
		{
		setup.N = 1000;
		setup.hAxisVals.assign(autoD, autoD + sizeof(autoD) / sizeof(autoD[0]));
		msg.print(1, "N = %d, sweep ", setup.N);
		msg.print(2, "\n");
		msg.print(1, "D = ");
		for (uint D : setup.hAxisVals)
		{
			msg.print(1, "%d ", D);
			if (!run(stats, setup, data_dir, compile, tests, D, setup.N, M_in, G_in, start_group, nullptr, flush_data, "t_DSWEEP_", msg)) return false;
			if (D != setup.hAxisVals.back()) msg.print(2, "\nD = ");
		}
		msg.print(1, "\n");
		}
		break;
	case '\0':
		{
		msg.print(1, "Running test ID: %s\n", run_type.c_str());
		for (TestGroup* test : tests)
		{
			TestParams params;
			if (test->unpack_test_id(params, run_type.c_str()))
			{
				vector<TestGroup*> single_test = { test };
				setup.N = params.N;
				setup.D = params.D;
				if (!run(stats, setup, data_dir, false, single_test, params.D, params.N, 1, 1, 0, &params, flush_data, "t_", msg)) return false;
				break;
			}
		}
		}
		break;
	}

	if (stats_out.stat != StatsOutput::kNone)
		output_results(setup, stats, stats_out, data_dir, msg);

	// Clean up
	for (TestGroup* test : tests) delete test;

	return true;
}


bool
parse_stats_cmd(StatsOutput& stats_out, const char* cmd)
{
	StatsOutput ret;

	if (cmd == nullptr) return false;

	const char* comma = strchr(cmd, ',');
	if (comma == nullptr)
	{
		cout << "Stat and deviation type must be comma-separated.\n";
		return false;
	}

	if (!strncmp(cmd, "mean", comma - cmd)) ret.stat = StatsOutput::kMean;
	else
	if (!strncmp(cmd, "median", comma - cmd)) ret.stat = StatsOutput::kMedian;
	else
	{
		cout << "Unknown stat: " << string(cmd).substr(0, comma - cmd) << "\n";
		return false;
	}

	if (!strcmp(comma + 1, "rms")) ret.dev = StatsOutput::kRMS;
	else
	if (!strcmp(comma + 1, "mad")) ret.dev = StatsOutput::kMAD;
	else
	{
		cout << "Unknown deviation: " << string(comma + 1) << "\n";
		return false;
	}

	stats_out = ret;
	return true;
}


int
main(int argc, char** argv)
{
	string run_type = "";
	string data_dir = "";
	StatsOutput stats_out;
	bool flush_data = false;
	bool exclusivelyHCP = false;
	uint M = 0;
	uint G = 0;
	uint start_group = 0;
	int verbosity = 1;

	if (argc <= 2)
	{
		cout << "\nUsage:\n\n";
		cout << "perf type dir [c stat,dev] [f] [x] [m M] [g G] [s S] [v verbosity]\n\n";
		cout << "type = 'N', 'D', or a test ID.  If 'N', the total number of half-spaces is\n";
		cout << "    swept through {" << autoN_str << "} while\n";
		cout << "    the number of spatial dimensions is fixed at D = 3.  If 'D', the number of\n";
		cout << "    spatial dimensions is swept through {" << autoD_str << "} while the total\n";
		cout << "    number of half-spaces is fixed at N = 1000.  Otherwise, this argument needs\n";
		cout << "    to be a test ID.  A test ID is output when an error occurs during a run. It\n";
		cout << "    will specify a test, objective, D, N, and arrangement to run.\n";
		cout << "dir = directory to write files containing all measurments.\n";
		cout << "stat,dev: compile statistics from data files.  The stat field must be 'mean' or\n";
		cout << "    'median'.  The 'dev' field must be 'rms' or 'mad' (for standard deviation or\n";
		cout << "    median absolute deviation'.  For example, 'c median,mad'.  If this command\n";
		cout << "    is given then all commands listed below are ignored.\n";
		cout << "f: causes data files to be written to after each measurement.  May impact\n";
		cout << "    performance, but does not affect measurements.\n";
		cout << "x: causes only HCP tests to be performed.\n";
		cout << "M = number of arrangments per test group.  If 0, M is set to (D,N)-dependent\n";
		cout << "    default values.  Default = 0.\n";
		cout << "G = number of test groups for calculating statistics.  If 0, G is set to\n";
		cout << "    (D,N)-dependent default values.  Default = 0.\n";
		cout << "S = starting group #.  Used for test restarts.  Default = 0.\n";
		cout << "o_dir = directory to write statistics files.  If none given, no files will be.\n";
		cout << "    written (stats are written to std out in any case).\n";
		cout << "verbosity = 0, 1, or 2. How much progress text is output.  Default = 1.\n";
		cout << "\n";
		return -1;
	}

	char* const * stop = argv + argc;
	run_type = *++argv;

	data_dir = *++argv;

	while (++argv < stop)
	{
		switch (**argv)
		{
		case 'c': if (++argv < stop && !parse_stats_cmd(stats_out, *argv)) return -1;	break;
		case 'f': flush_data = true;													break;
		case 'x': exclusivelyHCP = true;												break;
		case 'm': if (++argv < stop) M = (uint)(max(0, atoi(*argv)));					break;
		case 'g': if (++argv < stop) G = (uint)(max(0, atoi(*argv)));					break;
		case 's': if (++argv < stop) start_group = (uint)(max(0, atoi(*argv)));			break;
		case 'v': if (++argv < stop) verbosity = atoi(*argv);							break;
		default: cout << "Unknown commandline switch \"" << *argv << "\" ignored.\n";
		}
	}

	MessageFilter msg(verbosity);
	ThreadPriority().set(Priority::Very_High);
	if (!setup_and_run(run_type, data_dir, stats_out, M, G, start_group, exclusivelyHCP, flush_data, msg))
		msg.print(0, "Stopping due to error.\n");
	ThreadPriority().set(Priority::Normal);

	return 0;
}
