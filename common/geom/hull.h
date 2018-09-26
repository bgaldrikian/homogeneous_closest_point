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

#ifndef _HULL_H_
#define _HULL_H_

#include "math/vec_templates/vec_utils.h"
#include "math/vec_templates/geom_utils.h"
#include "math/vec_templates/math_ext.h"
#include "math/discrete.h"

#include <math.h>
#include <vector>
#include <limits>
using namespace std;


/**
	Bounds - axis-aligned bounding box.
*/
template<typename F, int D>
class Bounds
{
public:
					Bounds()							{ setEmpty(); }

	void			setEmpty()							{ extremum[0] = numeric_limits<F>().infinity(); extremum[1] = -numeric_limits<F>().infinity(); }

	bool			isEmpty() const						{ return extremum[0](0) > extremum[1](0); }

	void			include(const Vec<F,D+1>& point)
					{
						for (int i = 0; i < D; ++i)
						{
							const F e = point(i);
							if (e < extremum[0](i)) extremum[0](i) = e;
							if (e > extremum[1](i)) extremum[1](i) = e;
						}
					}

	void			include(const Bounds<F,D>& bounds)	{ include(extremum[0]); include(extremum[1]); }

	const Vec<F,D>&	getMinimum() const					{ return extremum[0]; }

	const Vec<F,D>&	getMaximum() const					{ return extremum[1]; }

	const Vec<F,D>&	operator [] (int index) const		{ return extremum[index&1]; }
	Vec<F,D>&		operator [] (int index)				{ return extremum[index&1]; }

private:
	Vec<F,D>	extremum[2];
};

template<typename F>
struct Extent
{
			Extent()						{ setEmpty(); }

	void	setEmpty()						{ extremum[0] = numeric_limits<F>().infinity(); extremum[1] = -numeric_limits<F>().infinity(); }

	bool	isEmpty() const					{ return extremum[0] > extremum[1]; }

	void	include(F point)				{ extremum[0] = min(extremum[0], point); extremum[1] = max(extremum[1], point); }

	void	include(const Extent& extent)	{ include(extremum[0]); include(extremum[1]); }

	F		getMinimum() const				{ return extremum[0]; }

	F		getMaximum() const				{ return extremum[1]; }

	F	extremum[2];
};


/**
	HullSeparation - output data for Hull classes (below)
*/
template<typename F, int D>
struct HullSeparation
{
	Vec<F,D>	normal;
	Extent<F>	extent0;
	Extent<F>	extent1;
	F			distance;
	int			type;	// 0 = face of hull0, 1 = face of hull1, 2 = edge-edge
};


/**
	HullBase - base for Hull (below)
*/
template<typename F, int D>
class HullBase
{
public:
						HullBase() { setEmpty(); }

	void				setEmpty()
						{
							vertices.resize(0);
							faces.resize(0);
							faceNormalIndices.resize(0);
							widths.resize(0);
							edgeDirections.resize(0);
							bounds.setEmpty();
						}

	bool				isEmpty() const { return vertices.size() == 0; }

	void				transform(const Mat<F,D+1,D+1>& tm);

	// If the distance between the hulls exceeds maxDistance, false is returned.
	// Otherwise, true is returned.  In this case, if 'separation' is not NULL, then separation normal
	// and projected extents are returned in *separation.
	static	bool		hullsInProximity(const HullBase& hull0, const HullBase& hull1, F maxDistance, HullSeparation<F,D>* separation = NULL);

	// Use GJK to determine if the hulls have an intersection
	static	bool		hullsIntersectGJK(const HullBase& hull0, const HullBase& hull1);

	// Returns the min and max dot product of the vertices with the given normal
	void				getExtent(Extent<F>& extent, const Vec<F,D+1>& normal) const
						{
							extent.setEmpty();
							for (unsigned i = 0; i < vertices.size(); ++i)
							{
								extent.include(normal|vertices[i]);
							}
						}

	// Returns the min and max dot product of the vertices in the given face normal direction
	void				getExtent(Extent<F>& extent, int faceNormalNum) const
						{
							extent.extremum[1] = -faces[faceNormalIndices[faceNormalNum]].plane_d();
							extent.extremum[0] = extent.extremum[1] - widths[faceNormalNum];
						}

	// Returns a vertex with the maximum dot product with the given direction.  If there are no vertices, returns the zero vector.
	const Vec<F,D+1>	support(const Vec<F,D+1>& dir) const
						{
							if (vertices.size() == 0) return Vec<F,D+1>((F)0, (F)1);
							F maxDot = -numeric_limits<F>().infinity();
							unsigned supportIndex = 0;
							for (unsigned i = 0; i < vertices.size(); ++i)
							{
								const F dot = vertices[i] | dir;
								if (dot > maxDot)
								{
									maxDot = dot;
									supportIndex = i;
								}
							}
							return vertices[supportIndex];
						}

	// AABB
	const Bounds<F,D>&	getBounds() const	{ return bounds; }

	// Accessors
	unsigned			getVertexCount() const { return vertices.size(); }
	unsigned			getFaceCount() const { return faces.size(); }

protected:
	static bool			hullsInProximity_D(const HullBase<F,D>& hull0, const HullBase<F,D>& hull1, F maxDistance, HullSeparation<F,D>* separation);

	vector< Vec<F,D> >		vertices;
	vector< Vec<F,D+1> >	faces;

	vector< unsigned >		faceNormalIndices;
	vector< F >				widths;			// Same size as faceNormalIndices.  Gives the width of the hull in face normal direction.
	vector< Vec<F,D> >		edgeDirections;
	Bounds<F,D>				bounds;
};

/**
	Hull - precomputed (redundant) information about a convex hull: vertices, hull planes, etc.

	D = 2 and D = 3 supported.
*/
template<typename F, int D>
class Hull : public HullBase<F,D>
{
public:
						Hull() : HullBase<F,D>()	{}

	void				setEmpty() { HullBase<F,D>::setEmpty(); }

	const Vec<F,D+1>&	support(const Vec<F,D+1>& dir) const { return Vec<F,D+1>((F)0, (F)1); }

	void				buildFromHalfspaces(const Vec<F,D+1>* surfaces, unsigned surfaceCount, F boundingSize, F eps);
};

// Specialization for D = 1
template<typename F>
class Hull<F,1> : public HullBase<F,1>
{
public:
					Hull() : HullBase<F,1>()	{}

	void			setEmpty() { HullBase<F,1>::setEmpty(); }

	const Vec<F,2>&	support(const Vec<F,2>& dir) const { return HullBase<F,1>::support(dir); }

	void			buildFromHalfspaces(const Vec<F,2>* surfaces, unsigned surfaceCount, F boundingSize, F eps);
};

// Specialization for D = 2
template<typename F>
class Hull<F,2> : public HullBase<F,2>
{
public:
					Hull() : HullBase<F,2>()	{}

	void			setEmpty() { HullBase<F,2>::setEmpty(); }

	const Vec<F,3>&	support(const Vec<F,3>& dir) const { return HullBase<F,2>::support(dir); }

	void			buildFromHalfspaces(const Vec<F,3>* surfaces, unsigned surfaceCount, F boundingSize, F eps);
};

// Specialization for D = 3
template<typename F>
class Hull<F,3> : public HullBase<F,3>
{
public:
					Hull() : HullBase<F,3>()	{}

	void			setEmpty() { HullBase<F,3>::setEmpty(); }

	const Vec<F,4>&	support(const Vec<F,2>& dir) const { return HullBase<F,3>::support(dir); }

	void			buildFromHalfspaces(const Vec<F,4>* surfaces, unsigned surfaceCount, F boundingSize, F eps);
};


//////////////////////////////////////////////////////
// Function definitions
//////////////////////////////////////////////////////

template<typename F>
inline F
extentDistance(const Extent<F>& extent0, const Extent<F>& extent1)
{
	return max(extent0.getMinimum() - extent1.getMaximum(), extent1.getMinimum() - extent0.getMaximum());
}

template<typename F>
inline F
extentDistanceAndNormalDirection(const Extent<F>& extent0, const Extent<F>& extent1, bool& normalPointsFrom0to1)
{
	const F d01 = extent1.getMinimum() - extent0.getMaximum();
	const F d10 = extent0.getMinimum() - extent1.getMaximum();

	normalPointsFrom0to1 = d01 > d10;

	return normalPointsFrom0to1 ? d01 : d10;
}

template<typename F, int D>
void
HullBase<F, D>::transform(const Mat<F, D + 1, D + 1>& tm)
{
	const Mat<F, D + 1, D + 1> cofTM = cof(tm);

	// Faces
	const unsigned faceCount = (unsigned)faces.size();
	for (unsigned i = 0; i < faceCount; ++i)
	{
		faces[i] = cofTM*faces[i];
	}

	// Face normals and widths
	const unsigned faceNormalCount = (unsigned)faceNormalIndices.size();
	for (unsigned i = 0; i < faceNormalCount; ++i)
	{
		widths[i] /= faces[faceNormalIndices[i]].plane_n().length();
	}

	// Now we can normalize the faces
	for (unsigned i = 0; i < faceCount; ++i)
	{
		faces[i].plane_normalize();
	}

	// Vertices and bounds
	bounds.setEmpty();
	const unsigned vertexCount = (unsigned)vertices.size();
	for (unsigned i = 0; i < vertexCount; ++i)
	{
		vertices[i] = tm*vertices[i];
		bounds.include(vertices[i]);
	}

	// Edge directions
	const unsigned edgeCount = (unsigned)edgeDirections.size();
	for (unsigned i = 0; i < edgeCount; ++i)
	{
		edgeDirections[i] = tm*edgeDirections[i];
	}
}

// Templated "3D" tests test cases only found in 3 (or greater) dimensions.  The default implementation will return true, as no test needs to be performed.
template<typename F, int D>
bool
hulls3DTests(F&, const HullBase<F, D>&, const vector< Vec<F, D> >&, const HullBase<F, D>&, const vector< Vec<F, D> >&, F, HullSeparation<F, D>*)
{
	return true;
}

template<typename F>
bool
hulls3DTests(F& result, const HullBase<F, 3>& hull0, const vector< Vec<F, 3> >& edgeDirections0, const HullBase<F, 3>& hull1, const vector< Vec<F, 3> >& edgeDirections1, F maxDistance, HullSeparation<F, 3>* separation)
{
	// Test hulls against cross-edge planes
	const unsigned edgeCount0 = (unsigned)edgeDirections0.size();
	const unsigned edgeCount1 = (unsigned)edgeDirections1.size();
	for (unsigned edge0Index = 0; edge0Index < edgeCount0; ++edge0Index)
	{
		const Vec<F, 3>& edge0 = edgeDirections0[edge0Index];
		for (unsigned edge1Index = 0; edge1Index < edgeCount1; ++edge1Index)
		{
			const Vec<F, 3>& edge1 = edgeDirections1[edge1Index];
			Vec<F, 3> n = edge0^edge1;
			const F n2 = n.length_squared();
			if (n2 < std::numeric_limits<F>().epsilon()*std::numeric_limits<F>().epsilon())
			{
				continue;
			}
			n *= (F)1 / sqrt(n2);
			Extent<F> extent0, extent1;
			hull0.getExtent(extent0, n);
			hull1.getExtent(extent1, n);
			bool normalPointsFrom0to1;
			const F dist = extentDistanceAndNormalDirection(extent0, extent1, normalPointsFrom0to1);
			if (dist > result)
			{
				if (dist > maxDistance)
				{
					return false;
				}
				result = dist;
				if (separation != NULL)
				{
					if (normalPointsFrom0to1)
					{
						separation->normal = n;
						separation->extent0 = extent0;
						separation->extent1 = extent1;
					}
					else
					{
						separation->normal = -n;
						separation->extent0.extremum[0] = -extent0.extremum[1];
						separation->extent0.extremum[1] = -extent0.extremum[0];
						separation->extent1.extremum[0] = -extent1.extremum[1];
						separation->extent1.extremum[1] = -extent1.extremum[0];
					}
					separation->distance = dist;
					separation->type = 2;
				}
			}
		}
	}

	return true;
}


template<typename F, int D>
bool
findProjectedOverlap(F& result, const Vec<F, D>& n, const HullBase<F, D>& hull0, const HullBase<F, D>& hull1, F maxDistance, HullSeparation<F, D>* separation)	// n must be normalized
{
	Extent<F> extent0, extent1;
	hull0.getExtent(extent0, n);
	hull1.getExtent(extent1, n);
	bool normalPointsFrom0to1;
	const F dist = extentDistanceAndNormalDirection(extent0, extent1, normalPointsFrom0to1);
	if (dist > result)
	{
		if (dist > maxDistance)
		{
			return false;
		}
		result = dist;
		if (separation != NULL)
		{
			if (normalPointsFrom0to1)
			{
				separation->normal = n;
				separation->extent0 = extent0;
				separation->extent1 = extent1;
			}
			else
			{
				separation->normal = -n;
				separation->extent0.extremum[0] = -extent0.extremum[1];
				separation->extent0.extremum[1] = -extent0.extremum[0];
				separation->extent1.extremum[0] = -extent1.extremum[1];
				separation->extent1.extremum[1] = -extent1.extremum[0];
			}
			separation->distance = dist;
		}
	}

	return true;
}


template<typename F, int D>
bool
HullBase<F, D>::hullsInProximity_D(const HullBase<F, D>& hull0, const HullBase<F, D>& hull1, F maxDistance, HullSeparation<F, D>* separation)
{
	F result = -std::numeric_limits<F>().infinity();

	int indices[D + 1];
	Mat<F, D + 1, D + 1> S;

	// Find exhaustive set of normals n, and use findProjectedOverlap on them all
	for (int N0 = 1; N0 <= D; ++N0)	// Number of faces to use from hull0
	{
		const int N1 = D + 1 - N0;	// Number of faces to use from hull1
		for (Choose<int> choose0((int)hull0.getFaceCount(), N0, indices); choose0; ++choose0)
		{
			for (int i = 0; i < N0; ++i) S(i) = hull0.faces[indices[i]];
			for (Choose<int> choose1((int)hull1.getFaceCount(), N1, indices + N0); choose1; ++choose1)
			{
				for (int i = N0; i < D + 1; ++i) S(i) = hull1.faces[indices[i]];
				Vec<F, D> n((F)0);
				for (int i = 0; i < N0; ++i) n += cof(S, D, i)*S(i).plane_n();
				if (n.vector_normalize() >= std::numeric_limits<F>().epsilon())
				{
					if (!findProjectedOverlap(result, n, hull0, hull1, maxDistance, separation)) return false;
				}
			}
		}
	}

	return true;
}


template<typename F, int D>
bool
HullBase<F, D>::hullsInProximity(const HullBase& hull0, const HullBase& hull1, F maxDistance, HullSeparation<F, D>* separation)
{
	if (D > 3) return hullsInProximity_D(hull0, hull1, maxDistance, separation);

	F result = -std::numeric_limits<F>().infinity();

	// Test hull1 against faces of hull0
	const unsigned faceCount0 = (unsigned)hull0.faceNormalIndices.size();
	for (unsigned face0Num = 0; face0Num < faceCount0; ++face0Num)
	{
		const Vec<F, D>& normal = hull0.faces[hull0.faceNormalIndices[face0Num]].plane_n();
		Extent<F> extent0;
		hull0.getExtent(extent0, face0Num);
		Extent<F> extent1;
		hull1.getExtent(extent1, normal);
		const F dist = extentDistance(extent1, extent0);
		if (dist > result)
		{
			if (dist > maxDistance)
			{
				return false;
			}
			result = dist;
			if (separation != NULL)
			{
				separation->normal = normal;
				separation->extent0 = extent0;
				separation->extent1 = extent1;
				separation->distance = dist;
				separation->type = 0;
			}
		}
	}

	if (D == 1)
	{
		return true;	// Only need to perform one test in 1D, as there's only one direction to test
	}

	// Test hull0 against faces of hull1
	const unsigned faceCount1 = (unsigned)hull1.faceNormalIndices.size();
	for (unsigned face1Num = 0; face1Num < faceCount1; ++face1Num)
	{
		const Vec<F, D>& normal = hull1.faces[hull1.faceNormalIndices[face1Num]].plane_n();
		Extent<F> extent0;
		hull0.getExtent(extent0, normal);
		Extent<F> extent1;
		hull1.getExtent(extent1, face1Num);
		const F dist = extentDistance(extent0, extent1);
		if (dist > result)
		{
			if (dist > maxDistance)
			{
				return false;
			}
			result = dist;
			if (separation != NULL)
			{
				separation->normal = normal;
				separation->extent0 = extent0;
				separation->extent1 = extent1;
				separation->distance = dist;
				separation->type = 1;
			}
		}
	}

	return hulls3DTests(result, hull0, hull0.edgeDirections, hull1, hull1.edgeDirections, maxDistance, separation);
}


// GJK utility class
template<typename F, int D>
struct GJK
{
	bool update(Vec<F, D> simplex[D + 1], unsigned&, Vec<F, D>&) { (void)simplex; return true; }
};

// update GJK state 1D
template<typename F>
struct GJK<F, 1>
{
	bool
		update(Vec<F, 1> simplex[2], unsigned& simplexSize, Vec<F, 1>& dir)
	{
		dir = simplex[1](0) - simplex[0](0);
		if ((simplex[1] | dir) <= (F)0) return true;
		simplex[0] = simplex[--simplexSize];	// simplexSize should have been 2
		return false;
	}
};

// update GJK state 2D
template<typename F>
struct GJK<F, 2>
{
	bool
		update(Vec<F, 2> simplex[3], unsigned& simplexSize, Vec<F, 2>& dir)
	{
		switch (simplexSize)
		{
		case 2:
		{
			dir = ~Vec<F, 2>(simplex[1] - simplex[0]);
			if ((dir | simplex[0]) > (F)0) dir = -dir;
		} break;
		case 3:
		{
			F greatestOriginArea = -std::numeric_limits<F>().infinity();
			unsigned vertexToRemove = 2;
			for (unsigned i = 0; i < 2; ++i)
			{
				Vec<F, 2> testDir = ~Vec<F, 2>(simplex[i ^ 1] - simplex[2]);
				if ((testDir | (simplex[i] - simplex[2])) >(F)0) testDir = -testDir;
				const F originArea = -(simplex[2] | testDir);
				if (originArea > greatestOriginArea)
				{
					greatestOriginArea = originArea;
					dir = testDir;
					vertexToRemove = i;
				}
			}
			if (greatestOriginArea <= (F)0) return true;
			simplex[vertexToRemove] = simplex[--simplexSize];
		} break;
		}
		return false;
	}
};

// update GJK state 3D
template<typename F>
struct GJK<F, 3>
{
	bool
		update(Vec<F, 3> simplex[4], unsigned& simplexSize, Vec<F, 3>& dir)
	{
		switch (simplexSize)
		{
		case 2:
		{
			const Vec<F, 3> line = simplex[1] - simplex[0];
			dir = line ^ (line^simplex[0]);
		} break;
		case 3:
		{
			dir = Vec<F, 3>(simplex[1] - simplex[0]) ^ Vec<F, 3>(simplex[2] - simplex[0]);
			if ((dir | simplex[0]) > (F)0) dir = -dir;
		} break;
		case 4:
		{
			F greatestOriginDistance = -std::numeric_limits<F>().infinity();
			unsigned vertexToRemove = 3;
			for (unsigned i = 0; i < 3; ++i)
			{
				Vec<F, 3> testDir = Vec<F, 3>(simplex[(i + 1) % 3] - simplex[3]) ^ Vec<F, 3>(simplex[(i + 2) % 3] - simplex[3]);
				testDir.vector_normalize();
				if ((testDir | (simplex[i] - simplex[3])) >(F)0) testDir = -testDir;
				const F originDistance = -(simplex[3] | testDir);
				if (originDistance > greatestOriginDistance)
				{
					greatestOriginDistance = originDistance;
					dir = testDir;
					vertexToRemove = i;
				}
			}
			if (greatestOriginDistance <= (F)0) return true;
			simplex[vertexToRemove] = simplex[--simplexSize];
		} break;
		}
		return false;
	}
};


template<typename F, int D>
bool
HullBase<F, D>::hullsIntersectGJK(const HullBase& hull0, const HullBase& hull1)
{
	Vec<F, D> simplex[D + 1];
	unsigned simplexSize = 0;
	Vec<F, D> dir = hull0.bounds[0] + hull0.bounds[1] - hull1.bounds[0] - hull1.bounds[1];
	if (dir.length_squared() == (F)0) dir(0) = (F)1;
	simplex[simplexSize++] = hull0.support(dir) - hull1.support(-dir);
	dir = -dir;
	for (;;)
	{
		simplex[simplexSize++] = hull0.support(dir) - hull1.support(-dir);
		if ((simplex[simplexSize - 1] | dir) < (F)0) return false;
		if (GJK<F, D>().update(simplex, simplexSize, dir)) return true;
	}
}


template<typename F, int D>
void
Hull<F, D>::buildFromHalfspaces(const Vec<F, D + 1>* surfaces, unsigned surfaceCount, F boundingSize, F eps)
{
	// For general D, only compute vertices
	setEmpty();

	// Reserve faces for the bounding box and the input surfaces
	HullBase<F, D>::faces.reserve(D * 2 + surfaceCount);

	// Fill faces array
	for (unsigned faceN = 0; faceN < D * 2; ++faceN)
	{
		Vec<F, D> normal((F)0);
		normal(faceN >> 1) = (faceN & 1) ? (F)1 : -(F)1;
		HullBase<F, D>::faces.push_back(Vec<F, D + 1>(normal, -boundingSize));
	}
	for (unsigned surfaceIndex = 0; surfaceIndex < surfaceCount; ++surfaceIndex)
	{
		HullBase<F, D>::faces.push_back(surfaces[surfaceIndex]);
	}

	unsigned faceCount = HullBase<F, D>::faces.size();
	HullBase<F, D>::vertices.reserve(faceCount);
	vector<bool> usePlane(faceCount, false);

	// Iterate through all (N chooose D) combinations of planes
	unsigned indices[D];
	Vec<F, D + 1> v[D];
	for (int i = 0; i < D; ++i)
	{
		indices[i] = i;
		v[i] = HullBase<F, D>::faces[i];
	}
	bool done = false;
	do
	{
		Vec<F, D + 1> pos;
		if (perp_D_safe<F, D>(pos, v))
		{
			pos.point_normalize();
			bool keepVertex = true;
			for (unsigned i = 0; keepVertex && i < faceCount; ++i)
			{
				bool ignore = false;
				for (unsigned j = 0; !ignore && j < D; ++j) ignore = (i == indices[j]);
				if (ignore) continue;
				keepVertex = (HullBase<F, D>::faces[i] | pos) < eps;
			}
			if (keepVertex)
			{
				HullBase<F, D>::vertices.push_back(pos);
				HullBase<F, D>::bounds.include(pos);
			}
		}
		done = true;
		for (int i = D; i--;)
		{
			if (++indices[i] < faceCount - (D - i - 1))
			{
				v[i] = HullBase<F, D>::faces[indices[i]];
				for (int j = i + 1; j < D; ++j)
				{
					indices[j] = indices[j - 1] + 1;
					v[j] = HullBase<F, D>::faces[indices[j]];
				}
				done = false;
				break;
			}
		}
	} while (!done);
}

// Specialization for D = 1
template<typename F>
void
Hull<F, 1>::buildFromHalfspaces(const Vec<F, 2>* surfaces, unsigned surfaceCount, F boundingSize, F /*eps*/)
{
	setEmpty();
	for (unsigned surfaceN = 0; surfaceN < surfaceCount; ++surfaceN)
	{
		HullBase<F, 1>::bounds.include(project<F, 1>((F)0, surfaces[surfaceN]));
	}

	if (HullBase<F, 1>::bounds.isEmpty())
	{
		return;
	}

	HullBase<F, 1>::vertices.resize(2);
	HullBase<F, 1>::faces.resize(2);
	HullBase<F, 1>::faceNormalIndices.resize(1);
	HullBase<F, 1>::widths.resize(1);
	for (int i = 0; i < 2; ++i)
	{
		HullBase<F, 1>::vertices[i] = HullBase<F, 1>::bounds[i];
		HullBase<F, 1>::faces[i] = Vec<F, 2>((i ? (F)1 : -(F)1), (i ? -(F)1 : (F)1)*HullBase<F, 1>::bounds[i](0));
	}
	HullBase<F, 1>::faceNormalIndices[0] = 0;
	HullBase<F, 1>::widths[0] = HullBase<F, 1>::vertices[1](0) - HullBase<F, 1>::vertices[0](0);
}

// Specialization for D = 2
template<typename F>
void
Hull<F, 2>::buildFromHalfspaces(const Vec<F, 3>* surfaces, unsigned surfaceCount, F boundingSize, F eps)
{
	const F eps2 = eps*eps;

	setEmpty();

	// Reserve faces for the bounding box and the input surfaces
	HullBase<F, 2>::faces.reserve(4 + surfaceCount);

	// Fill faces array, but redundant ones will be removed
	for (unsigned faceN = 0; faceN < 4; ++faceN)
	{
		Vec<F, 2> normal((F)0);
		normal(faceN >> 1) = (faceN & 1) ? (F)1 : -(F)1;
		HullBase<F, 2>::faces.push_back(Vec<F, 3>(normal, -boundingSize));
	}
	for (unsigned surfaceIndex = 0; surfaceIndex < surfaceCount; ++surfaceIndex)
	{
		HullBase<F, 2>::faces.push_back(surfaces[surfaceIndex]);
	}

	unsigned faceCount = HullBase<F, 2>::faces.size();
	HullBase<F, 2>::vertices.reserve(faceCount);
	for (unsigned i = faceCount; i--;)
	{
		Vec<F, 3>& planeI = HullBase<F, 2>::faces[i];
		// Create line from the intersection of planeI
		const Vec<F, 3> pos(project<F, 2>((F)0, planeI), (F)1);
		const Vec<F, 3> dir(~Vec<F, 2>(planeI.plane_n()), (F)0);
		// Intersect line against all other planes (j != i) to find edge
		F maxS = std::numeric_limits<F>().infinity();
		F minS = -std::numeric_limits<F>().infinity();
		bool edgeFound = true;	// Until proven otherwise
		for (unsigned j = 0; j < HullBase<F, 2>::faces.size(); ++j)
		{
			if (j == i)
			{
				continue;
			}
			const Vec<F, 3>& planeJ = HullBase<F, 2>::faces[j];
			const F num = -(planeJ | pos);
			const F den = (planeJ | dir);
			if (fabs(den) > std::numeric_limits<F>().epsilon())
			{
				const F s = num / den;
				if (den > (F)0)
				{
					maxS = min(maxS, s);
				}
				else
				{
					minS = max(minS, s);
				}
				if (minS > maxS)
				{
					edgeFound = false;
					break;
				}
			}
			else
				if (num < -eps)
				{
					edgeFound = false;	// Outside of planeJ
					break;
				}
		}
		if (edgeFound)
		{
			HullBase<F, 2>::vertices.push_back(pos + minS*dir);
			HullBase<F, 2>::vertices.push_back(pos + maxS*dir);
		}
		else
		{
			planeI = HullBase<F, 2>::faces[--faceCount];	// planeI was redundant
		}
	}
	HullBase<F, 2>::faces.resize(faceCount);

	if (HullBase<F, 2>::vertices.size() == 0 || HullBase<F, 2>::faces.size() == 0)
	{
		setEmpty();
		return;
	}

	// Cull vertices and create bounds
	unsigned vertexCount = HullBase<F, 2>::vertices.size();
	for (unsigned i = 0; i < vertexCount; ++i)
	{
		const Vec<F, 2>& vertex = HullBase<F, 2>::vertices[i];
		HullBase<F, 2>::bounds.include(vertex);
		for (int j = vertexCount; --j >(int)i;)
		{
			if ((HullBase<F, 2>::vertices[j] - vertex).length_squared() <= eps2)
			{
				swap(HullBase<F, 2>::vertices[j], HullBase<F, 2>::vertices[--vertexCount]);
				HullBase<F, 2>::vertices.pop_back();
			}
		}
	}

	// Now reduce the set of face normals
	HullBase<F, 2>::faceNormalIndices.reserve(faceCount);
	for (unsigned i = 0; i < faceCount; ++i)
	{
		const Vec<F, 2> testNormal = HullBase<F, 2>::faces[i].plane_n();
		bool normalFound = false;
		for (unsigned j = 0; j < HullBase<F, 2>::faceNormalIndices.size(); ++j)
		{
			if (fabs(testNormal | ~Vec<F, 2>(HullBase<F, 2>::faces[HullBase<F, 2>::faceNormalIndices[j]].plane_n())) < eps)
			{
				normalFound = true;
				break;
			}
		}
		if (!normalFound)
		{
			HullBase<F, 2>::faceNormalIndices.push_back(i);
		}
	}

	// Now create widths
	HullBase<F, 2>::widths.resize(HullBase<F, 2>::faceNormalIndices.size(), (F)0);
	for (unsigned i = 0; i < HullBase<F, 2>::faceNormalIndices.size(); ++i)
	{
		for (unsigned j = 0; j < HullBase<F, 2>::vertices.size(); ++j)
		{
			HullBase<F, 2>::widths[i] = max(HullBase<F, 2>::widths[i], -(HullBase<F, 2>::faces[HullBase<F, 2>::faceNormalIndices[i]] | Vec<F, 3>(HullBase<F, 2>::vertices[j], (F)1)));
		}
	}
}

// Specialization for D = 3
template<typename F>
void
Hull<F, 3>::buildFromHalfspaces(const Vec<F, 4>* surfaces, unsigned surfaceCount, F boundingSize, F eps)
{
	const F eps2 = eps*eps;

	setEmpty();

	// Reserve faces for the bounding box and the input surfaces
	HullBase<F, 3>::faces.reserve(6 + surfaceCount);

	// Fill faces array, but redundant ones will be removed
	for (unsigned faceN = 0; faceN < 6; ++faceN)
	{
		Vec<F, 3> normal((F)0);
		normal(faceN >> 1) = (faceN & 1) ? (F)1 : -(F)1;
		HullBase<F, 3>::faces.push_back(Vec<F, 4>(normal, -boundingSize));
	}
	for (unsigned surfaceIndex = 0; surfaceIndex < surfaceCount; ++surfaceIndex)
	{
		HullBase<F, 3>::faces.push_back(surfaces[surfaceIndex]);
	}

	size_t faceCount = HullBase<F, 3>::faces.size();
	HullBase<F, 3>::vertices.reserve(faceCount);
	vector<bool> usePlane(faceCount, false);
	for (unsigned i = 0; i < faceCount; ++i)
	{
		Vec<F, 4>& planeI = HullBase<F, 3>::faces[i];
		for (unsigned j = i + 1; j < faceCount; ++j)
		{
			const Vec<F, 4>& planeJ = HullBase<F, 3>::faces[j];
			// Create line from the intersection of planeI and planeJ
			Vec<F, 3> pos;
			Vec<F, 3> dir;
			if (!intersect_planes(pos, dir, planeI, planeJ))
			{
				continue;
			}
			const Vec<F, 4> posW(pos, (F)1);
			const Vec<F, 4> dirW(dir, (F)0);
			// Intersect line against all other planes (k != i && k != j) to find edge
			F maxS = std::numeric_limits<F>().infinity();
			F minS = -std::numeric_limits<F>().infinity();
			bool edgeFound = true;	// Until proven otherwise
			for (unsigned k = 0; k < faceCount; ++k)
			{
				if (k == i || k == j)
				{
					continue;
				}
				const Vec<F, 4>& planeK = HullBase<F, 3>::faces[k];
				const F num = -(planeK | posW);
				const F den = (planeK | dirW);
				if (fabs(den) > std::numeric_limits<F>().epsilon())
				{
					const F s = num / den;
					if (den > (F)0)
					{
						if (s < maxS) maxS = s;
					}
					else
					{
						if (s > minS) minS = s;
					}
					if (minS > maxS)
					{
						edgeFound = false;
						break;
					}
				}
				else
					if (num < -eps)
					{
						edgeFound = false;	// Outside of planeK
						break;
					}
			}
			if (edgeFound)
			{
				HullBase<F, 3>::vertices.push_back(pos + minS*dir);
				HullBase<F, 3>::vertices.push_back(pos + maxS*dir);
				HullBase<F, 3>::edgeDirections.push_back(dir);
				usePlane[i] = true;
				usePlane[j] = true;
			}
		}
	}

	// Pack faces with just ones that are relevant
	faceCount = 0;
	for (unsigned i = 0; i < HullBase<F, 3>::faces.size(); ++i)
	{
		if (usePlane[i])
		{
			HullBase<F, 3>::faces[faceCount++] = HullBase<F, 3>::faces[i];
		}
	}
	HullBase<F, 3>::faces.resize(faceCount);

	if (HullBase<F, 3>::vertices.size() == 0 || HullBase<F, 3>::faces.size() == 0)
	{
		setEmpty();
		return;
	}

	// Cull vertices and create bounds
	size_t vertexCount = HullBase<F, 3>::vertices.size();
	for (size_t i = 0; i < vertexCount; ++i)
	{
		const Vec<F, 3>& vertex = HullBase<F, 3>::vertices[i];
		HullBase<F, 3>::bounds.include(vertex);
		for (size_t j = vertexCount; j-- > i + 1;)
		{
			if ((HullBase<F, 3>::vertices[j] - vertex).length_squared() <= 1000 * eps2)
			{
				swap(HullBase<F, 3>::vertices[j], HullBase<F, 3>::vertices[--vertexCount]);
				HullBase<F, 3>::vertices.pop_back();
			}
		}
	}

	// Now reduce the set of face normals
	HullBase<F, 3>::faceNormalIndices.reserve(faceCount);
	for (unsigned i = 0; i < faceCount; ++i)
	{
		const Vec<F, 3> testNormal = HullBase<F, 3>::faces[i].plane_n();
		bool normalFound = false;
		for (unsigned j = 0; j < HullBase<F, 3>::faceNormalIndices.size(); ++j)
		{
			if ((testNormal^Vec<F, 3>(HullBase<F, 3>::faces[HullBase<F, 3>::faceNormalIndices[j]]).plane_n()).length_squared() < eps2)
			{
				normalFound = true;
				break;
			}
		}
		if (!normalFound)
		{
			HullBase<F, 3>::faceNormalIndices.push_back(i);
		}
	}

	// Now reduce the set of edge directions
	unsigned edgeDirectionCount = 0;
	for (unsigned i = 0; i < HullBase<F, 3>::edgeDirections.size(); ++i)
	{
		const Vec<F, 3>& testDir = HullBase<F, 3>::edgeDirections[i];
		bool dirFound = false;
		for (unsigned j = 0; j < edgeDirectionCount; ++j)
		{
			if ((testDir^HullBase<F, 3>::edgeDirections[j]).length_squared() < eps2)
			{
				dirFound = true;
				break;
			}
		}
		if (!dirFound)
		{
			HullBase<F, 3>::edgeDirections[edgeDirectionCount++] = testDir;
		}
	}
	HullBase<F, 3>::edgeDirections.resize(edgeDirectionCount);

	// Now create widths
	HullBase<F, 3>::widths.resize(HullBase<F, 3>::faceNormalIndices.size(), (F)0);
	for (unsigned i = 0; i < HullBase<F, 3>::faceNormalIndices.size(); ++i)
	{
		for (unsigned j = 0; j < HullBase<F, 3>::vertices.size(); ++j)
		{
			const F width = -(HullBase<F, 3>::faces[HullBase<F, 3>::faceNormalIndices[i]] | Vec<F, 4>(HullBase<F, 3>::vertices[j], (F)1));
			if (width > HullBase<F, 3>::widths[i]) HullBase<F, 3>::widths[i] = width;
		}
	}
}


#endif // #ifndef _HULL_H_
