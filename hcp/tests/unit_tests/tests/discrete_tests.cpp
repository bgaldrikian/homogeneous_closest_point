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


#include "util/unit_test.h"

#include "math/discrete.h"

#include <vector>
#include <algorithm>
using namespace std;


static inline unsigned factorial(unsigned N) { return N ? N*factorial(N - 1) : 1; }


UNIT_TEST(Choose_Iterator)
{
	const unsigned testMax = 10;
	vector<unsigned> indices(testMax);
	for (unsigned N = 0; N <= testMax; ++N)
	{
		for (unsigned M = 0; M <= testMax; ++M)
		{
			unsigned count = 0;
			for (Choose<unsigned> c(N, M, indices.data()); c; ++c) ++count;
			EXPECT(M > N IMPLIES count == 0);
			EXPECT(M <= N IMPLIES count == factorial(N) / (factorial(N - M)*factorial(M)));
		}
	}
}


UNIT_TEST(Complement_Iterator)
{
	const unsigned iterMax = 5;
	const unsigned indicesMax = 10;
	vector<unsigned> indices(indicesMax);
	uint32_t covered = 0;

	for (unsigned indicesCount = 0; indicesCount < indicesMax / 2; ++indicesCount)
	{
		for (Choose<unsigned> c(indicesMax, indicesCount, indices.data()); c; ++c)
		{
			for (unsigned j = 0; j < iterMax; ++j)
			{
				if (find(indices.begin(), indices.end(), j) != indices.end()) covered |= (1 << j);
			}

			for (Complement<unsigned> j(indices.data(), indicesCount); j < iterMax; ++j)
			{
				EXPECT(j < iterMax);
				EXPECT(find(indices.begin(), indices.begin() + indicesCount, j) == indices.begin() + indicesCount);
				covered |= (1 << j);
			}
		}

		EXPECT(covered == (1 << iterMax) - 1);
	}
}
