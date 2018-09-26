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

#ifndef _DISCRETE_H_
#define _DISCRETE_H_


/*** Discrete math utilities ***/


/**
Iterator over all M-element subsets of set {0, 1, ..., N-1}.

Each subset's elements will be ordered.

Example usage:

	int indices[2];
	for (Choose<int> it(4,2,indices); it; ++it) printf("{%d,%d} ", indices[0], indices[1]);

The code above will output:
{0,1} {0,2} {0,3} {1,2} {1,3} {2,3} 
*/
template<typename IndexType>
class Choose
{
public:
	/**
	ctor.

	\param[in] N		The size of the source set.  The source set is {0, 1, ..., N-1}.
	\param[in] M		The size of the choice subsets (must be less than or equal to N).
	\param[out] indices	User-supplied array of M IndexType elements, where the choice subsets will be written.

	Upon return, the first valid subset is written if M <= N.
	*/
	Choose(IndexType N, IndexType M, IndexType* indices) : m_N(N), m_M(M), m_indices(N >= M ? indices : nullptr)
	{
		if (*this) for (IndexType i = 0; i < M; ++i) m_indices[i] = i;
	}

	/**
	Casting this iterator to a bool type tests for end of iteration.
	\return true iff the user-supplied IndexType array is filled with a new valid subset.
	*/
	operator bool() const { return m_indices != nullptr; }

	/**
	Increment operator to generate the next subset.  If all subsets have been written, changes internal
	state such that the bool() cast returns false.
	*/
	void operator ++ ()
	{
		if (!(*this)) return;
		for (IndexType i = 0; i < m_M; ++i)
			if ((i + 1 < m_M && m_indices[i] < m_indices[i + 1] - 1) || (i + 1 == m_M && m_indices[i] + 1 < m_N))
			{
				++m_indices[i];
				for (IndexType j = 0; j < i; ++j) m_indices[j] = j;
				return;
			}
		m_indices = nullptr;
	}

private:
	IndexType	m_N;
	IndexType	m_M;
	IndexType*	m_indices;
};


/**
Iterator over all non-negative integers of type IndexType except those in a given set of indices.

Example usage:

	const int skip_set[4] = { 2, 3, 5, 7 };

	for (Complement<int> i(skip_set, 4); i < 10; ++i) printf("%d ", i);

The code above will output:
0 1 4 6 8 9 
*/
template<typename IndexType>
class Complement
{
public:
	/**
	ctor.

	\param[in] indices		User-supplied IndexType array, the indices to skip.  They MUST BE SORTED in ascending order.
	\param[in] index_count	The size of the indices array.
	*/
	Complement(const IndexType* indices, size_t index_count) : m_current(0), m_indices(indices), m_indices_stop(indices + index_count)
	{
		skip();
	}

	/** Casting this iterator to IndexType gives the current integer in the complement set. */
	operator IndexType() const { return m_current; }

	/** Increment operator to find the next integer in the complement set.  After this the IndexType() cast will return the next integer. */
	IndexType operator ++ ()
	{
		++m_current;
		skip();
		return m_current;
	}

private:
	void	skip()
	{
		while (m_indices < m_indices_stop && *m_indices == m_current)
		{
			++m_indices;
			++m_current;
		}
	}

	IndexType			m_current;
	const IndexType*	m_indices;
	const IndexType*	m_indices_stop;
};


#endif // #ifndef _DISCRETE_H_
