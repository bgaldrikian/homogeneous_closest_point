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

#ifndef _SAMPLE_UTILS_H_
#define _SAMPLE_UTILS_H_

#include <stddef.h>
#include <stdint.h>


// For windows-specific implementations
#undef WINDOWS
#ifdef _MSC_VER
#if defined(_M_IX86) || defined(_M_X64)
#define WINDOWS	1
#endif
#endif
#ifndef WINDOWS
#define WINDOWS 0
#endif


// Timing utility

#if WINDOWS
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include "windows.h"
#else
#include <time.h>
#endif

class Clock
{
public:
#if WINDOWS
						Clock() { QueryPerformanceFrequency((LARGE_INTEGER*)&m_ticksPerSecond); }

	inline long long	ticks() const
						{
							long long result;
							QueryPerformanceCounter((LARGE_INTEGER*)&result);
							return result;
						}
#else	// WINDOWS
						Clock() : m_ticksPerSecond((long long)CLOCKS_PER_SEC) {}

	inline long long	ticks() { return (long long)clock(); }
#endif	// WINDOWS

	inline long long	ticks_per_second() const { return m_ticksPerSecond; }

private:
	long long	m_ticksPerSecond;
};


// Thread priority
struct Priority
{
	enum Value { Very_Low, Low, Normal, High, Very_High };
};

class ThreadPriority
{
public:
	bool set(Priority::Value priority)
	{
#if WINDOWS
		int nPriority = 0;
		switch (priority)
		{
		case Priority::Very_Low:	nPriority = -2;	break;
		case Priority::Low:			nPriority = -1;	break;
		case Priority::Normal:		nPriority = 0;	break;
		case Priority::High:		nPriority = 1;	break;
		case Priority::Very_High:	nPriority = 2;	break;
		default: return false;
		}
		return TRUE == SetThreadPriority(GetCurrentThread(), nPriority);
#endif	// To do: implement for other environments
	}
};


// Memory
template<typename T>
class Buffer
{
public:
	Buffer(size_t size)			{ m_data = size ? new T[size] : nullptr; }
	~Buffer()					{ delete[] m_data; }
	operator T* ()				{ return m_data; }
	operator const T* () const	{ return m_data; }
private:
	T* m_data;
};


template<typename T, size_t Alignment>
class MemBuf
{
public:
	MemBuf(size_t size = 0) : m_mem(nullptr) { resize(size); }
	~MemBuf() { delete[] m_mem; }
	void resize(size_t size)
	{
		size *= sizeof(T);
		delete[] m_mem;
		m_mem = size ? new char[size + Alignment - 1] : nullptr;
		const uintptr_t addr = reinterpret_cast<uintptr_t>(m_mem) + Alignment - 1;
		m_buf = reinterpret_cast<T*>(addr - addr % Alignment);
	}
	operator T*	()				{ return m_buf; }
	operator const T* () const	{ return m_buf; }
private:
	char*	m_mem;
	T*		m_buf;
};


// Misc.
#define ArrayWithString(_type, _name, ...) \
static const _type _name[] = {__VA_ARGS__}; \
static const char _name##_str[] = #__VA_ARGS__


#endif // #ifndef _SAMPLE_UTILS_H_
