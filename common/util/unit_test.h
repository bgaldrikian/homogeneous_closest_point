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

#ifndef _UNIT_TEST_H_
#define _UNIT_TEST_H_


#if defined(_MSC_VER)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#define SetConsoleTextAttribute(_x, _y)	(void)0
#endif

#include <iostream>
#include <vector>
#include <string.h>


/*** Simple framework for unit testing ***/

// To disable a test, make the first character in its name an underscore: '_'.

struct UnitTest
{
	UnitTest(const char* testname) : m_next_test(nullptr)
	{
		if (s_test_list == nullptr) s_test_list = this;
		else s_test_list_end->m_next_test = this;
		s_test_list_end = this;
	}

	bool run()
	{
		std::cout << "*** Starting " << name() << "\n";
		uint32_t error_count = 0;
		exec(error_count);
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | (error_count ? FOREGROUND_RED : FOREGROUND_GREEN));
		std::cout << name() << " completed with " << error_count << " error" << pl_s(error_count) << ".\n\n";
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
		return error_count == 0;
	}

	virtual	const char*	name() const = 0;
	virtual void		exec(uint32_t&) = 0;

	static void run_all(std::vector<const char*> matches = std::vector<const char*>(1, ""))
	{
		uint32_t test_count = 0;
		uint32_t disabled_count = 0;
		for (UnitTest* t = s_test_list; t; t = t->m_next_test) for (const char* match : matches) if (strstr(t->name(), match)) { if (*t->name() != '_') ++test_count; else ++disabled_count; break; }
		std::cout << "\nRunning " << test_count << " test" << pl_s(test_count);
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_GREEN | FOREGROUND_RED);
		if (disabled_count) std::cout << " (" << disabled_count << " disabled)";
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
		std::cout << ".\n\n";
		uint32_t fail_count = 0;
		for (UnitTest* t = s_test_list; t; t = t->m_next_test) for (const char* match : matches) if (strstr(t->name(), match)) { if (*t->name() != '_') fail_count += (uint32_t)!t->run(); break; }
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | (fail_count != 0 ? FOREGROUND_RED : FOREGROUND_GREEN));
		if (!fail_count) std::cout << "All tests PASSED.\n"; else std::cout << fail_count << " test" << pl_s(fail_count) << " FAILED.\n";
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
		std::cout << "\n";
	}

	static const char* pl_s(uint32_t n) { return n == 1 ? "" : "s"; }

	static UnitTest*	s_test_list;
	static UnitTest*	s_test_list_end;

private:
	UnitTest*			m_next_test;
};


#define UNIT_TEST(_testname) \
struct _testname : public UnitTest \
{ \
				_testname()	: UnitTest(#_testname)	{} \
	const char*	name() const override { return #_testname; } \
	void		exec(uint32_t&) override; \
} g_##_testname; \
void _testname::exec(uint32_t& _error_count)


#define EXPECT(_x) \
if (!(_x)) \
{ \
	std::cout << __FILE__ << "(" << __LINE__ << "):\n"; \
	std::cout << "Failed: " << #_x << "\n"; \
	++_error_count; \
} (void)0


#define EXPECT_NEAR(_expected, _x, _tol) \
if (std::abs((_expected) - (_x) > (_tol))) \
{ \
	std::cout << __FILE__ << "(" << __LINE__ << "):\n"; \
	std::cout << "Failed: " << #_x << " (" << _x << ") not within " << #_tol << " (" << _tol << ") of " << #_expected << " (" << _expected << ").\n"; \
	++_error_count; \
} (void)0


#define IMPLIES	== false ||


#define UNIT_TEST_PROGRAM()														\
int																				\
main(int argc, char** argv)														\
{																				\
	if (argc <= 1) UnitTest::run_all();											\
	else UnitTest::run_all(std::vector<const char*>(argv + 1, argv + argc));	\
	return 0;																	\
}																				\
																				\
UnitTest*	UnitTest::s_test_list = nullptr;									\
UnitTest*	UnitTest::s_test_list_end = nullptr


#endif // #ifndef _UNIT_TEST_H_
