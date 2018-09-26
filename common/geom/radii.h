// Copyright (c) 2018 Bryan Galdrikian
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

#ifndef _RADII_H_
#define _RADII_H_

#include<map>
#include<utility>


// Returns 0 if the value is not found
inline float
find_polytope_radius(unsigned D, unsigned N_B)
{
	static const std::map<std::pair<unsigned, unsigned>, float>
	s_polytope_radii =
	{ //    D    N_B     R
		{ { 1,   50   }, 1.0f      },
		{ { 2,   50   }, 1.0031f   },
		{ { 3,   50   }, 1.057f    },
		{ { 4,   50   }, 1.163f    },
		{ { 5,   50   }, 1.291f    },
		{ { 6,   50   }, 1.430f    },
		{ { 7,   50   }, 1.574f    },
		{ { 8,   50   }, 1.722f    },
		{ { 9,   50   }, 1.874f    },
		{ { 10,  50   }, 2.03f     },
		{ { 3,   5    }, 2.7f      },
		{ { 3,   10   }, 1.40f     },
		{ { 3,   25   }, 1.124f    },
//		{ { 3,   50   }, 1.057f    },
		{ { 3,   100  }, 1.028f    },
		{ { 3,   250  }, 1.0108f   },
		{ { 3,   500  }, 1.0054f   },
		{ { 3,   1000 }, 1.0027f   },
		{ { 3,   2500 }, 1.00107f  },
		{ { 3,   5000 }, 1.00054f  },
		{ { 2,   500  }, 1.000031f },
//		{ { 3,   500  }, 1.0054f   },
		{ { 5,   500  }, 1.072f    },
		{ { 6,   500  }, 1.124f    },
		{ { 7,   500  }, 1.183f    },
		{ { 10,  500  }, 1.370f    },
		{ { 15,  500  }, 1.69f     },
		{ { 20,  500  }, 1.99f     },
		{ { 30,  500  }, 2.56f     },
		{ { 50,  500  }, 3.62f     },
		{ { 60,  500  }, 4.13f     },
		{ { 70,  500  }, 4.61f     },
		{ { 100, 500  }, 6.10f     },
	};

	const auto& it = s_polytope_radii.find(std::make_pair(D, N_B));
	return it != s_polytope_radii.end() ? it->second : 0;
}


#endif // #ifndef _RADII_H_
