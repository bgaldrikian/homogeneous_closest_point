../bin/unit_tests.exe: unit_tests.cpp tests/discrete_tests.cpp tests/linalg_tests.cpp tests/hcp_tests.cpp tests/hcp_linalg_tests.cpp ../../src/general/hcp.cpp ../../src/general/hcp_toi.cpp
	@mkdir -p ../bin
	g++ -std=c++11 -O2 -o ../bin/unit_tests.exe unit_tests.cpp tests/discrete_tests.cpp tests/linalg_tests.cpp tests/hcp_tests.cpp tests/hcp_linalg_tests.cpp ../../src/general/hcp.cpp ../../src/general/hcp_toi.cpp -I../../../common -I../../../hcp
