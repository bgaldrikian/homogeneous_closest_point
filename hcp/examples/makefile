all: ./bin/ex1.exe ./bin/ex2.exe ./bin/ex3.exe ./bin/ex4.exe ./bin/ex5.exe ./bin/ex6.exe ./bin/ex7.exe ./bin/ex8.exe

./bin/ex1.exe: ex1.cpp ../src/general/hcp.cpp
	@mkdir -p ./bin
	g++ -std=c++14 -o ./bin/ex1.exe ex1.cpp ../src/general/hcp.cpp -I../../common -I../../hcp

./bin/ex2.exe: ex2.cpp ../src/general/hcp.cpp
	@mkdir -p ./bin
	g++ -std=c++14 -o ./bin/ex2.exe ex2.cpp ../src/general/hcp.cpp -I../../common -I../../hcp

./bin/ex3.exe: ex3.cpp ../src/general/hcp.cpp
	@mkdir -p ./bin
	g++ -std=c++14 -o ./bin/ex3.exe ex3.cpp ../src/general/hcp.cpp -I../../common -I../../hcp

./bin/ex4.exe: ex4.cpp ../src/general/hcp.cpp
	@mkdir -p ./bin
	g++ -std=c++14 -o ./bin/ex4.exe ex4.cpp ../src/general/hcp.cpp -I../../common -I../../hcp

./bin/ex5.exe: ex5.cpp ../src/general/hcp.cpp
	@mkdir -p ./bin
	g++ -std=c++14 -o ./bin/ex5.exe ex5.cpp ../src/general/hcp.cpp -I../../common -I../../hcp

./bin/ex6.exe: ex6.cpp ../src/general/hcp.cpp ../src/general/hcp_toi.cpp
	@mkdir -p ./bin
	g++ -std=c++14 -o ./bin/ex6.exe ex6.cpp ../src/general/hcp.cpp ../src/general/hcp_toi.cpp -I../../common -I../../hcp

./bin/ex7.exe: ex7.cpp
	@mkdir -p ./bin
	g++ -std=c++14 -o ./bin/ex7.exe ex7.cpp -I../../common -I../../hcp

./bin/ex8.exe: ex8.cpp
	@mkdir -p ./bin
	g++ -std=c++14 -o ./bin/ex8.exe ex8.cpp -I../../common -I../../hcp
