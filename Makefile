all:
	g++ -c -o vector.o vector.cpp -std=c++23
	g++ -c -o global_vars.o global_variables.cpp -std=c++23
	g++ -c -o main.o main.cpp -std=c++23
	g++ -o program.out main.o vector.o global_vars.o -lraylib -std=c++23
