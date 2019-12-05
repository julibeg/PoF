make: src/main.cpp
	g++ -o bin/PoFcalc src/main.cpp -lm -fopenmp -Wall -O3
