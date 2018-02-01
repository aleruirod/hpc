all: ./*.cpp
	mpicxx -std=c++14 -lpthread -W ./*.cpp -o ./bin/HPC
compile:
	make all
	make run
run:
	mpiexec -np 5 ./bin/HPC

