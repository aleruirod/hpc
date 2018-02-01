#include <iostream>
#include "ForwardTimeCentralSpace.h"
#include "LaasonenMethod.h"
#include "CrankNicholsonMethod.h"
#include "AnalyticalSolution.h"
#include "Output.h"
#include "Tools.h"
#include <mpi.h>


int main() {

	MPI_Init(NULL, NULL);

	int npes, myrank;
	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// FTCS

	ForwardTimeCentralSpace forwTimeCentSpacSol = ForwardTimeCentralSpace();
	forwTimeCentSpacSol.compute();

	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == npes-1) {

		std::cout << "\n";
		std::cout << "Forward Time Central Space \n";

		Output::printSolution(forwTimeCentSpacSol.getAllSolutions());
		//Output::exportSolution(forwTimeCentSpacSol, 0.1, "output/FTCS_");
	}
	
	// LAASONEN

	LaasonenMethod laasonenSol = LaasonenMethod();
	laasonenSol.compute();

	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == npes-1) {

		std::cout << "\n";
	 	std::cout << "Simple Laasonen\n";

	 	Output::printSolution(laasonenSol.getAllSolutions());
	 	//Output::exportSolution(laasonenSol, 0.1, "output/Simple Laasonen Method_");
	}

	// CRANK NICHOLSON

	CrankNicholsonMethod crankNicholsonSol = CrankNicholsonMethod();
	crankNicholsonSol.compute();

	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == npes-1) {

		std::cout << "\n";
		std::cout << "Crank Nicholson \n";
	
		Output::printSolution(crankNicholsonSol.getAllSolutions());
		//Output::exportSolution(crankNicholsonSol, 0.1, "output/Crank Nicholson_");
	}

	// ANALYTICAL SOLUTION

	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == npes-1) {
	 	std::cout << "\n";
		std::cout << "Analytical Solution \n";

	 	AnalyticalSolution analyt = AnalyticalSolution();
	 	analyt.compute();
	 	Output::printSolution(analyt.getAllSolutions());
	 	//Output::exportSolution(analyt, 0.1, "output/Analytical Solution_");
	}

	MPI_Finalize();

}