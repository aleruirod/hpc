#include "ForwardTimeCentralSpace.h"
#include "Tools.h"

#include <iostream>

// CONSTRUCTORS
/*
* Default constructor (empty solution)
*/

ForwardTimeCentralSpace::ForwardTimeCentralSpace() {}

/*
* Constructor from a vector of double vectors.
*/

ForwardTimeCentralSpace::ForwardTimeCentralSpace(std::vector<std::vector<double>> sols) {
	ForwardTimeCentralSpace::allSolutions = sols;
}

// COMPUTATION
/**
* Computes and stores the values for the Solution using the FTCS method.
*/

void ForwardTimeCentralSpace::compute() {

	std::vector<double> prevSol = Tools::createT0Vector(n+1); // we start by creating a solution vector for the first timestep which will follow the boundary conditions.
	prevSol.reserve(n+1);

	addToAllSolutions(prevSol);// we need to store the values for t = 0.0.

	int npes, myrank;
	MPI_Status status;
	double commTime1 = 0.0, commTime2 = 0.0, compTime1 = 0.0, compTime2 = 0.0, totalCommTime = 0.0, totalCompTime = 0.0;


	if ((n-1) != 1) {
		MPI_Comm_size(MPI_COMM_WORLD, &npes);
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	}else {
		npes = 1;
		myrank = 0;
	}

	int numberPosPerProcess = rint(((double)n - 1)/((double)npes));
	int *counts = new int[npes];
    int *disps  = new int[npes];

	int lastSpaces = 0;

	if (npes != 1 && numberPosPerProcess*npes < (n-1)) {
		if (myrank == npes-1)
			lastSpaces = numberPosPerProcess*npes - (n-1);

	} else if (npes != 1 && numberPosPerProcess*npes > (n-1)) {
		if (myrank == npes-1)
			lastSpaces = 1;

	}

	disps[0] = 0;
	counts[0] = numberPosPerProcess;
	for (int mr = 1; mr < npes; mr++) {
		disps[mr] = disps[mr-1] + counts[mr-1];
		counts[mr] = counts[mr-1];
	}

	counts[npes-1] = numberPosPerProcess - lastSpaces;

	compTime1 = MPI_Wtime();

	for (int j = 1; j < (0.5 / deltaT) + 1; j++) { // the limit for this loop is the number of timesteps needed to get to t = 0.5.

		std::vector<double> compSol = Tools::createT0Vector(numberPosPerProcess - lastSpaces);
		compSol.reserve(numberPosPerProcess - lastSpaces);

		for (int i = 0; i < numberPosPerProcess - lastSpaces; i++) {

			compSol[i] = prevSol[myrank*numberPosPerProcess + i + 1] + ((D*deltaT) / (pow(deltaX, 2))) *
			 (prevSol[myrank*numberPosPerProcess + i + 2] - 2 * prevSol[myrank*numberPosPerProcess + i + 1] +
			 prevSol[myrank*numberPosPerProcess + i]); // this is the formula for the forward time central space scheme.

		}

		commTime1 = MPI_Wtime();
		if ((n-1) != 1)
			MPI_Allgatherv(&compSol[0], numberPosPerProcess - lastSpaces, MPI_DOUBLE, &prevSol[1], counts, disps, MPI_DOUBLE, MPI_COMM_WORLD);
		else if ((n - 1) == 1 && myrank == 0)
			prevSol[1] = compSol[0];

		commTime2 = MPI_Wtime();

		if (myrank == npes-1) {
			if ((j % (int)(0.1 / deltaT)) == 0) // We store the values for t = 0.1, 0.2,...,0.5.
				addToAllSolutions(prevSol);
		}

	}

	compTime2 = MPI_Wtime();

	totalCommTime += commTime2 - commTime1;
	totalCompTime += (compTime2 - compTime1) - totalCommTime;

	std::cout << "\n Communication time for p" << myrank << ": " << totalCommTime << "\n";
	std::cout << "\n Computation time for p" << myrank << ": " << totalCompTime << "\n";


}

