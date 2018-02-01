#include "CrankNicholsonMethod.h"
#include "Tools.h"


// CONSTRUCTORS
/*
* Default constructor (empty solution)
*/

CrankNicholsonMethod::CrankNicholsonMethod() {}

/*
* Constructor from a vector of double vectors.
*/

CrankNicholsonMethod::CrankNicholsonMethod(std::vector<std::vector<double>> sols) {
	CrankNicholsonMethod::allSolutions = sols;
}

// COMPUTATION
/**
* Computes and stores the values for the Solution using the Crank Nicholson method, for this, we will carry out the Thomas Algorithm.
*/

void CrankNicholsonMethod::compute() {

	std::vector<double> previousX = Tools::createT0Vector(n + 1);
	previousX.reserve(n + 1); // we need to keep track for the solutions of the previous timestep to satisfy the Thomas algorithm.

	addToAllSolutions(previousX);

	int npes, myrank;
	MPI_Status status;
	double compTime1 = 0.0, compTime2 = 0.0, totalCommTime = 0.0, totalCompTime = 0.0;
	double commTime1 = 0.0, commTime2 = 0.0, commTime3 = 0.0, commTime4 = 0.0, commTime5 = 0.0, commTime6 = 0.0;
	double commTime7 = 0.0, commTime8 = 0.0, commTime9 = 0.0, commTime10 = 0.0;

	if ((n-1) != 1) {
		MPI_Comm_size(MPI_COMM_WORLD, &npes);
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	} else {
		npes = 1;
		myrank = 0;
	}
	
	int numberPosPerProcess = rint(((double)n + 1)/((double)npes));
	int *counts = new int[npes];
    int *disps  = new int[npes];

	int lastSpaces = 0;

	if (npes != 1 && numberPosPerProcess*npes < (n+1)) {
		if (myrank == npes-1)
			lastSpaces = numberPosPerProcess*npes - (n+1);

	} else if (npes != 1 && numberPosPerProcess*npes > (n+1)) {
		if (myrank == npes-1)
			lastSpaces = 1;
	}

	disps[0] = 0;
	counts[0] = numberPosPerProcess;
	for (int mr = 1; mr < npes; mr++) {
		disps[mr] = disps[mr-1] + counts[mr-1];
		counts[mr] = counts[mr-1];
	}

	counts[npes-1] = (n+1) - (numberPosPerProcess*(npes-1));

	

	std::vector<double> d(n - 1); // this is the d vector used in the Thomas algorithm
	d.reserve(n - 1);
	std::vector<double> cStar(n - 1);// These are the c* and d* vectors used during the forward and backward substitution phases of the algorithm.
	cStar.reserve(n - 1);
	std::vector<double> dStar(n - 1);
	dStar.reserve(n - 1);

	int isLast = 0;
	int isFirst = 0;


	if (myrank == npes-1)
		isLast = 1;

	if (myrank == 0)
		isFirst = 1;

	compTime1 = MPI_Wtime();

	for (int j = 1; j < (0.5/deltaT) + 1; j++) {// the limit for this loop is the number of timesteps needed to get to t = 0.5.

		std::vector<double> x = Tools::createT0Vector(numberPosPerProcess - lastSpaces);
		x.reserve(numberPosPerProcess - lastSpaces);

		double globalX = x[numberPosPerProcess - lastSpaces - 1];
		
		for (int i = 0; i < numberPosPerProcess - lastSpaces - (2*isLast); i++) {

			commTime1 = MPI_Wtime();

			if (i == 0 && myrank != 0 && npes != 1) {
				MPI_Recv(&d[myrank*numberPosPerProcess - 1], 1, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&cStar[myrank*numberPosPerProcess - 1], 1, MPI_DOUBLE, myrank-1, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(&dStar[myrank*numberPosPerProcess - 1], 1, MPI_DOUBLE, myrank-1, 2, MPI_COMM_WORLD, &status);
			}

			commTime2 = MPI_Wtime();

			if (i == 0 && myrank == 0) {// first value of the d vector calculation and simultaneously the forward substitution phase of the algorithm.

				d[myrank*numberPosPerProcess + i] = -aCoef*previousX[myrank*numberPosPerProcess + i] + (1 - C)*previousX[myrank*numberPosPerProcess + i + 1]
					+ -cCoef*previousX[myrank*numberPosPerProcess + i + 2] - aCoef*x[0];
				cStar[myrank*numberPosPerProcess + i] = cCoef / bCoef;
				dStar[myrank*numberPosPerProcess + i] = d[myrank*numberPosPerProcess + i] / bCoef;

			} else if (i == numberPosPerProcess - lastSpaces - (2*isLast) - 1 && myrank == npes-1) {// last value of the d vector calculation and simultaneously the forward substitution phase of the algorithm.

				d[myrank*numberPosPerProcess + i] = -aCoef*previousX[myrank*numberPosPerProcess + i] + (1 - C)*previousX[myrank*numberPosPerProcess + i + 1]
					+ -cCoef*previousX[myrank*numberPosPerProcess + i + 2] - aCoef*x[numberPosPerProcess - lastSpaces - 1];
				cStar[myrank*numberPosPerProcess + i] = 0 / (bCoef - cStar[myrank*numberPosPerProcess + i - 1] * aCoef);
				dStar[myrank*numberPosPerProcess + i] = (d[myrank*numberPosPerProcess + i] - dStar[myrank*numberPosPerProcess + i - 1] * aCoef) / 
					(bCoef - cStar[myrank*numberPosPerProcess + i - 1] * aCoef);

			} else {//intermediate values of the d vector calculation and simultaneously the forward substitution phase of the algorithm.
				
				d[myrank*numberPosPerProcess + i] = -aCoef*previousX[myrank*numberPosPerProcess + i] + (1-C)*previousX[myrank*numberPosPerProcess + i + 1]
					+ -cCoef*previousX[myrank*numberPosPerProcess + i + 2];
				cStar[myrank*numberPosPerProcess + i] = cCoef / (bCoef - cStar[myrank*numberPosPerProcess + i - 1] * aCoef);
				dStar[myrank*numberPosPerProcess + i] = (d[myrank*numberPosPerProcess + i] - dStar[myrank*numberPosPerProcess + i - 1] * aCoef) / 
					(bCoef - cStar[myrank*numberPosPerProcess + i - 1] * aCoef);
				
			}

			commTime3 = MPI_Wtime();

			if (i == numberPosPerProcess - lastSpaces - (2*isLast) - 1 && myrank != npes-1 && npes != 1) {
				MPI_Send(&d[myrank*numberPosPerProcess + i], 1, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD);
				MPI_Send(&cStar[myrank*numberPosPerProcess + i], 1, MPI_DOUBLE, myrank+1, 1, MPI_COMM_WORLD);
				MPI_Send(&dStar[myrank*numberPosPerProcess + i], 1, MPI_DOUBLE, myrank+1, 2, MPI_COMM_WORLD);

			}

			commTime4 = MPI_Wtime();
			
		}

		for (int k = numberPosPerProcess - lastSpaces - (1*isLast) - 1; k > -1 + (isFirst); k--) { // backwards substitution phase to get the final values of the solution for the current timestep.
			
			if (myrank == npes-1 && k == numberPosPerProcess - lastSpaces - (1*isLast) - 1) {

				x[k] = dStar[myrank*numberPosPerProcess + k - 1];
				globalX = x[k];

			} 
			else {

				commTime5 = MPI_Wtime();

				if (k == numberPosPerProcess - lastSpaces - (1*isLast) - 1 && myrank != npes-1 && npes != 1)
				 	MPI_Recv(&globalX, 1, MPI_DOUBLE, myrank+1, 3, MPI_COMM_WORLD, &status);

				commTime6 = MPI_Wtime();

				x[k] = dStar[myrank*numberPosPerProcess + k - 1] - cStar[myrank*numberPosPerProcess + k - 1] * globalX;
				globalX = x[k];

				commTime7 = MPI_Wtime();

				if (k == 0 + (isFirst) && myrank != 0 && npes != 1)
				 	MPI_Send(&globalX, 1, MPI_DOUBLE, myrank-1, 3, MPI_COMM_WORLD);

				commTime8 = MPI_Wtime();
			}
		}

		commTime9 = MPI_Wtime();

		if ((n-1) != 1)
			MPI_Allgatherv(&x[0], numberPosPerProcess - lastSpaces, MPI_DOUBLE, &previousX[0], counts, disps, MPI_DOUBLE, MPI_COMM_WORLD);
		else
			previousX[1] = x[1];

		commTime10 = MPI_Wtime();


		if (myrank == npes-1) {
			if ((j % (int) (0.1/deltaT)) == 0) // We store the values for t = 0.1, 0.2,...,0.5.
				addToAllSolutions(previousX);
		}

	}

	compTime2 = MPI_Wtime();


	totalCommTime += (commTime2 - commTime1) + (commTime4 - commTime3) + (commTime6 - commTime5) +
		 (commTime8 - commTime7) + (commTime10 - commTime9);
	totalCompTime += (compTime2 - compTime1) - totalCommTime;

	
	std::cout << "\n Communication time for p" << myrank << ": " << totalCommTime << "\n";
	std::cout << "\n Computation time for p" << myrank << ": " << totalCompTime << "\n";
}