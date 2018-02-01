#include <iostream>
#include <fstream> 

#include "Output.h"

/**
* Print all solutions from a vector of double vectors.
*/

void Output::printSolution(std::vector<std::vector<double>> sols) {

	for (size_t i = 0; i < sols.size(); i++)
	{
		for (size_t j = 0; j < sols[i].size(); j++)
		{
			std::cout << sols[i][j] << " "; // for every value in a timestep, we introduce a blank space as separation.
		}
		std::cout << "\n"; // we introduce a line break between each different solution vector.
	}
}

/**
* Exports all results stored in a Solution object.
*/

void Output::exportSolution(Solution sol, double tSize, std::string name) {

	int npes;

	MPI_Comm_size(MPI_COMM_WORLD, &npes);

	for (size_t i = 0; i < sol.getAllSolutions().size(); i++) {

		std::string aux = name;
		aux += std::to_string(tSize*i)+"DeltaX_"+std::to_string(sol.getDeltaX())+"numProcesses_"+std::to_string(npes)+".dat";

		std::ofstream outfile(aux);

		for (size_t j = 0; j < sol.getAllSolutions()[0].size(); j++)
			outfile << sol.getDeltaX()*j << " " << sol.getAllSolutions()[i][j] << std::endl;
			
		outfile.close();

	}

}

/**
* Prints a concrete vector on screen.
*/

void Output::printVector(std::vector<double> v) {

	for (size_t i = 0; i < v.size(); i++)
		std::cout << v[i] << "\n";

}