#include "Tools.h"

/**
* Creates a vector at t = 0.
*/

std::vector<double> Tools::createT0Vector(int n) {

	std::vector<double> t0Vec(n, 100.0);

	t0Vec[0] = 300.0;
	t0Vec[n - 1] = 300.0;

	return t0Vec;

}

/**
* Applies the two norm on two solutions to see how they differ. One of them usually being the Analytical Solution.
*/

std::vector<double> Tools::twoNorm(Solution analytSol, Solution sol) {

	std::vector<std::vector<double>> allAnalyticalSols = analytSol.getAllSolutions();
	std::vector<std::vector<double>> allSolutionsToCompare = sol.getAllSolutions();

	std::vector<double> twoNormResult(allAnalyticalSols.size());

	for (std::size_t i = 0; i < allAnalyticalSols.size(); i++) {

		std::vector<double> analyticalSoln = allAnalyticalSols[i];
		std::vector<double> solutionToCompare = allSolutionsToCompare[i];

		double sum = 0.0;

		for (std::size_t j = 0; j < analyticalSoln.size(); j++) {

			//we add the square of the differences between each value of the analytical solution and the one we are comparing.
			sum += (solutionToCompare[j] - analyticalSoln[j]) * (solutionToCompare[j] - analyticalSoln[j]);
		}

		twoNormResult[i] = sqrt(sum);

	}
	
	return twoNormResult;
}