#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <cmath>
#include <iostream>
#include "Solution.h"

/**
* The Tools class contains functions to support the computation
* \n and processing of results.
*/

class Tools {

	public:

		/**
		* Creates a vector following the boundary conditions of the problem.
		* @param n Length of the vector.
		*/

		static std::vector<double> createT0Vector(int n /**int. Length of the vector. */);

		/**
		* Calculates the two norm for the difference between a certain Solution and the Analytical Solution.
		* @param analytSol Analytical Solution.
		* @param sol Solution object we want to compare.
		*/

		static std::vector<double> twoNorm(Solution analytSol /**Solution. Analytical Solution */,
			Solution sol /**Solution. Solution object we want to compare with the Analytical Solution. */);

};
#endif
