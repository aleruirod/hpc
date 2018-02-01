#ifndef ANALYTICAL_H
#define ANALYTICAL_H

#include "Solution.h"
#include <iostream>


/**
* the Analytical Solution of the problem gives us the exact right
* \n  values for our case. This solution will be used to compare
* \n the values we get from the other methods used.
*/

class AnalyticalSolution : public Solution {
	public:

		// CONSTRUCTORS
		/**
		* Creates an empty AnalyticalSolution object
		* @see AnalyticalSolution(Vector<Vector<double>> sols);
		*/

		AnalyticalSolution();

		/**
		* Creates an AnalyticalSolution object from a vector of double vectors.
		* @see AnalyticalSolution();
		* @param sols Vector of double vectors that the new Solution will use.
		*/

		AnalyticalSolution(std::vector<std::vector<double>> sols /**< std::vector<std::vector<double>>. Vector of double vectors that the new Solution will use. */);

		//COMPUTATION
		/**
		* Computes and stores the values for the Solution using this method.
		*/

		void compute();

	protected:
		const double PI = 3.141592; //!< Pi is needed in order to calculate the analytical solution of the problem.
};

#endif
