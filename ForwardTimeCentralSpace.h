#ifndef FTCS_H
#define FTCS_H

#include "Solution.h"
#include "Output.h"
#include <cmath>
#include <mpi.h>
#include <iostream>



/**
* ForwardTimeCentralSpace
*/

class ForwardTimeCentralSpace : public Solution {

	public:

		// CONSTRUCTORS
		/**
		* Creates an empty ForwardTimeCentralSpace object
		* @see ForwardTimeCentralSpace(Vector<Vector<double>> sols);
		*/

		ForwardTimeCentralSpace();

		/**
		* Creates an ForwardTimeCentralSpace object from a vector of double vectors.
		* @see ForwardTimeCentralSpace();
		* @param sols Vector of double vectors that the new Solution will use.
		*/

		ForwardTimeCentralSpace(std::vector<std::vector<double>> sols /**< std::vector<std::vector<double>>. Vector of double vectors that the new Solution will use. */);

		// COMPUTATION
		/**
		* Computes and stores the values for the Solution using this method.
		*/

		void compute();
};
#endif