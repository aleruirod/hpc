#ifndef CRANKNICHOLSON_H
#define CRANKNICHOLSON_H

#include "Solution.h"
#include <cmath>
#include <mpi.h>
#include <iostream>


class CrankNicholsonMethod : public Solution {

public:

	// CONSTRUCTORS
	/**
	* Creates an empty CrankNicholsonMethod object
	* @see CrankNicholsonMethod(Vector<Vector<double>> sols);
	*/

	CrankNicholsonMethod();

	/**
	* Creates an CrankNicholsonMethod object from a vector of double vectors.
	* @see CrankNicholsonMethod();
	* @param sols Vector of double vectors that the new Solution will use.
	*/

	CrankNicholsonMethod(std::vector<std::vector<double>> sols /**< std::vector<std::vector<double>>. Vector of double vectors that the new Solution will use. */);

	//COMPUTATION
	/**
	* Computes and stores the values for the Solution using this method.
	*/

	void compute();

protected:
	double C = (deltaT*D) / pow(deltaX, 2); //For the Laasonen method, the C coefficient is constant throughout the whole problem.
	double aCoef = -C/2; //!< the value for the a coefficient is defined when studying the formation of the tridiagonal matrix.
	double bCoef = (1 + C); //!< the value for the b coefficient is defined when studying the formation of the tridiagonal matrix.
	double cCoef = -C/2; //!< the value for the c coefficient is defined when studying the formation of the tridiagonal matrix.
};
#endif
