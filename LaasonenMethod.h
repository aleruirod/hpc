#ifndef LAASONEN_H
#define LAASONEN_H

#include "Solution.h"
#include "Output.h"
#include <cmath>
#include <iostream>
#include <mpi.h>


class LaasonenMethod : public Solution {

public:

	// CONSTRUCTORS
	/**
	* Creates an empty LaasonenMethod object
	* @see LaasonenMethod(Vector<Vector<double>> sols);
	*/

	LaasonenMethod();

	/**
	* Creates an LaasonenMethod object from a vector of double vectors.
	* @see LaasonenMethod();
	* @param sols Vector of double vectors that the new Solution will use.
	*/

	LaasonenMethod(std::vector<std::vector<double>> sols /**< std::vector<std::vector<double>>. Vector of double vectors that the new Solution will use. */);

	//COMPUTATION
	/**
	* Computes and stores the values for the Solution using this method.
	*/

	void compute();

protected:
	double C = (deltaT*D) / pow(deltaX, 2); //For the Laasonen method, the C coefficient is constant throughout the whole problem.
	double aCoef = -C; //!< the value for the a coefficient is defined when studying the formation of the tridiagonal matrix.
	double bCoef = (2 * C + 1); //!< the value for the b coefficient is defined when studying the formation of the tridiagonal matrix.
	double cCoef = -C; //!< the value for the c coefficient is defined when studying the formation of the tridiagonal matrix.
};
#endif