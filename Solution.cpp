#include "Solution.h"
#include "Tools.h"

// CONSTRUCTORS
/*
* Default constructor (empty solution)
*/

Solution::Solution() {}

/*
* Constructor from a vector of double vectors.
*/

Solution::Solution(std::vector<std::vector<double>> sols) {
	Solution::allSolutions = sols;
}

// GETTERS
/**
* Returns all solutions stored in a Solution object
*/

std::vector<std::vector<double>> Solution::getAllSolutions() {
	return Solution::allSolutions;
}

/**
* Returns the value of deltaX for the current solution.
*/

double Solution::getDeltaX() {
	return deltaX;
}

/**
* Returns the value of deltaT for the current solution.
*/

double Solution::getDeltaT() {
	return deltaT;
}

//SETTERS
/**
* Changes the value of deltaX for the current solution.
*/

void Solution::setDeltaX(double x) {
	deltaX = x;
}

/**
* Changes the value of deltaT for the current solution.
*/

void Solution::setDeltaT(double t) {
	deltaT = t;
}

//STORING
/**
* adds a new solution to allSolutions.
*/

void Solution::addToAllSolutions(std::vector<double> v) {
	Solution::allSolutions.resize(allSolPos + 1);// we change the size of the vector of vectors so a new solution can fit in it.
	Solution::allSolutions[allSolPos] = v;// we add the next solution we want to store.
	Solution::allSolPos++;// we increase the counter that tracks the size of the vector of vectors.
}