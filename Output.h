#ifndef OUTPUT_H
#define OUTPUT_H

#include <vector>
#include <string>
#include <mpi.h>
#include <iostream>
#include "Solution.h"//we will print Solution objects.


/**
* The Output class is utilised to export the content
* \n of the solutions as files or simply show the 
* \n results on the console.
*
*/

class Output {
	public:

		/**
		* Prints on console all the solutions stored in a Solution object.
		* @param sols This should be the allSolutions attribute from the Solution object we want to print.
		*/

		static void printSolution(std::vector<std::vector<double>> sols /**< std::vector<std::vector<double>>. Vector of double vectors that we want to print. */);

		/**
		* Exports the Solution in several .dat files.
		* @param sol this is the Solution we want to export the results from.
		* @param tSize this is the time difference between the different timesteps our solution has.
		* @param name this is the base name for the method used to get the solution we want to export.
		*/

		static void exportSolution(Solution sol /**Solution. this is the Solution we want to export the results from. */, 
			double tSize /**double. this is the time difference between the different timesteps our solution has. */,
			std::string name /**std::string. this is the base name for the method used to get the solution we want to export. */);

		/**
		* Prints a vector's contents on console.
		* @param v this is the vector we want to print.
		*/

		static void printVector(std::vector<double> v /**std::vector<double>. Vector we want to print on screen. */);

};
#endif
