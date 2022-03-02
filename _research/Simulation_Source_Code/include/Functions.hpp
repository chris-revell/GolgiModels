/*
 * Functions.hpp
 *
 *  Created on: Apr 12, 2017
 *      Author: jean-patrick
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_

#include <iostream>
#include <vector>
#include <numeric>

// Normalize a double vector transforming the exact values in ratio
std::vector<double> normalize(std::vector<double>& vectorToNormalize);

// Sum integers of a vector and return the sum as an integer
int sumIntegersInVector(std::vector<int>& vectorToSum);

// Sum doubles of a vector and return the sum as an double
double sumDoublesInVector(std::vector<double>& vectorToSum);

//Sum two integers vectors member to member
std::vector<int> sumIntVectors(std::vector<int> firstIntVector,
		std::vector<int>& secondIntVector);

//Multiply two double vectors member to member
std::vector<double> multiplyDoublesVectors(std::vector<double> firstDoublesVector,
		std::vector<double>& secondDoublesVecor);

// Print a vector of integers
void printIntsVector(std::vector<int>& vector_to_print);

// Print a vector of doubles
void printDoublesVector(std::vector<double>& vector_to_print);

// Use a random number to pick an index in a vector
int chooseIndex(double& inputDouble, std::vector<double> doubleVector);
void chooseIndexMatrix(double &inputDouble, double &inputDouble2, std::vector<std::vector<double> >& doubleMatrix, std::vector<double>& vectorTool, double &doubleTool, int &firstIndex, int &secondIndex);

// Calculate the homology between two normalized identity vectors
double homology(std::vector<double> firstVector, std::vector<double> secondVector);

// create a folder using the current time and give the path back
std::string createFolderWithDate(std::string string);


#endif /* FUNCTIONS_HPP_ */
