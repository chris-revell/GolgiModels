/*
 * Functions.cpp
 *
 *  Created on: Apr 12, 2017
 *      Author: jean-patrick
 */

#include "Functions.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <vector>
#include <algorithm>

// Normalize a double vector transforming the exact values in ratio
std::vector<double> normalize(std::vector<double>& vectorToNormalize){
	double normalizationDouble = sumDoublesInVector(vectorToNormalize);
	if (normalizationDouble == 0) {
		normalizationDouble = 1;
	};
	std::vector<double> normalizationVector (vectorToNormalize.size(), 1.0 / normalizationDouble);
	normalizationVector = multiplyDoublesVectors(vectorToNormalize, normalizationVector);
	return normalizationVector;
}

// Print a vector of integers
void printIntsVector(std::vector<int>& vector_to_print){
	for (std::vector<int>::iterator it = vector_to_print.begin() ; it != vector_to_print.end() ; ++it)
		{std::cout << *it << " ";};
	std::cout << std::endl;
}

// Print a vector of doubles
void printDoublesVector(std::vector<double>& vector_to_print){
	for (std::vector<double>::iterator it = vector_to_print.begin() ; it != vector_to_print.end() ; ++it)
		{std::cout << *it << " ";}
	std::cout << std::endl;
}

// Sum integers of a vector and return the sum as an integer
int sumIntegersInVector(std::vector<int>& vectorToSum){
	int sumOfIntegers = 0;
	for(std::vector<int>::iterator it = vectorToSum.begin() ; it != vectorToSum.end(); ++it)
		{sumOfIntegers += *it;};
	return sumOfIntegers;
}

// Sum doubles of a vector and return the sum as an double
double sumDoublesInVector(std::vector<double>& vectorToSum){
	return std::accumulate(vectorToSum.begin(), vectorToSum.end(), 0.0);
}

//Sum two integers vectors member to member
std::vector<int> sumIntVectors(std::vector<int> firstIntVector, std::vector<int>& secondIntVector){
	std::transform(firstIntVector.begin(), firstIntVector.end(), // range of first vector
			secondIntVector.begin(), 					 // beginning of second vector
			firstIntVector.begin(), std::plus<int>());    // where we store and what operation
	return firstIntVector;
}

//Multiply two double vectors member to member
std::vector<double> multiplyDoublesVectors(std::vector<double> firstDoublesVector, std::vector<double>& secondDoublesVector){
	std::transform(firstDoublesVector.begin(), firstDoublesVector.end(),
			secondDoublesVector.begin(),
			firstDoublesVector.begin(), std::multiplies<double>());
	return firstDoublesVector;
}

// Use a random number to pick an index in a vector
int chooseIndex(double& inputDouble, std::vector<double> doubleVector){
	doubleVector = normalize(doubleVector);
	size_t index = 0;
	while (index < doubleVector.size() && inputDouble >= doubleVector[index]) {
		inputDouble -= doubleVector[index];
		index++;};
	// We prevent some side effects
	if (index >= doubleVector.size()) {
		index = doubleVector.size() - 1;
		while (doubleVector[index] == 0) {
			index--;
		};
	};
	if (index < 0) {
		std::cout << "ERROR! cannot find a correct index" << std::endl;
	};
	return index;
}

void chooseIndexMatrix(double &inputDouble, double &inputDouble2, std::vector<std::vector<double> >& doubleMatrix, std::vector<double>& vectorTool, double &doubleTool, int &firstIndex, int &secondIndex){
	vectorTool.assign(doubleMatrix.size(), 0);
	for (size_t i = 0; i < doubleMatrix.size(); ++i) {
		vectorTool[i] = sumDoublesInVector(doubleMatrix[i]);
	}
	secondIndex = chooseIndex(inputDouble, vectorTool);
	vectorTool.assign(doubleMatrix[secondIndex].begin(), doubleMatrix[secondIndex].end());
	// We take the second element:
	vectorTool = doubleMatrix[secondIndex];
	firstIndex = chooseIndex(inputDouble2, vectorTool);
}

// Calculate the homology between two normalized identity vectors
double homology(std::vector<double> firstVector, std::vector<double> secondVector){
	if (firstVector.size() != secondVector.size()) {
		std::cout << "ERROR! I can't calculate the homology" << '\n';
	};
	firstVector = normalize(firstVector);
	secondVector = normalize(secondVector);
	firstVector = multiplyDoublesVectors(firstVector, secondVector);
	return sumDoublesInVector(firstVector);
}

// create a folder using the current time and give the path back
std::string createFolderWithDate(std::string string){
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[100];
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer,80,"%Y%m%d-%I.%M.%S",timeinfo);
	std::string str(buffer);
	struct stat st = {0};
	if (stat(string.c_str(), &st) == -1) {
		mkdir(string.c_str(), 0777);
	};
	std::string path = string + str + "/";
	if (stat(path.c_str(), &st) == -1) {
		mkdir(path.c_str(), 0777);};
	return path;
}
