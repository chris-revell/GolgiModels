/*
 * TimeFlow.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: jean-patrick
 */

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "TimeFlow.hpp"
#include "OrganelleNetwork.hpp"
#include "Functions.hpp"

namespace timing {

TimeFlow::TimeFlow(std::string parentDirectory,
		double influxRateSimu,
		double buddingRateSimu,
		std::vector<double> maturationRatesSimu,
		std::vector<std::vector<double> > exitIdentitiesSimu,
		std::vector<std::vector<double> > membraneCargoesAffinitySimu,
		std::vector<std::vector<double> > buddingCargoesAffinitySimu,
		int doWeFollowExitSimu,
		std::vector<int> cargoesCompositionSimu
		):organelleNetwork(influxRateSimu,
				buddingRateSimu,
				maturationRatesSimu,
				exitIdentitiesSimu,
				membraneCargoesAffinitySimu,
				buddingCargoesAffinitySimu,
				doWeFollowExitSimu){

	simulationPath = createFolderWithDate(parentDirectory);
	influxRate = influxRateSimu;
	buddingRate = buddingRateSimu;
	maturationRates = maturationRatesSimu;
	exitIdentities = exitIdentitiesSimu;
	membraneCargoesAffinity = membraneCargoesAffinitySimu;
	buddingCargoesAffinity = buddingCargoesAffinitySimu;
	doWeFollowExit = doWeFollowExitSimu;
	doWeFollowCargoes = 0;
	doWeInjectCargoes = 0;
	cargoesComposition = cargoesCompositionSimu;
	accumulatedSize = 0;

	// RANDOM NUMBER GENERATOR:
	seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	mt_rand = std::mt19937(seed);
	disint = std::uniform_real_distribution<double>(0, 1);

	// We save the parameters file:
	std::ofstream outputFile (simulationPath + "parameters.txt");
	std::stringstream seedString;
	seedString << std::setprecision(1000);
	seedString << seed;
	outputFile << "seed = " << seedString.str() << std::endl;
	outputFile << "budding rate = " << buddingRate << std::endl;
	outputFile << "exit identities = |";
	if (exitIdentities.size() != 0) {
		for (size_t i = 0; i < exitIdentities.size() - 1; i++) {
			outputFile << "|";
			for (size_t j = 0; j < exitIdentities[i].size() - 1; j++) {
				outputFile << exitIdentities[i][j] << ", ";
			};
			outputFile << exitIdentities[i].back() << "|,";
		};
		outputFile << "|";
		for (size_t j = 0; j < exitIdentities.back().size() - 1; j++) {
			outputFile << exitIdentities.back()[j] << ", ";
		};
		outputFile << exitIdentities.back().back() << "||" << std::endl;
	} else {
		outputFile << " |" << std::endl;
	}
	outputFile << "maturation rates = |";
	if (maturationRates.size() != 0) {
		for (size_t i = 0; i < maturationRates.size() - 1; i++) {
			outputFile << maturationRates[i] << ", ";
		}
		outputFile << maturationRates.back() << "|" << std::endl;
	} else {
		outputFile << " |" << std::endl;
	}
	outputFile << "influx rate = " << influxRate << std::endl;
	outputFile << "membrane cargoes affinity = |";
	if (membraneCargoesAffinity.size() != 0) {
		for (size_t i = 0; i < membraneCargoesAffinity.size() - 1; i++) {
			outputFile << "|";
			for (size_t j = 0; j < membraneCargoesAffinity[i].size() - 1; j++) {
				outputFile << membraneCargoesAffinity[i][j] << ", ";
			};
			outputFile << membraneCargoesAffinity[i].back() << "|,";
		};
		outputFile << "|";
		for (size_t j = 0; j < membraneCargoesAffinity.back().size() - 1; j++) {
			outputFile << membraneCargoesAffinity.back()[j] << ", ";
		};
		outputFile << membraneCargoesAffinity.back().back() << "||" << std::endl;
	} else {
		outputFile << " |" << std::endl;
	}
	outputFile << "budding cargoes affinity = |";
	if (buddingCargoesAffinity.size() != 0) {
		for (size_t i = 0; i < buddingCargoesAffinity.size() - 1; i++) {
			outputFile << "|";
			for (size_t j = 0; j < buddingCargoesAffinity[i].size() - 1; j++) {
				outputFile << buddingCargoesAffinity[i][j] << ", ";
			};
			outputFile << buddingCargoesAffinity[i].back() << "|,";
		};
		outputFile << "|";
		for (size_t j = 0; j < buddingCargoesAffinity.back().size() - 1; j++) {
			outputFile << buddingCargoesAffinity.back()[j] << ", ";
		};
		outputFile << buddingCargoesAffinity.back().back() << "||" << std::endl;
	} else {
		outputFile << " |" << std::endl;
	}


	nextEvent = 0;
	currentTime = 0;
	numberOfTimePoints = 0;
	eventsCounting = {1, 0, 0, 0, 0};
	isThereACargo.assign(membraneCargoesAffinity.size(), 0);
	organelleNetwork.getGlobalRates(globalNetworkRates);

	intToolOne = 0;
	intVectorToolOne = {};
	intVectorToolTwo = {};
	doubleVectorToolOne = {};
	doubleVectorToolTwo = {};
	intMatrixToolOne = {{}};
}

TimeFlow::~TimeFlow() {}

// METHODS:

			/********************************/
			/*		SETTERS AND GETTERS		*/
			/********************************/

// TIME following:
void TimeFlow::getNumberOfTimePoints(int &outputInt) const{
	outputInt = numberOfTimePoints;
}
void TimeFlow::getCurrentTime(double &outputDouble) const{
	outputDouble = currentTime;
}

// EVENTS counting:
void TimeFlow::getCountingOfEvent(std::vector<int>& outputIntVector) const{
	outputIntVector = eventsCounting;
}

// CARGOES following:
void TimeFlow::setDoWeFollowCargoes(int inputInt){
	doWeFollowCargoes = inputInt;
	std::ofstream outputFile;
	std::string temporaryPath = simulationPath + "cargoLife.txt";
	outputFile.open(temporaryPath.c_str(), std::ios::out | std::ios::app);
	if (inputInt == 1) {
		outputFile << currentTime << " ; ON \n";
	} else {
		outputFile << currentTime << " ; OFF \n";
	}
	outputFile.close();
}
void TimeFlow::setDoWeInjectCargoes(int inputInt){
	doWeInjectCargoes = inputInt;
	std::ofstream outputFile;
	std::string temporaryPath = simulationPath + "snapshotCargoes.txt";
	outputFile.open(temporaryPath.c_str(), std::ios::out | std::ios::app);
	if (inputInt == 1) {
		outputFile << currentTime << " ; ON\n";
	} else {
		outputFile << currentTime << " ; OFF\n";
	}
	outputFile.close();
}

			/****************************/
			/*			EVENTS			*/
			/****************************/

// Time updating using a Gillespie algorithm
int TimeFlow::updateTime(){
	if (sumDoublesInVector(globalNetworkRates) == 0) {
		return 1;
	} else {
		currentTime += -log(disint(mt_rand))/sumDoublesInVector(globalNetworkRates);
		numberOfTimePoints++;
		return 0;
	}
}

// Choose what event occures at this time point
void TimeFlow::chooseTypeOfEvent(){
	auto a = disint(mt_rand);
	organelleNetwork.getGlobalRates(globalNetworkRates);
	nextEvent = chooseIndex(a, globalNetworkRates);
	eventsCounting[nextEvent]++;
}
void TimeFlow::setTypeOfEvent(int inputInt){
	nextEvent = inputInt;
}
void TimeFlow::getTypeOfEvent(int &outputInt) const{
	outputInt = nextEvent;
}
void TimeFlow::applyModification(){
	if (nextEvent == 0) {
		intToolOne = isThereACargo.size();
		if (doWeFollowCargoes == 1 && sumIntegersInVector(isThereACargo) != intToolOne) {
			intVectorToolOne.assign(isThereACargo.size(), 0);
			intVectorToolTwo.assign(isThereACargo.size(), 0);
			for (size_t i = 0; i < isThereACargo.size(); ++i) {
				if (isThereACargo[i] == 0) {
					intVectorToolTwo[i] = 1;
				}
			}
			if (doWeInjectCargoes == 1) {
				intToolOne = organelleNetwork.injectCargo(cargoesComposition, intVectorToolTwo);
			} else {
				intToolOne = organelleNetwork.injectCargo(intVectorToolOne, intVectorToolTwo);
			}
			isThereACargo.assign(isThereACargo.size(), 1);
			organelleNetwork.followCargoLife(nextEvent, intToolOne, simulationPath + "cargoLife.txt", currentTime);
		} else {
			if (doWeInjectCargoes == 1) {
				intVectorToolOne.assign(cargoesComposition.size(), 0);
				intToolOne = organelleNetwork.injectCargo(cargoesComposition, intVectorToolOne);
			} else {
				organelleNetwork.vesicleInjection();
			}
		}
	}
	else if (nextEvent == 1) {
		intToolOne = organelleNetwork.elementsFusion(disint(mt_rand), disint(mt_rand));
		if (intToolOne != -1 && doWeFollowCargoes == 1) {
			organelleNetwork.followCargoLife(nextEvent, intToolOne, simulationPath + "cargoLife.txt", currentTime);
		}
	}
	else if (nextEvent == 2) {
		intMatrixToolOne = organelleNetwork.elementExit(disint(mt_rand), disint(mt_rand), simulationPath + "snapshotExit.txt", currentTime);
		if (doWeFollowCargoes == 1 && intMatrixToolOne[0][0] != -1) {
			for (size_t i = 0; i < isThereACargo.size(); ++i) {
				isThereACargo[i] -= intMatrixToolOne[1][i];
			}
			std::ofstream outputFile;
			outputFile.open((simulationPath + "cargoLife.txt").c_str(), std::ios::out | std::ios::app);
			outputFile << currentTime << " ; " << nextEvent << " ; ";
			for (size_t i = 0; i < maturationRates.size(); ++i) {
				outputFile << exitIdentities[intMatrixToolOne[0][0]][i] << " | ";

			}
			outputFile << exitIdentities[intMatrixToolOne[0][0]].back() << " ; ";
			for (size_t i = 0; i < intMatrixToolOne[1].size() - 1; ++i) {
				outputFile << intMatrixToolOne[1][i] << " | ";
			}
			outputFile <<intMatrixToolOne[1].back();
			outputFile << "\n";
			outputFile.close();
		}
	}
	else if (nextEvent == 3) {
		doubleVectorToolOne.assign(membraneCargoesAffinity.size(), 0);
		doubleVectorToolTwo.assign(membraneCargoesAffinity.size(), 0);
		for (size_t i = 0; i < doubleVectorToolOne.size(); ++i) {
			doubleVectorToolOne[i] = disint(mt_rand);
			doubleVectorToolTwo[i] = disint(mt_rand);
		}
		intToolOne = organelleNetwork.vesicleBudding(disint(mt_rand), disint(mt_rand), doubleVectorToolOne, doubleVectorToolTwo);
		if (intToolOne != -1 && doWeFollowCargoes == 1) {
			organelleNetwork.followCargoLife(nextEvent, intToolOne, simulationPath + "cargoLife.txt", currentTime);
		}
	}
	else if (nextEvent == 4){
		intToolOne = organelleNetwork.specyMaturation(disint(mt_rand), disint(mt_rand));
		if (intToolOne != -1 && doWeFollowCargoes == 1) {
			organelleNetwork.followCargoLife(nextEvent, intToolOne, simulationPath + "cargoLife.txt", currentTime);
		}
	}
	else {
		std::cout << "ERROR! cannot choose event" << '\n';
	};
}

// Accumulated size promoter getter and reseter:
void TimeFlow::accumulateSize(){
	organelleNetwork.getElementsSize(intVectorToolOne);
	accumulatedSize += sumIntegersInVector(intVectorToolOne);
}
void TimeFlow::getAccumulatedSize(int &outputInt) const{
	outputInt = accumulatedSize;
}
void TimeFlow::resetAccumulatedSize(){
	accumulatedSize = 0;
}

// Follow parameters during the simulation:
// snapshot the total organization:
void TimeFlow::initializeSnapshotOrganisation(){
	std::ofstream outputFile (simulationPath + "snapshotOrganization.txt");
	outputFile << "Time ; Organization" << std::endl;
	outputFile.close();
}
void TimeFlow::snapshotOrganisation(){
	organelleNetwork.snapshotOrganisation(simulationPath + "snapshotOrganization.txt", currentTime);
}
// follow the size distribution:
void TimeFlow::initializeFollowSizes(){
	std::ofstream outputFile (simulationPath + "snapshotSizeDistribution.txt");
	outputFile << "Time ; Sizes" << std::endl;
	outputFile.close();
}
void TimeFlow::followSizes(){
	organelleNetwork.followSizes(simulationPath + "snapshotSizeDistribution.txt", currentTime);
}
// follow the network rates:
void TimeFlow::initializeFollowRates(){
	std::ofstream outputFile (simulationPath + "snapshotRates.txt");
	outputFile << "Time ; Injection ; Fusion ; Exit ; Budding ; Maturation ; Purity" << std::endl;
	outputFile.close();
}
void TimeFlow::followRates(){
	organelleNetwork.followRates(simulationPath + "snapshotRates.txt", currentTime);
}
// follow the global network identities:
void TimeFlow::initializeFollowIdentity(){
	std::ofstream outputFile (simulationPath + "snapshotIdentities.txt");
	std::string ABC = "ABCDEFGHIJKLM";
	outputFile << "Time ; ";
	for (size_t i = 0; i < maturationRates.size() + 1; i++) {
		outputFile << ABC[i] << " ; ";
	}
	outputFile << "Purity" << std::endl;
	outputFile.close();
}
void TimeFlow::followIdentity(){
	organelleNetwork.followIdentity(simulationPath + "snapshotIdentities.txt", currentTime);
}
// follow the exit composition:
void TimeFlow::initializeFollowExit(){
	std::ofstream outputFile (simulationPath + "snapshotExit.txt");
	outputFile << "Time ; Exit ; Composition" << std::endl;
	outputFile.close();
}
// follow the cargoes:
void TimeFlow::initializeFollowCargoLife(){
	std::ofstream outputFile (simulationPath + "cargoLife.txt");
	outputFile << "Time ; Event ; Composition" << std::endl;
	outputFile.close();
}
void TimeFlow::initializeFollowCargoes(){
	std::ofstream outputFile (simulationPath + "snapshotCargoes.txt");
	outputFile << "Time ; Cargoes" << std::endl;
	outputFile.close();
}
void TimeFlow::followCargoes(){
	organelleNetwork.followCargoes(simulationPath + "snapshotCargoes.txt", currentTime);
}
} /* namespace timing */
