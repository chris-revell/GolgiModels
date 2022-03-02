/*
 * OrganelleNetwork.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: jean-patrick
 */

#include <fstream>

#include "OrganelleNetwork.hpp"
#include "OrganelleElement.hpp"
#include "Functions.hpp"



namespace network {

OrganelleNetwork::OrganelleNetwork(double influxRateSimu,
		double buddingRateSimu,
		std::vector<double> maturationRatesSimu,
		std::vector<std::vector<double> > exitIdentitiesSimu,
		std::vector<std::vector<double> > membraneCargoesAffinitySimu,
		std::vector<std::vector<double> > buddingCargoesAffinitySimu,
		int doWeFollowExitSimu
		){

	influxRate = influxRateSimu;
	buddingRate = buddingRateSimu;
	maturationRates = maturationRatesSimu;
	exitIdentities = exitIdentitiesSimu;
	membraneCargoesAffinity = membraneCargoesAffinitySimu;
	buddingCargoesAffinity = buddingCargoesAffinitySimu;
	doWeFollowExit = doWeFollowExitSimu;

	networkOrganisation = {};
	elementsSizes.assign(0, 0);
	elementsIdentities.assign(0, {});
	cargoesOrganization.assign(0, {});
	followedCargoOrganization.assign(0, {});
	elementsPurities.assign(0, 0);
	totalNetworkPurity = 1;

	elementsFusionRates.assign(0, {});
	totalNetworkFusionRate = 0;
	elementsExitRates.assign(0, 0);
	totalNetworkExitRate = 0;
	elementsBuddingRates.assign(0, 0);
	totalNetworkBuddingRate = 0;
	elementsMaturationRates.assign(0, 0);
	totalNetworkMaturationRate = 0;
	globalNetworkRates = {influxRateSimu, 0, 0, 0, 0};

	intVectorToolOne = {};
	intVectorToolTwo = {};
	doubleVectorToolOne = {};
	doubleVectorToolTwo = {};
	doubleMatrixTool = {{}};
	firstIndex = 0;
	secondIndex = 0;
	intVectorToolOneFunc = {};
	intVectorToolTwoFunc = {};
	doubleVectorToolOneFunc = {};
	doubleVectorToolTwoFunc = {};
	doubleToolOneFunc = 0;
	doubleToolTwoFunc = 0;
	intMatrixToolFunc = {{}};
	doubleMatrixToolFunc= {{}};
	firstIndexFunc = 0;
	secondIndexFunc = 0;
	intToolOneFunc = 0;
	intToolTwoFunc = 0;

	// We calculate Na, Nb and Nc:
	doubleToolOne = influxRate / (maturationRates[0] + exitIdentities[0][0]);
	doubleToolTwo = doubleToolOne * maturationRates[0] / exitIdentities.back().back();
	intToolOne = doubleToolOne;
	intToolTwo = doubleToolTwo;

	// We generate a vesicular Golgi with the good identities:
	for (int i = 0; i < intToolTwo; ++i) {
		vesicleInjection();
		specyMaturation(1, 1);
		specyMaturation(1, 1);
	}
	for (int i = 0; i < intToolOne; ++i) {
		vesicleInjection();
		specyMaturation(1, 1);
	}
	for (int i = 0; i < intToolOne; ++i) {
		vesicleInjection();
	}
}

OrganelleNetwork::~OrganelleNetwork() {}

// METHODS:
				/********************************/
				/*		SETTERS AND GETTERS		*/
				/********************************/


// SIZE updater (if needed) and getter.
void OrganelleNetwork::updateElementsSize(){
	elementsSizes.assign(networkOrganisation.size(), 0);
	for (size_t i = 0; i < networkOrganisation.size(); i++) {
		networkOrganisation[i].getElementSize(elementsSizes[i]);
	}
}
void OrganelleNetwork::getElementsSize(std::vector<int>& outputIntVector) const{
	outputIntVector = elementsSizes;
}

// IDENTITIES getter and updater:
void OrganelleNetwork::updateElementsIdentities(){
	elementsIdentities.assign(networkOrganisation.size(), {});
	for (size_t i = 0; i < networkOrganisation.size(); i++) {
		networkOrganisation[i].getElementIdentity(elementsIdentities[i]);
	}
}
void OrganelleNetwork::getElementsIdentities(std::vector<std::vector<int> >& outputIntMatrix) const{
	outputIntMatrix = elementsIdentities;
}

// CARGOES getter and setter:
void OrganelleNetwork::updateCargoesOrganization(){
	cargoesOrganization.assign(networkOrganisation.size(), {});
	for (size_t i = 0; i < networkOrganisation.size(); ++i) {
		networkOrganisation[i].getElementCargoes(cargoesOrganization[i]);
	}
}
void OrganelleNetwork::getCargoesOrganization(std::vector<std::vector<int> >& outputIntMatrix) const{
	outputIntMatrix = cargoesOrganization;
}
void OrganelleNetwork::updateFollowedCargoesOrganization(){
	followedCargoOrganization.assign(networkOrganisation.size(), {});
	for (size_t i = 0; i < networkOrganisation.size(); ++i) {
		networkOrganisation[i].getFollowedCargoes(followedCargoOrganization[i]);
	}
}
void OrganelleNetwork::getFollowedCargoesOrganization(std::vector<std::vector<int> >& outputIntMatrix) const{
	outputIntMatrix = followedCargoOrganization;
}

// PURITY of every element:
void OrganelleNetwork::updateElementsPurities(){
	elementsPurities.assign(networkOrganisation.size(), 0);
	for (size_t i = 0; i < networkOrganisation.size(); ++i) {
		networkOrganisation[i].updateElementRates();
		networkOrganisation[i].getElementPurity(elementsPurities[i]);
	}
}
void OrganelleNetwork::getElementsPurities(std::vector<double>& outputDoubleVector) const{
	outputDoubleVector = elementsPurities;
}
void OrganelleNetwork::updateTotalNetworkPurity(){
	totalNetworkPurity = 0;
	doubleToolTwoFunc = 0;
	for (size_t i = 0; i < elementsPurities.size(); ++i) {
		doubleToolOneFunc = elementsSizes[i];
		if (doubleToolOneFunc >= 2) {
			totalNetworkPurity += elementsPurities[i] * doubleToolOneFunc;
			doubleToolTwoFunc += doubleToolOneFunc;
		}
	}
	if (doubleToolTwoFunc != 0) {
		totalNetworkPurity /= doubleToolTwoFunc;
	} else{
		totalNetworkPurity = 1;
	}
}
void OrganelleNetwork::getTotalNetworPurity(double &outputDouble) const{
	outputDouble = totalNetworkPurity;
}

// FUSION rates between two elements.
void OrganelleNetwork::updateElementsFusionRates(){
	elementsFusionRates.assign(networkOrganisation.size(), {});
	for (size_t i = 0; i < networkOrganisation.size(); i++) {
		elementsFusionRates[i].assign(i, {});
	}
	for (size_t i = 0; i < networkOrganisation.size(); i++) {
		for (size_t j = 0; j < networkOrganisation.size(); j++) {
			if (i > j) {
				networkOrganisation[i].getElementIdentity(intVectorToolOneFunc);
				networkOrganisation[j].getElementIdentity(intVectorToolTwoFunc);
				doubleVectorToolOneFunc.assign(intVectorToolOneFunc.begin(), intVectorToolOneFunc.end());
				doubleVectorToolTwoFunc.assign(intVectorToolTwoFunc.begin(), intVectorToolTwoFunc.end());
				elementsFusionRates[i][j] = homology(doubleVectorToolOneFunc, doubleVectorToolTwoFunc);
			}
		}
	}
}
void OrganelleNetwork::getElementsFusionRates(std::vector<std::vector<double> >& outputMatrix) const{
	outputMatrix = elementsFusionRates;
}
void OrganelleNetwork::updateTotalNetworkFusionRate(){
	totalNetworkFusionRate = 0;
	for (size_t i = 0; i < elementsFusionRates.size(); i++) {
		totalNetworkFusionRate += sumDoublesInVector(elementsFusionRates[i]);
	}
}
void OrganelleNetwork::getTotalNetworkFusionRate(double &outputDouble){
	outputDouble = totalNetworkFusionRate;
}

// EXIT rates of every element.
void OrganelleNetwork::updateElementsExitRates(){
	elementsExitRates.assign(networkOrganisation.size(), 0);
	for (size_t i = 0; i < networkOrganisation.size(); i++) {
		networkOrganisation[i].updateElementRates();
		networkOrganisation[i].getTotalElementExitRate(elementsExitRates[i]);
	}
}
void OrganelleNetwork::getElementsExitRates(std::vector<double>& outputDoubleVector) const{
	outputDoubleVector = elementsExitRates;
}
void OrganelleNetwork::updateTotalNetworkExitRate(){
	totalNetworkExitRate = sumDoublesInVector(elementsExitRates);
}
void OrganelleNetwork::getTotalNetworkExitRate(double &outputDouble) const{
	outputDouble = totalNetworkExitRate;
}

// BUDDING rates of every element.
void OrganelleNetwork::updateElementsBuddingRates(){
	elementsBuddingRates.assign(networkOrganisation.size(), 0);
	for (size_t i = 0; i < networkOrganisation.size(); i++) {
		networkOrganisation[i].updateElementRates();
		networkOrganisation[i].getElementBuddingRate(elementsBuddingRates[i]);
	}
}
void OrganelleNetwork::getElementsBuddingRates(std::vector<double>& outputDoubleVector) const{
	outputDoubleVector = elementsBuddingRates;
}
void OrganelleNetwork::updateTotalNetworkBuddingRate(){
	totalNetworkBuddingRate = sumDoublesInVector(elementsBuddingRates);
}
void OrganelleNetwork::getTotalNetworkBuddingRate(double &outputDouble) const{
	outputDouble = totalNetworkBuddingRate;
}

// MATURATION rates of every element.
void OrganelleNetwork::updateElementsMaturationRates(){
	elementsMaturationRates.assign(networkOrganisation.size(), 0);
	for (size_t i = 0; i < networkOrganisation.size(); i++) {
		networkOrganisation[i].updateElementRates();
		networkOrganisation[i].getTotalElementMaturationRate(elementsMaturationRates[i]);
	}
}
void OrganelleNetwork::getElementsMaturationRates(std::vector<double>& outputDoubleVector) const{
	outputDoubleVector = elementsMaturationRates;
}
void OrganelleNetwork::updateTotalNetworkMaturationRate(){
	totalNetworkMaturationRate = sumDoublesInVector(elementsMaturationRates);
}
void OrganelleNetwork::getTotalNetworkMaturationRate(double &outputDouble) const{
	outputDouble = totalNetworkMaturationRate;
}

// GLOBAL : every total network rates.
void OrganelleNetwork::updateGlobalRates(){
	globalNetworkRates[0] = influxRate;
	globalNetworkRates[1] = totalNetworkFusionRate;
	globalNetworkRates[2] = totalNetworkExitRate;
	globalNetworkRates[3] = totalNetworkBuddingRate;
	globalNetworkRates[4] = totalNetworkMaturationRate;
}
void OrganelleNetwork::getGlobalRates(std::vector<double>& outputDoubleVector) const{
	outputDoubleVector = globalNetworkRates;
}

// MASSIVER UPDATER VON HÖLLE!!
void OrganelleNetwork::updateElementRates(int &index){
	networkOrganisation[index].updateElementRates();
	networkOrganisation[index].getElementSize(elementsSizes[index]);
	networkOrganisation[index].getElementIdentity(elementsIdentities[index]);
	networkOrganisation[index].getElementCargoes(cargoesOrganization[index]);
	networkOrganisation[index].getFollowedCargoes(followedCargoOrganization[index]);
	networkOrganisation[index].getElementPurity(elementsPurities[index]);
	updateTotalNetworkPurity();
	updateMatrix(index);
	updateTotalNetworkFusionRate();
	networkOrganisation[index].getTotalElementExitRate(elementsExitRates[index]);
	updateTotalNetworkExitRate();
	networkOrganisation[index].getElementBuddingRate(elementsBuddingRates[index]);
	updateTotalNetworkBuddingRate();
	networkOrganisation[index].getTotalElementMaturationRate(elementsMaturationRates[index]);
	updateTotalNetworkMaturationRate();
}

// update the fusion matrix when an element is modified:
void OrganelleNetwork::updateMatrix(int &indexOfElement){
	doubleVectorToolOneFunc.assign(elementsIdentities[indexOfElement].begin(), elementsIdentities[indexOfElement].end());
	for (int j = 0; j < indexOfElement; j++) {
		doubleVectorToolTwoFunc.assign(elementsIdentities[j].begin(), elementsIdentities[j].end());
		elementsFusionRates[indexOfElement][j] = homology(doubleVectorToolOneFunc, doubleVectorToolTwoFunc);
	}
	for (size_t i = indexOfElement + 1; i < elementsFusionRates.size(); i++) {
		doubleVectorToolTwoFunc.assign(elementsIdentities[i].begin(), elementsIdentities[i].end());
		elementsFusionRates[i][indexOfElement] = homology(doubleVectorToolOneFunc, doubleVectorToolTwoFunc);
	}
}

				/****************************/
				/*			EVENTS			*/
				/****************************/

// INJECTION of a new element (vesicle) from the entry.
void OrganelleNetwork::vesicleInjection(){
	firstIndex = networkOrganisation.size();
	networkOrganisation.push_back(element::OrganelleElement(buddingRate, maturationRates,
			exitIdentities,membraneCargoesAffinity, buddingCargoesAffinity));
	elementsSizes.push_back(0);
	elementsIdentities.push_back({});
	elementsFusionRates.push_back(std::vector<double> (elementsFusionRates.size(), 0));
	elementsExitRates.push_back(0);
	elementsBuddingRates.push_back(0);
	elementsMaturationRates.push_back(0);
	elementsPurities.push_back(0);
	cargoesOrganization.push_back({});
	followedCargoOrganization.push_back({});
	updateElementRates(firstIndex);
	updateGlobalRates();
}

// FUSION of two elements. (random number to choose the elements that fuse)
int OrganelleNetwork::elementsFusion(double inputDouble, double inputDouble2){
	// We project the fusion rates matrix:
	chooseIndexMatrix(inputDouble, inputDouble2, elementsFusionRates, doubleVectorToolOne, doubleToolOne, firstIndex, secondIndex);
	// We verify that the two elements are different:
	if (firstIndex >= secondIndex){
		std::cout << "ERROR! I cannot fuse an element with itself" << std::endl;
	}
	// We modify the objects:
	networkOrganisation[firstIndex].getElementIdentity(intVectorToolOne);
	networkOrganisation[secondIndex].getElementIdentity(intVectorToolTwo);
	networkOrganisation[firstIndex].setElementIdentity(sumIntVectors(intVectorToolOne, intVectorToolTwo));
	networkOrganisation[firstIndex].getElementSize(intToolOne);
	networkOrganisation[secondIndex].getElementSize(intToolTwo);
	networkOrganisation[firstIndex].setElementSize(intToolOne + intToolTwo);
	networkOrganisation[firstIndex].getElementCargoes(intVectorToolOne);
	networkOrganisation[secondIndex].getElementCargoes(intVectorToolTwo);
	networkOrganisation[firstIndex].setElementCargoes(sumIntVectors(intVectorToolOne, intVectorToolTwo));
	networkOrganisation[firstIndex].getFollowedCargoes(intVectorToolOne);
	networkOrganisation[secondIndex].getFollowedCargoes(intVectorToolTwo);
	networkOrganisation[firstIndex].setFollowedCargoes(sumIntVectors(intVectorToolOne, intVectorToolTwo));
	networkOrganisation.erase(networkOrganisation.begin() + secondIndex);
	networkOrganisation[firstIndex].updateElementRates();

	elementsSizes.erase(elementsSizes.begin() + secondIndex);
	elementsIdentities.erase(elementsIdentities.begin() + secondIndex);
	cargoesOrganization.erase(cargoesOrganization.begin() + secondIndex);
	followedCargoOrganization.erase(followedCargoOrganization.begin() + secondIndex);
	for (size_t i = secondIndex + 1; i < elementsFusionRates.size(); i++) {
		elementsFusionRates[i].erase(elementsFusionRates[i].begin() + secondIndex);
	}
	elementsFusionRates.erase(elementsFusionRates.begin() + secondIndex);
	elementsExitRates.erase(elementsExitRates.begin() + secondIndex);
	elementsBuddingRates.erase(elementsBuddingRates.begin() + secondIndex);
	elementsMaturationRates.erase(elementsMaturationRates.begin() + secondIndex);
	elementsPurities.erase(elementsPurities.begin() + secondIndex);

	updateElementRates(firstIndex);
	updateGlobalRates();
	if (sumIntegersInVector(intVectorToolOne) != 0 || sumIntegersInVector(intVectorToolTwo) != 0) {
		return firstIndex;
	} else {
		return -1;
	}
}

// EXIT of an element by fusing with the entry or an exit. (random number to
// choose the element that exit)
std::vector<std::vector<int> > OrganelleNetwork::elementExit(double firstInputDouble, double secondInputDouble, std::string path, double currentTime){
	// Element that exits:
	firstIndex = chooseIndex(firstInputDouble, elementsExitRates);
	// the chosen exit:
	networkOrganisation[firstIndex].getElementExitRates(doubleVectorToolOne);
	secondIndex = chooseIndex(secondInputDouble, doubleVectorToolOne);
	if (doWeFollowExit == 1) {
		followExit(firstIndex, secondIndex, path, currentTime);
	}
	networkOrganisation[firstIndex].getFollowedCargoes(intVectorToolOne);

	networkOrganisation.erase(networkOrganisation.begin() + firstIndex);
	elementsSizes.erase(elementsSizes.begin() + firstIndex);
	elementsIdentities.erase(elementsIdentities.begin() + firstIndex);
	cargoesOrganization.erase(cargoesOrganization.begin() + firstIndex);
	followedCargoOrganization.erase(followedCargoOrganization.begin() + firstIndex);
	for (size_t i = firstIndex + 1; i < elementsFusionRates.size(); i++) {
		elementsFusionRates[i].erase(elementsFusionRates[i].begin() + firstIndex);
	}
	elementsFusionRates.erase(elementsFusionRates.begin() + firstIndex);
	elementsExitRates.erase(elementsExitRates.begin() + firstIndex);
	elementsBuddingRates.erase(elementsBuddingRates.begin() + firstIndex);
	elementsMaturationRates.erase(elementsMaturationRates.begin() + firstIndex);
	elementsPurities.erase(elementsPurities.begin() + firstIndex);

	updateTotalNetworkPurity();
	updateTotalNetworkFusionRate();
	updateTotalNetworkExitRate();
	updateTotalNetworkBuddingRate();
	updateTotalNetworkMaturationRate();

	updateGlobalRates();

	if (sumIntegersInVector(intVectorToolOne) != 0) {
		return {{secondIndex},{intVectorToolOne}};
	} else {
		return {{-1}};
	}
}

// BUD a vesicle from an element (2 random numbers to choose the element
// from which it buds and the identity contained in the new vesicle)
int OrganelleNetwork::vesicleBudding(double firstInputDouble, double secondInputDouble, std::vector<double> inputDoubleVectorOne, std::vector<double> inputDoubleVectorTwo){
	// We modify the organization:
	firstIndex = chooseIndex(firstInputDouble, elementsBuddingRates);
	secondIndex = networkOrganisation.size();
	networkOrganisation.push_back(element::OrganelleElement(buddingRate, maturationRates,
			exitIdentities,membraneCargoesAffinity, buddingCargoesAffinity));

	// We bud:
	intMatrixToolFunc = networkOrganisation[firstIndex].budVesicle(secondInputDouble,inputDoubleVectorOne, inputDoubleVectorTwo);

	// We generate the correct vesicle:
	intVectorToolOne.assign(maturationRates.size() + 1, 0);
	intVectorToolOne[intMatrixToolFunc[0][0]]++;
	networkOrganisation.back().setElementIdentity(intVectorToolOne);
	networkOrganisation.back().setElementCargoes(intMatrixToolFunc[1]);
	networkOrganisation.back().setFollowedCargoes(intMatrixToolFunc[2]);
	networkOrganisation[firstIndex].updateElementRates();
	networkOrganisation.back().updateElementRates();

	// We initialize the followers:
	elementsSizes.push_back(0);
	elementsIdentities.push_back({});
	elementsFusionRates.push_back(std::vector<double> (elementsFusionRates.size(), 0));
	elementsExitRates.push_back(0);
	elementsBuddingRates.push_back(0);
	elementsMaturationRates.push_back(0);
	elementsPurities.push_back(0);
	cargoesOrganization.push_back({});
	followedCargoOrganization.push_back({});

	// We need to get the second identity before:
	networkOrganisation.back().getElementIdentity(elementsIdentities.back());
	updateElementRates(firstIndex);
	updateElementRates(secondIndex);
	updateGlobalRates();

	if (sumIntegersInVector(intMatrixToolFunc[2]) != 0) {
		return secondIndex;
	} else {
		return -1;
	}
}

// MATURE a specy inside an element (2 random numbers, one to choose the
// element, one to chose the specy that mature)
int OrganelleNetwork::specyMaturation(double firstInputDouble, double secondInputDouble){
	firstIndex = chooseIndex(firstInputDouble, elementsMaturationRates);
	networkOrganisation[firstIndex].matureElement(secondInputDouble);
	networkOrganisation[firstIndex].updateElementRates();
	updateElementRates(firstIndex);
	updateGlobalRates();
	if (sumIntegersInVector(followedCargoOrganization[firstIndex]) != 0) {
		return firstIndex;
	} else {
		return -1;
	}
}

// CARGO injection:
int OrganelleNetwork::injectCargo(std::vector<int> firstIntVector, std::vector<int> secondIntVector){
	firstIndex = networkOrganisation.size();
	vesicleInjection();
	networkOrganisation.back().setElementCargoes(firstIntVector);
	networkOrganisation.back().setFollowedCargoes(secondIntVector);
	networkOrganisation.back().updateElementRates();
	updateElementRates(firstIndex);
	updateGlobalRates();
	if (sumIntegersInVector(secondIntVector) != 0) {
		return firstIndex;
	} else {
		return -1;
	}
}

				/********************************/
				/*			FOLLOWERSS			*/
				/********************************/

// snapshot the total organization:
void OrganelleNetwork::snapshotOrganisation(std::string path, double currentTime){
	std::ofstream outputFile;
	outputFile.open(path.c_str(), std::ios::out | std::ios::app);
	outputFile << currentTime << " || ";
	for (size_t i = 0; i < elementsIdentities.size(); i++) {
		for (size_t j = 0; j < elementsIdentities[i].size(); ++j) {
			outputFile << elementsIdentities[i][j] << " ; ";
		}
		outputFile << " || ";
	}
	outputFile << "\n";
	outputFile.close();
}

// follow the size distribution:
void OrganelleNetwork::followSizes(std::string path, double currentTime){
	std::ofstream outputFile;
	outputFile.open(path.c_str(), std::ios::out | std::ios::app);
	outputFile << currentTime << " ; ";
	for (size_t i = 0; i < elementsSizes.size(); i++) {
		outputFile << elementsSizes[i] << " ; ";
	}
	outputFile << "\n";
	outputFile.close();
}

// follow the network rates:
void OrganelleNetwork::followRates(std::string path, double currentTime){
	std::ofstream outputFile;
	outputFile.open(path.c_str(), std::ios::out | std::ios::app);
	outputFile << currentTime << " ; ";
	for (size_t i = 0; i < globalNetworkRates.size(); i++) {
		outputFile << globalNetworkRates[i] << " ; ";
	}
	outputFile << totalNetworkPurity << "\n";
	outputFile.close();
}

// follow the global network identities:
void OrganelleNetwork::followIdentity(std::string path, double currentTime){
	std::ofstream outputFile;
	outputFile.open(path.c_str(), std::ios::out | std::ios::app);
	outputFile << currentTime << " ; ";
	intVectorToolOneFunc.assign(maturationRates.size() + 1, 0);
	for (size_t i = 0; i < elementsIdentities.size(); i++) {
		for (size_t j = 0; j < maturationRates.size() + 1; ++j) {
			intVectorToolOneFunc[j] += elementsIdentities[i][j];
		}
	}
	for (size_t i = 0; i < intVectorToolOneFunc.size(); ++i) {
		outputFile << intVectorToolOneFunc[i] << " ; ";
	}
	outputFile << totalNetworkPurity << "\n";
	outputFile.close();
}

// follow the exit identities:
void OrganelleNetwork::followExit(int indexToExit, int indexOfExit, std::string path, double currentTime){
	std::ofstream outputFile;
	outputFile.open(path.c_str(), std::ios::out | std::ios::app);
	outputFile << currentTime << " ; " << indexOfExit << " ; ";
	for (size_t i = 0; i < maturationRates.size(); ++i) {
		outputFile << elementsIdentities[indexToExit][i] << " | ";
	}
	outputFile << elementsIdentities[indexToExit].back() << "\n";
	outputFile.close();
}
// follow the cargoes:
void OrganelleNetwork::followCargoes(std::string path, double currentTime){
	std::ofstream outputFile;
	outputFile.open(path.c_str(), std::ios::out | std::ios::app);
	outputFile << currentTime << " ; ";
	for (size_t i = 0; i < cargoesOrganization.size(); i++) {
		for (size_t j = 0; j < cargoesOrganization[i].size() - 1; ++j) {
			outputFile << cargoesOrganization[i][j] << " | ";
		}
		outputFile << cargoesOrganization[i].back() << " ; ";
	}
	outputFile << "\n";
	outputFile.close();
}

// follow the chosen one:
void OrganelleNetwork::followCargoLife(int event, int index, std::string path, double currentTime){
	std::ofstream outputFile;
	outputFile.open(path.c_str(), std::ios::out | std::ios::app);
	outputFile << currentTime << " ; " << event << " ; ";
	for (size_t i = 0; i < maturationRates.size(); ++i) {
		outputFile << elementsIdentities[index][i] << " | ";
	}
	outputFile << elementsIdentities[index].back() << " ; ";
	for (size_t i = 0; i < followedCargoOrganization[index].size() - 1; ++i) {
		outputFile << followedCargoOrganization[index][i] << " | ";
	}
	outputFile << followedCargoOrganization[index].back();
	outputFile << "\n";
	outputFile.close();
}

} /* namespace network */
