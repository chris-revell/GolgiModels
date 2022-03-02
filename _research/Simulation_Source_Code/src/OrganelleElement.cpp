/*
 * OrganelleElement.cpp
 *
 *  Created on: Apr 12, 2017
 *      Author: jean-patrick
 */

#include "OrganelleElement.hpp"
#include "Functions.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/binomial.hpp>

namespace element {

OrganelleElement::OrganelleElement(double buddingRateSimu,
		std::vector<double> maturationRatesSimu,
		std::vector<std::vector<double> > exitIdentitiesSimu,
		std::vector<std::vector<double> > membraneCargoesAffinitySimu,
		std::vector<std::vector<double> > buddingCargoesAffinitySimu
		) {

	buddingRate = buddingRateSimu;
	maturationRates = maturationRatesSimu;
	exitIdentities = exitIdentitiesSimu;
	membraneCargoesAffinity = membraneCargoesAffinitySimu;
	buddingCargoesAffinity = buddingCargoesAffinitySimu;

	elementSize = 1;
	elementIdentity.assign(maturationRates.size() + 1, 0);
	elementIdentity.front()++;
	numberOfIdentities = 0;
	elementExitRates.assign(exitIdentities.size(), 0);
	totalElementExitRate = 0;
	elementBuddingRate = 0;
	elementMaturationRates.assign(maturationRates.size(), 0);
	totalElementMaturationRate = 0;
	elementPurity = 0;
	elementCargoes.assign(membraneCargoesAffinity.size(), 0);
	followedCargo.assign(membraneCargoesAffinity.size(), 0);

	intVectorToolOne = {};
	intVectorToolOne = {};
	doubleVectorToolOne = {};
	doubleVectorToolTwo = {};
	doubleVectorToolFunc = {};
	intToolOne = 0;
	intToolTwo = 0;
	doubleToolOne = 0;
	doubleToolTwo = 0;
	intmatrixTool = {{}};

	updateElementRates();
}

OrganelleElement::~OrganelleElement() {}


// SETTER AND GETTERS, and updaters, in lowcaps, because they deserve it.
// SIZE:
void OrganelleElement::setElementSize(int inputInt){
	elementSize = inputInt;
}
void OrganelleElement::getElementSize(int &outputInt) const{
	outputInt = elementSize;
}
// ELEMENT IDENTITY:
void OrganelleElement::setElementIdentity(std::vector<int> inputIntVector){
	elementIdentity = inputIntVector;
}
void OrganelleElement::getElementIdentity(std::vector<int>& outputIntVector) const{
	outputIntVector = elementIdentity;
}
void OrganelleElement::updateNumberOfIdentities(){
	numberOfIdentities = 0;
	for (size_t i = 0; i < elementIdentity.size(); ++i) {
		if (elementIdentity[i] > 0) {
			numberOfIdentities++;
		}
	}
}
// EXIT RATE:
void OrganelleElement::updateElementExitRates(){
	doubleVectorToolOne.assign(elementIdentity.begin(), elementIdentity.end());
	doubleVectorToolOne = normalize(doubleVectorToolOne);
	for (size_t i = 0; i < exitIdentities.size(); ++i) {
		doubleVectorToolTwo = multiplyDoublesVectors(doubleVectorToolOne, exitIdentities[i]);
		elementExitRates[i] = sumDoublesInVector(doubleVectorToolTwo);
	}
}
void OrganelleElement::getElementExitRates(std::vector<double>& outputVector) const{
	outputVector = elementExitRates;
}
void OrganelleElement::updateTotalElementExitRate(){
	totalElementExitRate = sumDoublesInVector(elementExitRates);
}
void OrganelleElement::getTotalElementExitRate(double &outputDouble) const{
	outputDouble = totalElementExitRate;
}
// ELEMENT BUDDING RATE:
void OrganelleElement::updateElementBuddingRate(){
	if (elementSize <= 1) {
		elementBuddingRate = 0;
	} else {
		doubleToolOne = elementSize;
		doubleToolTwo = numberOfIdentities;
		elementBuddingRate = doubleToolOne * doubleToolTwo * buddingRate;
	}
}
void OrganelleElement::getElementBuddingRate(double &outputDouble) const{
	outputDouble = elementBuddingRate;
}
// ELEMENT MATURATION RATE:
void OrganelleElement::updateElementMaturationRates(){
	doubleVectorToolOne.assign(elementIdentity.begin(), elementIdentity.end() - 1);
	elementMaturationRates = multiplyDoublesVectors(doubleVectorToolOne, maturationRates);
}
void OrganelleElement::getElementMaturationRates(std::vector<double>& outputDoubleVector) const{
	outputDoubleVector = elementMaturationRates;
}
void OrganelleElement::updateTotalElementMaturationRate(){
	totalElementMaturationRate = sumDoublesInVector(elementMaturationRates);
}
void OrganelleElement::getTotalElementMaturationRate(double &outputDouble) const{
	outputDouble = totalElementMaturationRate;
}
// ELEMENT PURITY:
void OrganelleElement::updateElementPurity(){
	doubleToolOne = maturationRates.size() + 1;
	if (doubleToolOne > 1) {
		elementPurity = 0;
		doubleVectorToolOne.assign(elementIdentity.begin(), elementIdentity.end());
		doubleVectorToolOne = normalize(doubleVectorToolOne);
		for (int i = 0; i < doubleToolOne; i++) {
			elementPurity += pow(doubleVectorToolOne[i] - (1.0/doubleToolOne), 2);};
		elementPurity *= doubleToolOne / (doubleToolOne - 1.0);
		elementPurity = sqrt(elementPurity);
	} else {
		elementPurity = 1;
	}

}
void OrganelleElement::getElementPurity(double &outputDouble) const{
	outputDouble = elementPurity;
}
// CARGOES:
void OrganelleElement::setElementCargoes(std::vector<int> inputIntVector){
	elementCargoes = inputIntVector;
}
void OrganelleElement::getElementCargoes(std::vector<int>& outputIntVector) const{
	outputIntVector = elementCargoes;
}
void OrganelleElement::setFollowedCargoes(std::vector<int> inputIntVector){
	followedCargo = inputIntVector;
}
void OrganelleElement::getFollowedCargoes(std::vector<int>& outputIntVector) const{
	outputIntVector = followedCargo;
}
// MASSSSSSSIVE UPDATER OF HELL:
void OrganelleElement::updateElementRates(){
	updateNumberOfIdentities();
	updateElementExitRates();
	updateTotalElementExitRate();
	updateElementBuddingRate();
	updateElementMaturationRates();
	updateTotalElementMaturationRate();
	updateElementPurity();
}

// EVENTS:
void OrganelleElement::matureElement(double inputDouble){
	intToolOne = chooseIndex(inputDouble, elementMaturationRates);
	elementIdentity[intToolOne]--;
	elementIdentity[intToolOne + 1]++;
}
std::vector<std::vector<int> > OrganelleElement::budVesicle(double inputDouble, std::vector<double> inputDoubleVectorOne, std::vector<double> inputDoubleVectorTwo){
	// We choose the specy that buds:
	doubleVectorToolOne.assign(elementIdentity.size(), 0);
	for (size_t i = 0; i < elementIdentity.size(); ++i) {
		if (elementIdentity[i] >= 1) {
			doubleVectorToolOne[i] = 1;
		}
	}
	intToolOne = chooseIndex(inputDouble, doubleVectorToolOne);
	// We generate a vector of cargoes that exit:
	intVectorToolOne.assign(membraneCargoesAffinity.size(), 0);
	for (size_t i = 0; i < membraneCargoesAffinity.size(); ++i) {
		// We iterate for every class of cargoes
		if (elementCargoes[i] > 0) {
			// We first calculate the normalization
			doubleToolOne = buddingCargoesAffinity[i][intToolOne] - membraneCargoesAffinity[i][intToolOne];
			doubleVectorToolOne.assign(elementIdentity.begin(), elementIdentity.end());
			for (size_t j = 0; j < membraneCargoesAffinity[i].size(); ++j) {
				doubleToolOne += doubleVectorToolOne[j] * membraneCargoesAffinity[i][j];
			}
			// We calculate the probability for every cargo to be in the vesicle:
			if (doubleToolOne != 0) {
				doubleToolTwo = buddingCargoesAffinity[i][intToolOne] / doubleToolOne;
			} else {
				doubleToolTwo = 1.0 / sumDoublesInVector(doubleVectorToolOne);
			}
			// We then generate the probability distribution for one, two or more cargoes in the vesicle:
			doubleVectorToolOne.assign(elementCargoes[i] + 1, 0);
			for (size_t j = 0; j < doubleVectorToolOne.size(); ++j) {
				doubleToolOne =  boost::math::binomial_coefficient<double>(elementCargoes[i], j);
				doubleVectorToolOne[j] = doubleToolOne * pow(doubleToolTwo, j) * pow(1 - doubleToolTwo, elementCargoes[i] - j);
			}
			// We select the correct number
			intVectorToolOne[i] = chooseIndex(inputDoubleVectorOne[i], doubleVectorToolOne);
		} else {
			intVectorToolOne[i] = 0;
		}
	}
	intVectorToolTwo.assign(membraneCargoesAffinity.size(), 0);
	for (size_t i = 0; i < membraneCargoesAffinity.size(); ++i) {
		// We iterate for every class of cargoes
		if (followedCargo[i] > 0) {
			// We calculate the normalization
			doubleToolOne = buddingCargoesAffinity[i][intToolOne] - membraneCargoesAffinity[i][intToolOne];
			doubleVectorToolOne.assign(elementIdentity.begin(), elementIdentity.end());
			for (size_t j = 0; j < membraneCargoesAffinity[i].size(); ++j) {
				doubleToolOne += doubleVectorToolOne[j] * membraneCargoesAffinity[i][j];
			}
			// We calculate the probability for every cargo to be in the vesicle:
			if (doubleToolOne != 0) {
				doubleToolTwo = buddingCargoesAffinity[i][intToolOne] / doubleToolOne;
			} else {
				doubleToolTwo = 1.0 / sumDoublesInVector(doubleVectorToolOne);
			}
			// We then generate the probability distribution for one, two or more cargoes in the vesicle:
			doubleVectorToolOne = {0, 0};
			doubleVectorToolOne[0] = 1 - doubleToolTwo;
			doubleVectorToolOne[1] = doubleToolTwo;
			// We select the correct number
			intVectorToolTwo[i] = chooseIndex(inputDoubleVectorTwo[i], doubleVectorToolOne);
		} else {
			intVectorToolTwo[i] = 0;
		}
	}
	// We remove the identity and the cargoes from the compartment:
	elementIdentity[intToolOne]--;
	for (size_t i = 0; i < intVectorToolOne.size(); ++i) {
		elementCargoes[i] -= intVectorToolOne[i];
		followedCargo[i] -= intVectorToolTwo[i];
	}
	elementSize --;
	return {{intToolOne}, intVectorToolOne, intVectorToolTwo};
}
} /* namespace element */
