/*
 * OrganelleElement.hpp
 *
 *  Created on: Apr 12, 2017
 *      Author: jean-patrick
 */

#ifndef ORGANELLEELEMENT_HPP_
#define ORGANELLEELEMENT_HPP_

#include <iostream>
#include <vector>

namespace element {

class OrganelleElement {
private:
	// To build the object:
	/*
	 * BUDDING RATE:
	 * Budding rate per unit area
	 */
	double buddingRate;
	/*
	 * MATURATION RATES:
	 * Rates to mature an identity unit from one biochemical identity to another.
	 * Stored in a vector.
	 */
	std::vector<double> maturationRates;
	/*
	 * EXIT COMPARTMENT IDENTITIES:
	 * The fraction of every identity contained in an exit compartment is stored
	 * in a matrix with a number of lines defined by the number of exits (ER,
	 * exit 1, exit 2, etc...), and a number of columns defined by the number of
	 * possible identities in an element (A, B, C ...).
	 */
	std::vector<std::vector<double> > exitIdentities;
	/*
	 * CARGOES AFFINITY:
	 * Two vectors, one for the membrane identity, the other for the budding
	 * membrane.
	 */
	std::vector<std::vector<double> > membraneCargoesAffinity;
	std::vector<std::vector<double> > buddingCargoesAffinity;

	// That can change during the simulation:
	/*
	 * SIZE:
	 * The total number of fused vesicles inside the element.
	 */
	int elementSize;
	/*
	 * ELEMENT IDENTITY:
	 * Stored in a vector. Every biochemical identity has a place in the vector
	 * which is associated to an integer : the number of unit of this identity
	 * inside the element. Its size is larger than the maturation rates vector
	 * size. Number of identities helps to choose which identity buds.
	 */
	std::vector<int> elementIdentity;
	double numberOfIdentities;
	/*
	 * EXIT RATE:
	 * A vector containing  the flux of probability for an element to exit
	 * by fusing with every exit.
	 */
	std::vector<double> elementExitRates;
	double totalElementExitRate;
	/*
	 * ELEMENT BUDDING RATE:
	 * A double where the budding rate of a new vesicle is stored.
	 */
	double elementBuddingRate;
	/*
	 * ELEMENT MATURATION RATE:
	 * Rate of maturation for the element for every specie.
	 */
	std::vector<double> elementMaturationRates;
	double totalElementMaturationRate;
	/*
	 * ELEMENT PURITY:
	 * Equals to 1 if the element is pure, and 0 if it is perfectly mixed.
	 */
	double elementPurity;
	/*
	 * CARGOES:
	 * Stored in a vector of integers, one for every type of cargoes.
	 * We store also a maximum of 1 cargo in a second vector, to follow this guy.
	 */
	std::vector<int> elementCargoes;
	std::vector<int> followedCargo;
	/*
	 * USEFUL TOOLS:
	 * to prevent further declarations.
	 */
	std::vector<int> intVectorToolOne;
	std::vector<int> intVectorToolTwo;
	std::vector<double> doubleVectorToolOne;
	std::vector<double> doubleVectorToolTwo;
	std::vector<double> doubleVectorToolFunc;
	int intToolOne;
	int intToolTwo;
	double doubleToolOne;
	double doubleToolTwo;
	std::vector<std::vector<int> > intmatrixTool;

public:
	//CONSTRUCTOR, DESTRUCTOR:
	OrganelleElement(double buddingRateSimu,
			std::vector<double> maturationRatesSimu,
			std::vector<std::vector<double> > exitIdentitiesSimu,
			std::vector<std::vector<double> > membraneCargoesAffinitySimu,
			std::vector<std::vector<double> > buddingCargoesAffinitySimu
			);
	virtual ~OrganelleElement();

	// SETTER AND GETTERS, and updaters, in lowcaps, because they deserve it.
	// SIZE:
	void setElementSize(int inputInt);
	void getElementSize(int &outputInt) const;
	// ELEMENT IDENTITY:
	void setElementIdentity(std::vector<int> inputIntVector);
	void getElementIdentity(std::vector<int>& outputIntVector) const;
	void updateNumberOfIdentities();
	// EXIT RATE:
	void updateElementExitRates();
	void getElementExitRates(std::vector<double>& outputVector) const;
	void updateTotalElementExitRate();
	void getTotalElementExitRate(double &outputDouble) const;
	// ELEMENT BUDDING RATE:
	void updateElementBuddingRate();
	void getElementBuddingRate(double &outputDouble) const;
	// ELEMENT MATURATION RATE:
	void updateElementMaturationRates();
	void getElementMaturationRates(std::vector<double>& outputDoubleVector) const;
	void updateTotalElementMaturationRate();
	void getTotalElementMaturationRate(double &outputDouble) const;
	// ELEMENT PURITY:
	void updateElementPurity();
	void getElementPurity(double &outputDouble) const;
	// CARGOES:
	void setElementCargoes(std::vector<int> inputIntVector);
	void getElementCargoes(std::vector<int>& outputIntVector) const;
	void setFollowedCargoes(std::vector<int> inputIntVector);
	void getFollowedCargoes(std::vector<int>& outputIntVector) const;
	// MASSSSSSSIVE UPDATER OF HELL:
	void updateElementRates();

	// EVENTS:
	void matureElement(double inputDouble);
	std::vector<std::vector<int> > budVesicle(double inputDouble, std::vector<double> inputDoubleVectorOne, std::vector<double> inputDoubleVectorTwo);
};

} /* namespace element */

#endif /* ORGANELLEELEMENT_HPP_ */
