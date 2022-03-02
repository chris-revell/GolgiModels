/*
 * OrganelleNetwork.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: jean-patrick
 */

#ifndef ORGANELLENETWORK_HPP_
#define ORGANELLENETWORK_HPP_

#include "OrganelleElement.hpp"

namespace network {

class OrganelleNetwork {

private:
	// To build the object:
	/*
	 * INFLUX:
	 * The probability flux for the vesicle injection.
	 */
	double influxRate;
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
	int doWeFollowExit; // 0: nope, 1: yep!

	// That change during the simulation:
	// The followers:
	/*
	 * ORGANISATION:
	 * Every organelle element is stored in a vector.
	 */
	std::vector<element::OrganelleElement> networkOrganisation;
	/*
	 * SIZE:
	 * The size of every element is stored in a vector at the same index than the
	 * corresponding element.
	 */
	std::vector<int> elementsSizes;
	/*
	 * IDENTITIES:
	 * A vector with every identities of every element, to follow all this.
	 */
	std::vector<std::vector<int> > elementsIdentities;
	/*
	 * CARGOES:
	 * A vector where cargoes are stored to know where there are.
	 */
	std::vector<std::vector<int> > cargoesOrganization;
	std::vector<std::vector<int> > followedCargoOrganization;
	/*
	 * PURITY:
	 * The purity of every element is stored in a double vector. We also calculate
	 * the total purity, normalized by the size.
	 */
	std::vector<double> elementsPurities;
	double totalNetworkPurity;

	// The rates:
	/*
	 * FUSION:
	 * Every fusion event between two compartments is stored in a matrix. We also
	 * compute the sum of all rates.
	 */
	std::vector<std::vector<double> > elementsFusionRates;
	double totalNetworkFusionRate;
	/*
	 * EXIT:
	 * The exit rate of every element is stored in a vector. We also calculate
	 * the total rate.
	 */
	std::vector<double> elementsExitRates;
	double totalNetworkExitRate;
	/*
	 * BUDDING:
	 * The budding rate of every element is stored in a vector. We also calculate
	 * the total rate.
	 */
	std::vector<double> elementsBuddingRates;
	double totalNetworkBuddingRate;
	/*
	 * MATURATION:
	 * The total maturation rate for every element is stored in a vector. We also
	 * calculate the total rate of maturation.
	 */
	std::vector<double> elementsMaturationRates;
	double totalNetworkMaturationRate;
	/*
	 * GLOBAL:
	 * Every total network rate is stored in a vector. At index 0 we store the
	 * injection rate, at 1 the fusion rate, at 2 the exit rate, at 3 the
	 * budding rate, at 4 the maturation rate.
	 */
	std::vector<double> globalNetworkRates;

	/*
	 * USEFUL TOOLS:
	 * To prevent further declarations.
	 */
	std::vector<int> intVectorToolOne;
	std::vector<int> intVectorToolTwo;
	std::vector<double> doubleVectorToolOne;
	std::vector<double> doubleVectorToolTwo;
	double doubleToolOne;
	double doubleToolTwo;
	std::vector<std::vector<double> > doubleMatrixTool;
	int firstIndex;
	int secondIndex;
	int intToolOne;
	int intToolTwo;
	std::vector<int> intVectorToolOneFunc;
	std::vector<int> intVectorToolTwoFunc;
	std::vector<double> doubleVectorToolOneFunc;
	std::vector<double> doubleVectorToolTwoFunc;
	double doubleToolOneFunc;
	double doubleToolTwoFunc;
	std::vector<std::vector<int> > intMatrixToolFunc;
	std::vector<std::vector<double> > doubleMatrixToolFunc;
	int firstIndexFunc;
	int secondIndexFunc;
	int intToolOneFunc;
	int intToolTwoFunc;


public:
	OrganelleNetwork(double influxRateSimu,
			double buddingRateSimu,
			std::vector<double> maturationRatesSimu,
			std::vector<std::vector<double> > exitIdentitiesSimu,
			std::vector<std::vector<double> > membraneCargoesAffinitySimu,
			std::vector<std::vector<double> > buddingCargoesAffinitySimu,
			int doWeFollowExitSimu
			);
	virtual ~OrganelleNetwork();

	// METHODS:
					/********************************/
					/*		SETTERS AND GETTERS		*/
					/********************************/


	// SIZE updater (if needed) and getter.
	void updateElementsSize();//if needed
	void getElementsSize(std::vector<int>& outputIntVector) const;

	// IDENTITIES getter and updater:
	void updateElementsIdentities();
	void getElementsIdentities(std::vector<std::vector<int> >& outputIntMatrix) const;

	// CARGOES getter and setter:
	void updateCargoesOrganization();
	void getCargoesOrganization(std::vector<std::vector<int> >& outputIntMatrix) const;
	void updateFollowedCargoesOrganization();
	void getFollowedCargoesOrganization(std::vector<std::vector<int> >& outputIntMatrix) const;

	// PURITY of every element:
	void updateElementsPurities(); //if needed
	void getElementsPurities(std::vector<double>& outputDoubleVector) const;
	void updateTotalNetworkPurity();
	void getTotalNetworPurity(double &outputDouble) const;

	// FUSION rates between two elements.
	void updateElementsFusionRates(); //if needed
	void getElementsFusionRates(std::vector<std::vector<double> >& outputMatrix) const;
	void updateTotalNetworkFusionRate();
	void getTotalNetworkFusionRate(double &outputDouble);

	// EXIT rates of every element.
	void updateElementsExitRates(); //if needed
	void getElementsExitRates(std::vector<double>& outputDoubleVector) const;
	void updateTotalNetworkExitRate();
	void getTotalNetworkExitRate(double &outputDouble) const;

	// BUDDING rates of every element.
	void updateElementsBuddingRates(); //if needed
	void getElementsBuddingRates(std::vector<double>& outputDoubleVector) const;
	void updateTotalNetworkBuddingRate();
	void getTotalNetworkBuddingRate(double &outputDouble) const;

	// MATURATION rates of every element.
	void updateElementsMaturationRates(); //if needed
	void getElementsMaturationRates(std::vector<double>& outputDoubleVector) const;
	void updateTotalNetworkMaturationRate();
	void getTotalNetworkMaturationRate(double &outputDouble) const;

	// GLOBAL : every total network rates.
	void updateGlobalRates();
	void getGlobalRates(std::vector<double>& outputDoubleVector) const;

	// MASSIVER UPDATER VON HÖLLE!!
	void updateElementRates(int &index);

	// update the fusion matrix when an element is modified:
	void updateMatrix(int &indexOfElement);

					/****************************/
					/*			EVENTS			*/
					/****************************/

	// INJECTION of a new element (vesicle) from the entry.
	void vesicleInjection();

	// FUSION of two elements. (random number to choose the elements that fuse)
	int elementsFusion(double inputDouble, double inputDouble2);

	// EXIT of an element by fusing with the entry or an exit. (random number to
	// choose the element that exit)
	std::vector<std::vector<int> > elementExit(double firstInputDouble, double secondInputDouble, std::string path, double currentTime);

	// BUD a vesicle from an element (2 random numbers to choose the element
	// from which it buds and the identity contained in the new vesicle)
	int vesicleBudding(double firstInputDouble, double secondInputDouble, std::vector<double> inputDoubleVectorOne, std::vector<double> inputDoubleVectorTwo);

	// MATURE a specy inside an element (2 random numbers, one to choose the
	// element, one to chose the specy that mature)
	int specyMaturation(double firstInputDouble, double secondInputDouble);

	// CARGO injection:
	int injectCargo(std::vector<int> firstIntVector, std::vector<int> secondIntVector);

					/********************************/
					/*			FOLLOWERSS			*/
					/********************************/

	// snapshot the total organization:
	void snapshotOrganisation(std::string path, double currentTime);

	// follow the size distribution:
	void followSizes(std::string path, double currentTime);

	// follow the network rates:
	void followRates(std::string path, double currentTime);

	// follow the global network identities:
	void followIdentity(std::string path, double currentTime);

	// follow the exit identities and cargoes:
	void followExit(int indexToExit, int indexOfExit, std::string path, double currentTime);

	// follow the cargoes:
	void followCargoes(std::string path, double currentTime);

	// follow the chosen one:
	void followCargoLife(int event, int index, std::string path, double currentTime);

};

} /* namespace network */

#endif /* ORGANELLENETWORK_HPP_ */
