/*
 * TimeFlow.hpp
 *
 *  Created on: Apr 18, 2017
 *      Author: jean-patrick
 */

#ifndef TIMEFLOW_HPP_
#define TIMEFLOW_HPP_

#include "OrganelleNetwork.hpp"

#include <iostream>
#include <vector>
#include <random>

namespace timing {

class TimeFlow {
private:
	// To build the object:
	/*
	 * PATH:
	 * Where the simulation is saved.
	 */
	std::string simulationPath;
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
	int doWeInjectCargoes;
	int doWeFollowCargoes;
	std::vector<int> cargoesComposition;

	// That change during the simulation:
	/*
	 * RANDOM NUMBER GENERATOR:
	 * We store the see, the generator, and the uniform distribution generator.
	 */
	double seed;
	std::mt19937 mt_rand;
	std::uniform_real_distribution<double> disint;
	/*
	 * ORGANELLE NETWORK:
	 * as an object in an object, objectception stuff, yep, it's a dream.
	 */
	network::OrganelleNetwork organelleNetwork;
	/*
	 * NEXT EVENT:
	 * Stored in an integer : 0 for an injection, 1 fusion, 2 exit, 3 budding.
	 */
	int nextEvent;
	/*
	 * GLOBAL:
	 * Every total network rate is stored in a vector. At index 0 we store the
	 * injection rate, at 1 the fusion rate, at 2 the exit rate, at 3 the
	 * budding rate, at 4 the maturation rate.
	 */
	std::vector<double> globalNetworkRates;
		/*
	 * TIME:
	 * the current time of the simulation and the number of time points.
	 */
	double currentTime;
	int numberOfTimePoints;
	std::vector<int> eventsCounting;
	/*
	 * FOLLOWED CARGOES:
	 * 1 for yes, 0 for no.
	 */
	std::vector<int> isThereACargo;
	/*
	 * ACCUMULATED SIZE:
	 * To follow if we are in a transient state.
	 */
	int accumulatedSize;
	/*
	 * USEFUL TOOLS:
	 * to prevent further declarations.
	 */
	int intToolOne;
	std::vector<int> intVectorToolOne;
	std::vector<int> intVectorToolTwo;
	std::vector<double> doubleVectorToolOne;
	std::vector<double> doubleVectorToolTwo;
	std::vector<std::vector<int> > intMatrixToolOne;

public:
	TimeFlow(std::string parentDirectory,
			double influxRateSimu,
			double buddingRateSimu,
			std::vector<double> maturationRatesSimu,
			std::vector<std::vector<double> > exitIdentitiesSimu,
			std::vector<std::vector<double> > membraneCargoesAffinitySimu,
			std::vector<std::vector<double> > buddingCargoesAffinitySimu,
			int doWeFollowExitSimu,
			std::vector<int> cargoesCompositionSimu
			);
	virtual ~TimeFlow();


	// METHODS:

				/********************************/
				/*		SETTERS AND GETTERS		*/
				/********************************/

	// TIME following:
	void getNumberOfTimePoints(int &outputInt) const;
	void getCurrentTime(double &outputDouble) const;

	// EVENTS counting:
	void getCountingOfEvent(std::vector<int>& outputIntVector) const;

	// CARGOES following:
	void setDoWeFollowCargoes(int inputInt);
	void setDoWeInjectCargoes(int inputInt);

				/****************************/
				/*			EVENTS			*/
				/****************************/

	// Time updating using a Gillespie algorithm
	int updateTime();

	// Choose what event occures at this time point
	void chooseTypeOfEvent();
	void setTypeOfEvent(int inputInt);
	void getTypeOfEvent(int &outputInt) const;
	void applyModification();

	// Accumulated size promoter getter and reseter:
	void accumulateSize();
	void getAccumulatedSize(int &outputInt) const;
	void resetAccumulatedSize();

	// Follow parameters during the simulation:
	// snapshot the total organization:
	void initializeSnapshotOrganisation();
	void snapshotOrganisation();
	// follow the size distribution:
	void initializeFollowSizes();
	void followSizes();
	// follow the network rates:
	void initializeFollowRates();
	void followRates();
	// follow the global network identities:
	void initializeFollowIdentity();
	void followIdentity();
	// follow the exit composition:
	void initializeFollowExit();
	// follow the cargoes:
	void initializeFollowCargoLife();
	void initializeFollowCargoes();
	void followCargoes();
};

} /* namespace timing */

#endif /* TIMEFLOW_HPP_ */
