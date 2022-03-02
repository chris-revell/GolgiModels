#include <iostream> // in order to use the std library
#include <vector>
#include <string>
#include <sstream>
#include "cmath"

#include "OrganelleElement.hpp"
#include "OrganelleNetwork.hpp"
#include "TimeFlow.hpp"
#include "Functions.hpp"

int main(int argc, char const *argv[]) {

    float time;
    clock_t t1, t2;

    using namespace timing;

    int a = 0;
    double mat = 1; // the maturation rate over the fusion rate
    double buddingRateSimu = 1; // the budding rate over the fusion rate
    std::vector<int> v;
    std::vector<std::vector<double>> exitIdentitiesSimu = {{1, 0, 0}, {0, 0, 1}}; // ER and TGN identities
    std::vector<double> maturationRatesSimu = {};
    double influxRateSimu = 0;
    std::vector<std::vector<double> > cargoesAffinitySimu = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
    std::vector<std::vector<double> > cargoesAffinitySimuBud = {{0, 0, 0}, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}};
    std::string parentDirectory = "";


    for (size_t j = 0; j < 1; j++) { // Change hier the number of maturation rates you want to investigate
    for (size_t i = 0; i < 1; i++) { // Change hier the number of budding rates you want to investigate
    influxRateSimu = 300.0 * 1.0 * (1.0 + mat) / (mat + 2.0);
    // Initialize the time to get the time of simulation and the path
    maturationRatesSimu = {mat, mat};
    parentDirectory = "Data/bud=" + std::to_string(buddingRateSimu) + "_exit=1-1_mat=" +
                      std::to_string(mat) + "_N=300/" ;
    std::cout << parentDirectory << '\n';
    TimeFlow simulation(parentDirectory,
        influxRateSimu,
        buddingRateSimu,
        maturationRatesSimu,
        exitIdentitiesSimu,
        cargoesAffinitySimu,
        cargoesAffinitySimuBud,
        0, {1, 1, 1});

    a = 0;

    for (int k = 0; k < 1000; k++) { // to reach steady-state, change the maximum value '; k < MAX_VALUE;' to set this as you like
          simulation.chooseTypeOfEvent();
          simulation.applyModification();
    }

    t1 = clock();
    simulation.initializeSnapshotOrganisation();
    simulation.initializeFollowRates();
    simulation.initializeFollowCargoLife();

    simulation.setDoWeFollowCargoes(1);
    while (a < 3000) { // simulation iterations, change the maximum value 'a < MAX_VALUE' to set this as you like
        simulation.chooseTypeOfEvent();
        simulation.updateTime();
        simulation.applyModification();
        if (a%1000 == 0) { // number of points between 2 snapshots
	         simulation.followRates();
           simulation.snapshotOrganisation();
        }
        a++;
    }
    simulation.getCountingOfEvent(v);
    printIntsVector(v);

    t2 = clock();
    time = (float)(t2-t1)/CLOCKS_PER_SEC;
    std::cout << "duration = " << time << " s" << std::endl;

    // We increase the budding rate value:
    buddingRateSimu = buddingRateSimu * pow(10.0, 1.0/10.0);
  };
  mat = mat *  pow(10.0, 1.0/5.0);
  buddingRateSimu = 0.001;
};

    return 0;
}
