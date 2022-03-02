Gillespie simulation of a self-organized Golgi.

Please contact jean-patrick.vrel@curie.fr if you are interested or if you can't make this program work.

This is the very first version of this simulation. A newer, more readable and user friendly will soon be published. 

This program is free to use, share and modify (MIT Licence), as long as you're including the MIT License in your own project. If you do, please consider citing our work. 

It has been written and used on a UNIX system (GNU Linux, but it should work on Mac), and you need both g++ and the boost library to compile it.

The goal of this README is not to explain the entire functioning of the code, but to provide sufficient knowledge to use it and reproduce the results presented in the paper. More information are given in the code as comments in the header files for each function. We tried our best to make this program readable, but as it is a very preliminary version (particularly for the cargoes) the implementation may be complicated to understand. If you're facing any difficulties, feel free to contact us.


#########################


Quick start:
- Create a folder and put all the .hpp and .cpp files in it, except program.cpp
- Create a second folder, in the first one, with program.cpp and the makefile. This is your simulation folder.
- Create a "Data" folder (by default, can be changed in program.cpp) in which all the simulation will be performed.
It should look like this:
 Folder1
    | Functions.hpp
    | Functions.cpp
    | OrganelleElement.hpp
    | OrganelleElement.cpp
    | OrganelleNetwork.hpp
    | OrganelleNetwork.cpp
    | TimeFlow.hpp
    | TimeFlow.cpp
    |SimulationFolder
        | makefile
        | program.cpp
        |Data
            |
            
- Then open program.cpp and set the following elements as you want:
    - line 20, set the maturation rate
    - line 21, the budding rate
    - line 23, the identities for the ER and the TGN (by default 100% Cis, {1,0,0} for the ER and 100% Trans {0,0,1} for the TGN}
    - line 26, the affinity for the compartment membrane of ALL cargoes you want to follow (by default, 4 that don't have any preference for any identity)
    - line 27, the affinity for the vesicle membrane of ALL cargoes you want to follow (by default, 2 that can't bud and 2 that don't have any preference for any identity)
    - Line 37, you can change the name of the parent directory in which simulation results are stored, make sure it corresponds the folder name you created to store data, and the ER/TGN identities
    - Line 51, the number of iterations you want to wait before starting snapshots of the system (to reach a steady-state)
    - line 62, the number of iterations you want to acquire
    - line 66, how many iterations between 2 snapshots (does not affect how you follow cargoes). If you want to wait 25 points between 2 snapshots, it'll be: if (a%25 == 0) {
    
- Open a terminal and type "make" and enter.
- then type "./program" and enter. It will run, and return the total counts of each mechanism it applies on the system (injection, fusion, exit, budding, maturation).

All data for this particular simulation will be stored in the Data/ folder. You will find, by default, sub-folders with the parameters. In these, you'll find folders with the date of simulations. They contain 4 files by default: a parameters file, snapshots of the systems and its rates (including its purity), and the cargoes we followed. To treat these data, please read the README file provided with the post-treatment scripts. 


#########################

Multiple simulations in one compilation:
You can increment the maturation and budding rate (set lines 20 and 21) by tuning the respective maximums for i (line 32) and j (line 31), and their increase factor line 80 and 82.
Then compile using make, and launch using "./program" as before.
