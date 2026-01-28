# Notes - Last Modified 2026/01/27

## Queries

I am currently working on testing with a new molecule as a trial. For this I am starting with 1,2-Ethanediol as a new molecule to work with in the goal of testing as many functional groups as possible. These tests aren't doing any of the simulations yet, just working on generating the input. For now I am testing the script with the Build_Box format so it can generate the grofiles

I have encountered a few problems that will require fixing:

- Some approach for making sure the gro files we generate are consistent with topology files
- Need to build in some approach for how to make sure we are differentiating molecules for systems were we have a different central atom to the surrounding solvent such as in mixtures or where the Central is AA and solvent is UA etc...
- Will likely need to skip the liquid calculation if they give us a box as we won't be able to do the vacuum and pcms without the NVT_Vacuum output, we will only be able to calculate $\mu_{SCEE}$.
- Look into idea of having the option to either loop through T and P, or to zip T and P and look at specific state points.


## To Dos 
As I am writing the how it works I am coming across things I need to add or do to make it work for the system. These are easy fixes in comparison to those in Queries.
- Add an option for telling us how many molecules are in the box for Yes_Gro
- Add in steps to make sure the Top and gro are being read in as files not just strings
- Figure out the mess that is mixtures eventually

## Changes/Fixes Made

- General typos and missing : etc have been fixed
- Made made modifications to how *Gaussian_Calculations.py* and *Gro_Simulations.py* did not work previously now a series of if and elif statements depending on system. I do think this is messier than how I had it before.
- The step for locating GMX will need to be updated to work with input Not sure how to go about modifying that just yet
- How #include is read for building the oniom generation and dat generation as it isn't quite working anymore outside the big if/else sequence we had. I am thinking we could concat the file contents.

## How the Scripts Work

- The first requirement for this process is to fill information required for the molecule system within the Settings.yml, where there are two main sections requiring information.
    - State conditions is where you can provide the temperature, pressure etc. We also require some information about the molecule from experimental values such as the gas phase dipole moment, dielectric constant and refractive index. You will need to find the closest solvent keyword using https://gaussian.com/scrf/ (I think we need to find something for molecules not on this list, but I will be honest and say I have no idea how that works)
    - In this section is also where you will tell us where your software is stored so we can run our calculations on it. (Think we could add something that locates it for the user instead somehow)
    - Next you can set the mode of the system and if it is or isn't a mixture. You can select from if you have already got the coordinates for the box you wish to use (Yes_Gro), if you have the coordinates of the molecule in the vacuum and you need the script to build the box (No_Gro should change to No_Box maybe) and if you would like our script to generate the coordinates using RDKit and then generates the box using a Smile sting in AA. Each option has a different number of requirements, for instance if you have already given us a generated box of molecules, we do not need the density and molecule of the molecule at ambient conditions to calculate the number of molecules in the box, but you will need to tell us how many molecules are in the box to generate our inhouse oniom script. Similarly if you have given us the gro file for a single molecule, you will not need to give us the Smile string to build the molecule with RDKit. If it is a mixture there is space to add what the solute is in comparison to pure solvents.
    - Advance_Settings is were we store the information about the box length, cluster radius and charge scaling. We do not recommend changing these unless you know what you are doing. Increasing cluster radius will need the box length to increase. 
    - If your system requires us to build the box and is a mixture it will run through the process of calculating the dipole moment for pure solvent first and then repeat most of the processes with the solute added. If you have given us a box that is a mixture, it will skip the pure solvent steps and focus on the mixture steps but you will need to give us some information about the solvent to perform the process. If your system is not a mixture, it will run through the pure solvent steps and ignore the mixture steps. This part of the script is work in progress and so should be avoided in testing for now.
    
- The next step is to initate the python script using Run_SCEE.py. The script will check if you are ready to start the process, if you say no it will give you the chance to make any changes to the simulation files or yaml. If you say yes the process will begin. The script will start by reading the Settings.yml, extracting the information like temperature and pressure, as well as using your mode and mixture settings. Based on the information you have given us the simulation will either start by generating the gro files and the box for your system, or will jump to the inital minimisation of your box. I think we will need to change how the input of the minimisation of the box works due to potentially being given a premade box.
- The script will generate the oniom input file for the system and perform 3 quick gaussian calculations using the NVT_Vacuum output. The output of these calculations are used in the calculation of induced dipole moment which is then used in the calculation of the dipole moment of the liquid. Our script will generate the oniom input file for the pure solvent but has not been designed for handling the difference between UA and AA systems and other molecules.
- If the system is a pure solvent, either because the box build or given is a pure solvent, it will start with the first loop in the script for performing the simulations and calculations. If the box given is a mixture, then it will   Currently the script will make a simulations folder where it will loop through the number of replicas you wish the system to perform, then through temperature and pressure. Each iteration in the loops is a MD simulation performed and the configurations extracted in the last 4 ns of the MD simulations for the purpose of using our SCEE dipole moment calculations. When iterations are completed the script will move onto performing the gaussian calculations on the extracted dipole moments which is then used for calculating the dielectric constant correction.
- If you go beyond this point you will find the mixture loop and analysis steps (still need to be added).
    
