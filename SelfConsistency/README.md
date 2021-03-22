This is the self consistency part of the program.
This program was written by Patrick Sémon (patrick.semon@gmail.com as of June 2019).
For license and use please advise him.
If you have trouble following the steps detailled here you can always contact me (nicolas.kowalski.2016@polytechnique.org) 

# Installation
In order to be able to install this program, you need to first follow the steps of the README.md file in the main directory
When this is all done and properly installed, come back in this directory and type the command : 

	make

That's it you are set up to use the self-consistency and full Green's function programs.

# Usage
Always use the command : 

	source ../scripts/export.sh

before using the program. This allows the computer to load the necessary libraries into your path as well as the right ComputeCanada modules.
When you installed the program, just before, you actually installed 3 different programs : 

## Self-consistency

The first program is called CDMFT and it allows to do the self-consistency equations to the CDMFT solution.
It uses a Green's function computed on the impurity (from the impurity solver) to get the self-energy and then the next hybridation function to feed into the impurity solver. 
To call the program, use :

	./CDMFT inputDirectory/ outputDirectory/ dataDirectory/ inputfilename iteration

(Don't forget the / after the directory names)

There are two ways of using this program : 
* if iteration equals 0, it uses the parameter file in inputDirectory (`inputfilename{iteration}.json` - for example `params0.json`) in order to create an initial hybridation file with an all-zero self-energy. It however initializes the anomal self-energy to delta/(1+i\omega_n^2) in order to permit supraconductivity (delta being defined in the params0.json file) and i\omega_n the Matsubara frequency. This case is only used if you want to start from scratch. The usual way of doing things is to reuse a Hyb file frome previous simulations.

* if iteration is non-zero, it uses the results/parameter file in `inputDirectory/` (`inputfilename{iteration}.meas.json` - for example `params0.meas.json`) AND the `Hyb{iteration}.json` file in `outputDirectory/` in order to compute the hybridation file for the next iteration. In this case to ensure a greater stability in the results, the new hybridation function is a weighted combination of the old hybridation function (weight w) and the hybridation function computed using the self-consistency relation (weight 1-w). This w is equal to weightR + i\*weightI defined in the results/parameter file.

In both cases, it creates a `inputfilename{iteration + 1}.json` (for example `params2.json` if iteration=1) file in `outputDirectory/` with the parameters of the next iteration. It also creates a `Hyb{iteration + 1}.json` file in this same directory. This allows to continue the self-consistency cycle.
This program also creates .dat and .json files in `dataDirectory/` in order to more easily access the results of the simulation. 
Notable files are:
* `N.dat` that provides the copper occupation per site per spin
* `D.dat` that provides the copper double occupation per site
* `self{iteration}.json` that contains the self-energy for the iteration
* `green{iteration}.json` that contains the copper Green's funcition for the iteration

Here I provide a minimal working example with : 

	./CDMFT IN/ OUT/ DATA/ params 0

and :

	./CDMFT IN/ OUT/ DATA/ params 70

## Including Oxygens

We never compute observables about the oxygens in the program. That's because we only look directly at copper sites. The effectof the oxygen sites is only included in the Hybridization function. In order to do look at what is happening on the oxygen, we use the GFULL program.
	
	./GFULL inputDirectory/ outputDirectory/ dataDirectory/ inputfilename iteration

(Don't forget the / after the directory names)

It computes the full Green's function on the 3 bands, the 2 spins and the 4 sites. It takes more memory and time than the self-consistency and this is why it is a separate program. Its inputs are the `inputDirectory/inputfilename{iteration}.meas.json` that comes out of the impurity-solver and the `dataDirectory/self{iteration}.json` that comes out of the CDMFT program. It then outputs many observables into dat files in the `dataDirectory/`.
Notable observables include :
* The oxygen occupation per site per spin in `pn.dat`
* The kinetic energy in `ekin.dat` 
* Some terms of the oxygen-oxygen green's function in pxgreen, pygreen and pxygreen files

Here I provide a minimal working example with : 

	./GFULL IN/ OUT/ DATA/ params 70

## Superfluid Stiffness

We provide a calculation of the superfluid stiffness using a self-energy on the copper cluster.

	./STIFFNESS inputDirectory/ outputDirectory/ dataDirectory/ inputfilename iteration

(Don't forget the / after the directory names)
The only output is a line in the `dataDirectory/stiffness.dat`.
On this line the terms are in order : 
* the iteration
* the real part of the stiffness
* the imaginary part of the stiffness (it should always be zero,I provide it just for verification purposes)

Here I provide a minimal working example with : 

	./STIFFNESS IN/ OUT/ DATA/ params 70
