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

there are two cases of use of this program : 
* if iteration equals 0, it uses the parameter file in inputDirectory ('inputfilename{iteration}.json' - for example 'params0.json') in order to create an initial hybridation file with an all-zero self-energy. It however initializes the anomal self-energy to delta/(1+i\omega_n^2) in order to permit supraconductivity (delta being defined in the params0.json file) and i\omega_n the Matsubara frequency. This case is only used if you want to start from scratch. The usual way of doing things is to reuse a Hyb file frome previous simulations.

* if iteration is non-zero, it uses the results/parameter file in inputDirectory (inputfilename{iteration}.meas.json - for example 'params0.json') AND the Hybiteration.json file in outputDirectory in order to compute the hybridation file for the next iteration. In this case to ensure a greater stability in the results, the new hybridation function is a weighted combination of the old hybridation function (weight w) and the hybridation function computed using the self-consistency relation (weight 1-w). This w is equal to weightR + i\*weightI defined in the results/parameter file.

	In both cases, it creates a inputfilename(iteration + 1).json (for example params1.json if iteration=0) file in outputDirectory with the parameters of the next iteration. It also creates a Hyb(iteration + 1).json file in this same directory. This allows to continue the self-consistency cycle.
	This program also creates some .dat files in dataDirectory (greeniteration.dat, selfiteration.dat and hybiteration.dat or others that you may inspect at your will) (TODO Maybe make a list in future versions of this readme)

	Here we provide a minimal working example with : 
$ ./CDMFT IN/ OUT/ DATA/ params 0

	and :
$ ./CDMFT IN/ OUT/ DATA/ params 1

	-------- Including oxygens -------
	We never compute observables about the oxygens in the program. In order to do so, we use the GFULL program.
	
$ ./GFULL inputFolder/ dataFolder/ filename iteration
	(Don't forget the / after the directory names)

	It computes the full Green's function on the 3 bands, the 2 spins and the 4 sites. It takes more memory and time than the self-consistency and this is why it is a separate program. Its input are the inputFolder/filenameiteration.meas.json that comes out of the impurity-solver and the dataFolder/selfiteration.dat that comes out of the CDMFT program. It then outputs many observables into dat files in the dataFolder/ (for example the oxygen occupation, the kinetic energy or the oxygen-oxygen and the copper-oxygen Green's function) 

	Here we provide a minimal working example with : 
$ ./GFULL IN/ DATA/ params 1
