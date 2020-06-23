This is the impurity solver program
This program was written by Patrick Sémon (patrick.semon@gmail.com as of June 2019).
For license and use please advise him.
If you have trouble following the steps detailled here you can always contact me (nicolas.kowalski.2016@polytechnique.org) 

%========= INSTALLATION ==========%






%========= USAGE ==========%
Always use the command : 
$ . export.sh
before using the program in order to load the libraries into your path as well as the right computecanada modules.
This program is actually composed of two different programs. 

	------ Self-consistency -----
	The first program called CDMFT allows to do the self-consistency equations and using a Green's function computed on the impurity to get the self-energy and then the next hybridation function to feed into the impurity solver. 
	To call the program, use :

$ ./CDMFT inputDirectory/ outputDirectory/ dataDirectory/ inputfilename iteration
	(Don't forget the / after the directory names)

	there are two cases of use of this program : 
	- if iteration equals 0, it uses the parameter file in inputDirectory (inputfilenameiteration.json - for example params0.json) in order to create an initial hybridation file with an all-zero self-energy. It however initializes the anomal self-energy to delta/(1+i\omega_n^2) in order to permit supraconductivity (delta being defined in the params0.meas.json file) and i\omega_n the Matsubara frequency. This case is only used if you want to start from scratch. The usual way of doing things is to reuse a Hyb file frome previous simulations.

	- if iteration is non-zero, it uses the results/parameter file in inputDirectory (inputfilenameiteration.meas.json) AND the Hybiteration.json file in outputDirectory in order to compute the hybridation file for the next iteration. In this case to ensure a greater stability in the results, the new hybridation function is a weighted combination of the old hybridation function (weight w) and the hybridation function computed using the self-consistency relation (weight 1-w). This w is equal to weightR + i*weightI defined in the results/parameter file.

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
