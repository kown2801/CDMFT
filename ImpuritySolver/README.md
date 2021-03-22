This is the impurity solver part of the programm
This programm was written by Patrick Sémon (patrick.semon@gmail.com as of June 2019).
For license and use please advise him.
If you have trouble following the steps you can always contact me (nicolas.kowalski.2016@polytechnique.org) 

# Installation
In order to be able to install this program, you need to first follow the steps of the README.md file in the main directory
When this is all done and properly installed, come back in this directory and type the command : 

	make

There will maybe be a warning message about destructors. Just ignore it.
That's it you are set up to use this impurity solver.

# Usage
Always use the command : 

	source ../scripts/export.sh

before using the program. This allows the computer to load the necessary libraries into your path as well as the right ComputeCanada modules.

## One processor

To call the program, use 

	./IS inputDirectory/ outputDirectory/ inputfilename

(Don't forget the / after the directory names)

I provide a minimal working example in the folder. Calling  

	./IS IN/ OUT/ params70

 will launch the program here.

In general, for the program to work, there should be :  
* `inputfilename.json` contains the parameters of the simulation (physical and Monte-Carlo parameters)
	I provide an example in `IN/params1.json`. 
* `hybfilename.json` is the **hybridization seed** of the program. It is the hybridisation function of the bath from which the program compute the Impurity Green's function. Inside the `inputfilename.json` file the `HYB` field should read the name of the hybridization from the `IN/` folder. In the example it is `Hyb1.json`
* Link files. Those represent the structure of the Hamiltonian. There is two possibilities : 
	* if "LINK" is defined in `inputfilename.json`, then the file at field "LINK" should contain a 2D array of size `(2*clusterSize,2*clusterSize)` to include the spin degrees of freedom. See the example in `IN/Link_antiferro.json` for the antiferromagntic phase.
	* if not, "LINKA" and "LINKN" should be defined in `inputfilename.json`. The file at field "LINKN" should contain a 2D array of size `(clusterSize,clusterSize)` representing the same spin Green's function. The file at field "LINKA" should contain a 2D array of size `(clusterSize,clusterSize)` representing Green's function for opposite spins. Examples can be found at `IN/LinkN.json` and `IN/LinkA.json`. Inside the program, the link matrix is created by tiling the LINKN and LINKA matrices : 
```
	LINK =  LINKN	LINKA
		LINKA	LINKN
```
The simulation will output in the outputDirectory/ a file inputfilename.meas.json
It will also output a config1.json file that represents the state of the operators in the segment picture for the processor at the end of the simulation. If such a file exists in the OUT folder at the beginning of a simulation, the programm will use this file as the starting point of the Markov Chain. This is useful to reduce the overall need for thermalization.
## Multiple processors 

In order to reduce the time needed for a simulation (because a simulation requires a lot of Monte-Carlo steps), you can parallelize the simulation. Programmers would call this kind of parallelization "Brain Dead Parallelization". What is actually done is running independant simulations simultaneously and then taking the average over all those simulations. In a usual simulation with this program at temperature 60 and for a sign of 0.1, 96 processors and 30 minutes of simulation are needed to reach convergence. 
In order to launch this kind of simulation, you can use srun on a slurm based cluster. This is the case for example on ComputeCanada cluster. I provide a batch file run.sh that allows one to run the program on 4 processors as an example and test. The `nodes` parameter is commented for this test but feel free to use it in you future calculations of course. Don't forget to change the other SBATCH parameters to fit your situation. On cedar for example the command to launch this job would be : 

	sbatch run.sh

Running a parallel job will not create one single config1.json file but configm.json files where m is the processor number (so 96 config file for 96 processors). Just as in the case for 1 processors, if you provide the same number of config files as processors you use, the program will load the segment picture in order to reduce the thermalization time needed. In order to make the simulations independent, thermalization is done on all processors independently.

On thing to remember here : never use a config file created for a different number of processors. This may result in crashes of the application.




