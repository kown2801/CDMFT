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
* `inputDirectory/inputfilename.json` should contain the parameters of the simulation (physical and Monte-Carlo parameters)
	I provide an example in `IN/params1.json`. 
* `inputDirectory/hybfilename.json` is the **hybridization seed** of the program. It is the hybridisation function of the bath from which the program compute the Impurity Green's function. Inside the `inputDirectory/inputfilename.json` file the `HYB` field should read the name of the hybridization from the `IN/` folder. In the example it is `Hyb1.json`
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


# Inputs and Outpus

## Inputs

`inputDirectory/inputFilename.json` is the main parameter file. It should contain at least:
* Physical parameter : 
	* U : double occupation energy)
	* mu : chemical potential
	* beta : inverse temperature)
	* EGreen : Energy cutoff for the Green's function, maximum energy we consider in the simulation)
	* EObs : Energy cutoff for the recorded observables (used for frequency dependant quantities, like Chi for example)
	* EHyb : Energy cutoff for the Hybridization function
	* EOrder : Maximum number of operators recorded on the time lines. This does not limit the number of operators inside the progam but only for the output.
* File parameters : 
	* HYB : name of the Hyb file in `inputDirectory/`
	* (LINKN and LINKA) or LINK : name of the Link files in `inputDirectory/`. If LINK is defined, it will use the file indicated to build the link matrix (correspondance between the cluster and the components defined in HYB). Else it will build the link matrix from the LINKN and LINKA files. See above for more info on this
* Monte-Carlo parameters:
	* THERMALIZATION_TIME in minutes
	* MEASUREMENT_TIME in minutes
	* SEED of the random number generator
	* CLEAN_EVERY_SWEEP : number of Monte-Carlo sweeps between every cleaning of the bath matrix. Cleaning the bath matrix means recomputing it from scratch. THis is done to avoid numercial instabilities in the Shermann-Morisson formula.
	* SAMPLE_EVERY_SWEEP : number of Monte-Carlo sweeps between every measurement (or sample)
	* STORE_EVERY_SAMPLE : number of measurements between every save in the binning procedure (used to avoid storing the results too often)
	* PROBFLIP : probability of a flip sweep. This is used to allow the program to go into the whole integration space
	
`inputDirectory/{inputDirectory/inputFilename.json["HYB"]}` is the hybridation file. The structure should be like the example given in the folder.
All the components indicated in (LINKN and LINKA) or LINK should exist (except `empty`)

(`inputDirectory/{inputDirectory/inputFilename.json["LINKN"]}` and `inputDirectory/{inputDirectory/inputFilename.json["LINKA"]}`) or `inputDirectory/{inputDirectory/inputFilename.json["LINK"]}` are the link file. They describe how the components in the HYB file describe the physical cluster. They are 2D arrays of string. You have multiple examples in the`IN/` folder.

It contains other physcial and simulation parameters. Those are used in the self-consistency part of the iteration cycle (see `SelfConsistency/`)

## Outputs

`outputDirectory/inputFilename.meas.json` is the output file. It contains the measured values of the observables and their errors. It also contains little information about the simulation itself. You can see an example in `OUT/`. Because of the sign problem, we sample separately sign\*observable and sign. In order to get the true value of the observables printed in the file, you should divided them by the value of the Sign field. This is true for the Errors and Measurements fields.

The fields inside the file are at this stage : 
* Errors contains the errors on the observables. This is the maximum value from the binning procedure that is used to decorrelate the measurements.
* Measurements contains the measured observables. They are at this stage : 
	* Chi as a function of matsubara frequencies : Site spin susceptibilities. This is the Fourrier transform with respect to time of <S^z_i(\tau) S^z_i>
	* Chi0 is the first component of Chi : Chi(\iq_n=0)
	* When measurements end with \_j, it means they are the observable on site j.
	* Chiij is <S^z_i(\tau) S^z_j>. It is flattened here. It is ordered first by sites i,j and then by matsubara frequencies. When reading it, you should then reshape it. For a 4 sites cluster, the actual shape should be (4,4,beta*EObs/(2*pi)+2)). As an example, the first beta*EObs/(2*pi)+2 numbers should be exactly equal to Chi_0.
	* D : double occupation `<n_\uparrow n_\downarrow>`
	* GreenI\_*component* (for example GreenI\_00) : Imaginary part of component *component* of the cluster Green's function. This observable is the main result of the impurity Solver. It is what allows the iteration cycle to continue.
	* GreenR\_*component* (for example GreenR\_00) : Real part of component *component* of the cluster Green's function. This observable is the main result of the impurity Solver. It is what allows the iteration cycle to continue.
	* N : occupation on the cluster per site per spin (betwwen 0 and 1)
	* Sign : sign of the simulation. The true value of all the observables is \frac{observable}{Sign}
	* Sz : spin
	* k : mean expansion order, mean number of operators on the time lines (this is linked to the kinetic energy of the cluster)
	* pK : observed probability of each order appearing on the time lines. pK[0] for example is the observed probability of having no operators on the time lines.
* Job Specifications. These are pretty self-explanatory
	* Measurement Sweeps per Processor
	* Number of Processors
	* Thermalization Sweeps per Processor
* Parameters. It simply contains a copy of `inputDirectory/inputFilename.json`. 

The program also outputs time line configuration files `outputDirectory/config_{processor_id}.json`. At the end of the simulation, the program saves the current time line configuration in order to decrease the necessary thermalization time for the next iteration. Those files are not meant to be used outside of the program and can be deleted once the solution is converged.