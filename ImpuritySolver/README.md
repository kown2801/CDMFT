This is the impurity solver programm
This programm was written by Patrick Sémon (patrick.semon@gmail.com as of June 2019).
For license and use please advise him.
If you have trouble following the steps you can always contact me (nicolas.kowalski.2016@polytechnique.org) 

%========= INSTALLATION ==========%
In order to be able to install this program, you need to first follow the steps of teh README.md file in the main directory
When this is all done and properly installed, come back in this directory and type the two commands : 
$ . export.sh
$ make
There will maybe be a warning message about destructors. Just ignore it.
That's it you are set up to use this impurity solver.

%========= USAGE ==========%
Always use the command : 
$ . export.sh
before using the program in order to load the libraries into your path as well as the right computecanada modules.

	-------- One processor ------------

	To call the program, use 

$ ./IS input_folder/ output_folder/ FILE_NAME
	(Don't forget the / after the directory names)

	We provide a minimal working example in the folder. Calling  

$ ./IS IN/ OUT/ params1 will launch the program here.

	in the input_folder folder there should be at least 4 files : 
	- FILE_NAME.json contains the parameters of the simulation (physical and monte-carlo parameters)
		Inside this file are defined the HYB, LINKA and LINKN variable that represent the path to files from the input_folder/ directory.
		Thos files are json files. We provide an example of such files at IN/Hyb1.json, IN/LinkA.json, IN/LinkA.json
	- HYB is the "hybridization seed" of the program. It is the hybridisation function of the bath needed for the program.
	- LINKA and LINKN shouldn't generally be changed. You can change them to change the structure of the Hamiltonian and the Symmetries. (here we allow for superconductivity but not anti-ferrmagnetism).
	The simulation will output in the output_folder folder a file FILE_NAME.meas.json
	It will also output a config1.json file that represents the state of the operators in the segment picture for the processor at the end of the simulation. If such a file exists in the OUT folder at the beginning of a simulation, the programm will use this file as the starting point of the Markov Chain. This is useful to reduce the overall need for thermalization.

	---------- Multiple processors -----------

	In order to reduce the time needed for a simulation (because a simulation requires a lot of Monte-Carlo steps), you can parallelize the simulation. Programmers would call this kind of parallelization "Brain Dead Parallelization". What is actually done is running independant simulations simultaneously and then taking the average over all those simulations. In a usual simulation with this program at temperature 60 and for a sign of 0.1, 96 processors and 30 minutes of simulation are needed to reach convergence. 
	In order to launch this kind of simulation, you can use srun on a slurm based cluster. This is the case for example on ComputeCanada cluster. We provide a batch file run.sh that allows one to run the program on 96 processors across 3 nodes (this works on ComputeCanada:cedar, but you should choose the number of processors and nodes according to your machine). Don't forget to change the other SBATCH parameters to fit your situation. On cedar for example the command to launch this job would be : 
	$ sbatch run.sh
	Running a parallel job will not create one single config1.json file but configm.json files where m is the processor number (so 96 config file for 96 processors). Just as in the case for 1 processors, if you provide the same number of config files as processors you use, the program will load the segment picture in order to reduce the thermalization time needed.
	In order to make the simulations independent, thermalization is done on all processors independently.




