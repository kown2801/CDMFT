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

# Inputs and outputs

## Self-consistency (CDMFT)

### If iteration!=0,

#### Inputs

`inputDirectory/{inputfilename}{iteration}.meas.json` is the main input file. 
It should contain certain fields for this part of the program : 
* Measurements
	* GreenI\_*component* (for example GreenI\_00) : Imaginary part of component *component* of the cluster Green's function. 
	* GreenR\_*component* (for example GreenR\_00) : Real part of component *component* of the cluster Green's function. 
	* Sign (the program divides all measured observables by this value at the start of the program. See more details in `ImpuritySolver/README.md`
	* All other observables that come out of the impuritysolver program. Those will be saved into dat files in order to have them already divided by the sign and access them easily.

* Parameters
	* Physical parameters
		* ep
		* mu
		* tpd
		* tppp (if not present, the program takes tppp=tpp)
		* tpp (Usually taken equal to 1)
		* beta
		* EGreen : Energy cutoff for the Green's function, maximum energy we consider in the simulation)
		* Optional parameters : S and n. When you define one you should define the other also. It is used to fix the occupation instead of the chemical potential mu. In order to fix n, the program changes mu by -S\\Delta n , \\Delta n being the difference between Measurements["N"] and n. Their use os however not recommended as convergence may be way slower.
	* File parameters. Beware here that the `inputDirectory/` of the ImpuritySolver in the `ouputDirectory/` of this program: 
		* HYB : name of the Hyb file in `outputDirectory/`
		* (LINKN and LINKA) or LINK : name of the Link files in `outputDirectory/`. If LINK is defined, it will use the file indicated to build the link matrix (correspondance between the cluster and the components defined in HYB). Else it will build the link matrix from the LINKN and LINKA files. See README.md in `ImpuritySolver/` for more info on this.


#### Outputs
* Files for the next iteration : 
	* `outputDirectory/{inputfilename}{iteration+1}.json`
	* `outputDirectory/Hyb{iteration+1}.json`
* Data files. Those files are : 
	* `dataDirectory/self{iteration}.json` Computed in `CDMFT`
	* `dataDirectory/green{iteration}.json` Computed in `CDMFT`
* Measurement files : 
	* `dataDirectory/ChiFull{iteration}.dat` result from the input file (Chi), divided by the sign and copied here
	* `dataDirectory/ChiFullSites{iteration}.dat` result from the input file (Chi_j), divided by the sign and copied here
	* `dataDirectory/pK{iteration}.dat` result from the input file (pK), divided by the sign and copied here
	* Chi0.dat Chi0Sites.dat D.dat DSites.dat k.dat kSites.dat N.dat NSites.dat sign.dat Sz.dat SzSites.dat. A line is added at the end of those files, all in `dataDirectory/` . This line has the structure : iteration_nb data1 data2 data3... This data is only the result taken from the input file, divided by the sign and copied here.


### Else if iteration==0,

#### Inputs

In this case the input file is `inputDirectory/{inputfilename}0.json`
The content of this file should be the same as the `Parameters` field in the case iteration!=0. There is then two possibilities : 
 * if `dataDirectory/` contains a self0.json file, it will load it and create the hybridization function with it
 * else it assumes the self-energy is zero for the normal part. For the anomalous part, it uses the `delta` parameter of the input file to initialize it to a non zero value (in order to allow superconductivity).
Its outputs are the same as for the iteration!=0 case except for the measurement files (as there are no results yet).


