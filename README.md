# CDMFT CT-QMC 3-band segment solver

This solver is used to solve the three band Hubbard Model (or Emery model) in the CDMFT approximation with a Continuous Time Quantum Monte-Carlo impurity solver.
This programm was written by Patrick Semon (patrick.semon[at]gmail[dot]com as of June 2019).
For license and use please ask him directly.
If you have trouble installing or using the program, you can contact the author of this REAMDE.md (nicolas.kowalski.2016[at]polytechnique[dot]org) 

## Installation

This file describes how to install the program on a Compute Canada cluster.
If you are not on such computers, you will need to install more things than explicited here.

### Requirements

This program relies on multiple libraries.
OpenBlas, pthread, nlohmann/json, openmpi

* pthread is used only by OpenBlas and not used directly in the present program.
* pthread and openmpi are already accessible from ComputeCanada clusters.
* Openblas is also accessible directly from ComputeCanada clusters but the version they have there doesn't seem to be able to compute the inverse of a matrix. Therefore, we need to recompile it by hand. I provide the instructions for installing OpenBlas
* The nlomann/json library is already installed in the project. It is only a single `.hpp` file. this file is called `nlohmann_json.hpp` in both the `SelfConsistency/` and `ImpuritySolver/` directories. It is a copy of the file `single_include/nlohmann/json.hpp` from commit n°176d8e261 of https://github.com/nlohmann/json.git. I just changed the dump function for arrays in order to make them oneline. It makes the `.json` files more readable.

You need to install one library. If you are lazy you can just execute the `./install-libraries.sh` script to install it. This will execute all command lines indicated here until `cp libopenblas.a ~/local/lib` included. It may take a while (installing openblas is a bit long).

Let's go.
First you need to create a "local" folder in your home directory. Then inside of this folder, create a `lib/` directory.
	
	cd
	mkdir local
	cd local
	mkdir lib


### Installing OpenBlas : 

go into the 'local' folder and clone the Openblas github repository

	git clone https://github.com/xianyi/OpenBLAS.git

then go into the directory and get the version 0.3.9 of Openblas (that is the last version I have seen that compiles and runs correctly on Cedar)

	cd OpenBLAS
	git checkout v0.3.9
	
and make the program

	make -j4

(The j options allows you to compile on multiple processors simultaneously, so that is goes faster)	
This step may be a bit long, so go grab a cup of coffee.

Then copy the freshly compiled library to your lib folder 

	cp libopenblas.a ~/local/lib 

The library installation part is finished. Now you still need to make the program. In order to do that, go to the Readme files inside the two folders (ImpuritySolver, SelfConsistency). Don't forget to do this much as you wouldn't be able to run the program without it.
This is the minimal structure of the ${HOME}/local folder you need.

	~/local
	├── lib
		├── libopenblas.a 

## Usage

CMDFT consists in iterating an system of equations. Here, you use `ImpuritySolver/IS` and `SelfConsistency/CDMFT` alternatively until convergence is achieved. The script in `scripts/launch.py` can be used to do the cycle automatically. The other files in the `scripts/` directory can be used to launch a new simulation or execute different actions on simulations. 

You can either use the programs using what is written in the [SelfConsistency/README.md](SelfConsistency/README.md) and [ImpuritySolver/README.md](ImpuritySolver/README.md) files or use the python library that I wrote in order to launch simulations. 
I describe here how to setup the python scripts on your local computer to start using the code now.
All the scripts I write about here are located in the `localScripts/` folder at the root of the repository. In order to use those, copy the content of this directory on your personal - or work - computer. Please make sure to copy all the content of the directory because the structure of the code is important for it to work.

### Setting up your ssh keys

Then you will need to create a ssh key in order to connect to the ComputeCanada computers. It is fairly easy don't worry.
This tutorial is a simplified version of https://docs.computecanada.ca/wiki/Using_SSH_keys_in_Linux.
First on your local compute use the command : 

	ssh-keygen -b 4096 -t rsa

These lines should appear : 

	Generating public/private rsa key pair.
	Enter file in which to save the key (/home/username/.ssh/id_rsa):

Press enter if it is ok. Then you can use a password to protect you connection if you wish (not mandatory):

	Enter passphrase (empty for no passphrase):
	Enter same passphrase again:

Finally : 

	Your identification has been saved in /home/username/.ssh/id_rsa.
	Your public key has been saved in /home/username/.ssh/id_rsa.pub.
	The key fingerprint is:
	ef:87:b5:b1:4d:7e:69:95:3f:62:f5:0d:c0:7b:f1:5e username@hostname
	The key's randomart image is:
	+--[ RSA 2048]----+
	|                 |
	|                 |
	|           .     |
	|            o .  |
	|        S    o o.|
	|         .  + +oE|
	|          .o O.oB|
	|         .. poo+*|
	|          ... o..|
	+-----------------+

There, you just created a pair of rsa keys. 
Then you need to tell your local computer to use your key when it wants to connect to a cluster via ssh. In order to do this, use those two commands : 

	eval `ssh-agent`
	ssh-add

The final step is to copy the public key on the compute canada computer. To do so create a '.ssh/' directory in your home and a 'authorized_keys' file inside this directory.

	mkdir ~/.ssh
	touch ~/.ssh/authorized_keys

Then you need to copy the content of the 'id_rsa.pub' local file that you just created using the ssh-keygen command inside this 'authorized_keys' file (it is one line). You should be able to connect to your supercomputer space without a password now. If not so, please check https://docs.computecanada.ca/wiki/Using_SSH_keys_in_Linux for help. This step is mandatory in order to be able to connect from the jupyter notebooks.

### Setup before using the python library

The `localScripts/` directory should be copied on your computer. It includes multiple jupyter notebook files that will help you launch simulations and get/analyse data at the end of simulations. Keep the `AllData/` and `transfered` directories as results of simulations will be downloaded into them.


You have multiple files in order to do multiple things.
Before doing anything, please read those lines and follows those simple steps

* You NEED to change the first part of the distant_consts.py file in order to reflect your credentials. 

  * The variable distant_main_dir should be the main directory of your simulation space (in this directory one should find the structure of the repository - or zip file)

* You will also need to change the account that submits the jobs. In order to do that, you will need to change it at multiple locations in the "scripts" and "scripts/BACKUP_PATH" folder (approximately all the .sh files)

* Don't forget to change the number of nodes and of processor per nodes in scripts/BACKUP_PATH/parameter_run.sh to suit the super-computer you use and your needs.

### Scripts description
We are all set.
Here is a quick description of all the files that you can use. More detailed descriptions are given inside the individual jupyter notebooks.
First the two programs that you will use the most : 

* Fast-distant allows one to monitor all the current simulations from one place. It allows to almost do any operation to a simulation.

* New simulation allows one to start a new simulation from a Hyb seed

Then some programs that may help you but that you will use less often : 

* Convergence distant allows you to monitor the convergence folder by folder
* Convergence local allows you to check the behavior of your old simulations
* Resume Simulation allows you to resume your old simulation

Finally, I give my program that allows one to study easily the results that come out of the program.
This file is called Order-Parameter and its decription is also located at the beginning of the file. 
It is very convenient to plot graphs and manipulate simulation data.



