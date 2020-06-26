This is the impurity solver programm
This programm was written by Patrick S�mon (patrick.semon@gmail.com as of June 2019).
For license and use please advise him.
If you have trouble following the steps you can always contact me (nicolas.kowalski.2016@polytechnique.org) 

%========= INSTALLATION ==========%
This file describes how to install the program in a compute canada computer. If you are not on such computers, you will need to install more things than explicited here.
Let's go.
First you need to create a "local" folder in your home directory. Then inside of this folder, create two directories "lib" and "include".
If at some point you don't understand the steps or the organisation of the folders, check the end of this installation tutorial where you have the finale structure of the local folder.
You need to install two libraries.
Installing OpenBlas : 

	go into the "local" folder and clone the Openblas github repository

$ git clone https://github.com/xianyi/OpenBLAS.git

	then go into the directory 

$ cd OpenBLAS
	
	and make the program

$ make

Then you have to locate the file called something like libopenblas_barcelonap-r0.3.10.dev.a and move it to the "lib" directory. Then rename it to libopenblas.a.
So you should have a file at ${HOME}/local/lib called libopenblas.a (where ${HOME} is your home directory).
Then you have to copy the "cblas.h", "f77blas.h", "lapacke_mangling.h", "openblas_config.h" files into a "openblas" folder in the include folder.
Those files are located at the root of the OpenBlas folder.

Installing Json_spirit

	Again in the local folder, type : 

$ git clone https://github.com/png85/json_spirit.git

	Then you have to prepare for compilation : 

$ cd json_spirit
$ mkdir build
$ cd build/
$ module load boost
$ cmake .. -DJSON_SPIRIT_MVALUE_ENABLED:BOOL=ON

	You can then make the library

$ make -j4

	Copy the libjson_spirit.so file from the json_spirit/build/json_spirit folder into your lib folder
	Copy all the .h files from the json_spirit/json_spirit folder into a json_spirit folder in your include folder (there should be 10 files).

	You then need to modify the json_spirit_utils.h file. Replace the lines : 
	
#if defined(_MSC_VER)
  // comment out the value types you don't need to reduce build times and intermediate file sizes
  #define JSON_SPIRIT_WVALUE_ENABLED
  #define JSON_SPIRIT_MVALUE_ENABLED
  #define JSON_SPIRIT_WMVALUE_ENABLED
#endif

by #define JSON_SPIRIT_MVALUE_ENABLED

The installation part is finished. Now you still need to make the program. In order to do that, go to the Readme files inside the two folders (ImpuritySolver, Autocoherence). Don't forget to do this much as you wouldn't be able to run the program without it.

${HOME} ==>	local ==>	lib 	==> libopenblas.a
								==>	libjson_spirit.so

						include ==> json_spirit ==>json_spirit.h
												==>json_spirit_error_position.h
												==>json_spirit_reader.h
												==>json_spirit_reader_template.h
												==>json_spirit_stream_reader.h
												==>json_spirit_utils.h
												==>json_spirit_value.h
												==>json_spirit_writer.h
												==>json_spirit_writer_options.h
												==>json_spirit_writer_template.h

								==> openblas 	==>	cblas.h
												==> f77blas.h 
												==> lapacke_mangling.h
												==> openblas_config.h

%========= USAGE ==========%

You can either use the programs using what is written in the Autocoherence/README.md and ImpuritySolver/README.md files or use the python library that I wrote in order to launch simulations. 
We will describe here the second part. 
All the scripts I write about here are located in the localScripts folder at the root of the repository. In order to use those, copy the content of this directory on your personal - or work - computer. Please make sure to copy all the content of the directory because the structure of the code is important for it to work.

You have multiple files in  order to do multiple things.
Before doing anything, please read those lines and follows those simple steps

	You NEED to change the first part of the distant_consts.py file in order to reflect your credentials. 
	The variable distant_main_dir should be the main directory of your simulation space (in this directory one should find the structure of the repository - or zip file)

	You will also need to change the account that submits the jobs. In order to do that, you will need to change it at multiple locations in the "scripts" and "scripts/BACKUP_PATH" folder (approximately all the .sh files)
	Don't forget to change the number of nodes and of processor per nodes in scripts/BACKUP_PATH/parameter_run.sh to suit the super-computer you use and your needs.
We are all set.
Here is a quick description of all the files that you can use. More detailed descriptions are given inside the individual jupyter notebooks.
First the two programs that you will use the most : 

	- Fast distant allows one to monitor all the current simulations from one place. It allows to almost do any operation to a simulation.

	- New simulation allows one to start a new simulation from a Hyb seed

Then some programs that may help you but that you will use less often : 

	- Convergence distant allows you to monitor the convergence folder by folder
	- Convergence local allows you to check the behavior of your old simulations
	- Resume Simulation allows you to resume your old simulation

Finally, I give my program that allows one to study easily the results that come out of the program.
This file is called Order-Parameter and its decription is also located at the beginning of the file. It is very convenient to plot 2D and 3D graphs.



