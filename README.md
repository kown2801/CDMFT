# Foobar

Foobar is a Python library for dealing with word pluralization.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install foobar
```

## Usage

```python
import foobar

foobar.pluralize('word') # returns 'words'
foobar.pluralize('goose') # returns 'geese'
foobar.singularize('phenomena') # returns 'phenomenon'
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)


# CDMFT CT-QMC 3-band segment solver

This solver is used to solve the three band Hubbard Model (or Emery model) in the CDMFT approximation with a Continuous Time Quantum Monte-Carlo impurity solver.
This programm was written by Patrick Semon (patrick.semon@gmail.com as of June 2019).
For license and use please advise him.
If you have trouble installing or using the program, you can contact the author of this REAMDE.md (nicolas.kowalski.2016@polytechnique.org) 

## Installation

This file describes how to install the program on a Compute Canada cluster.
If you are not on such computers, you will need to install more things than explicited here.
### Requirements

This program relies on multiple libraries.
Boost, OpenBlas, pthread, json_spirit, openmpi

pthread is used only by OpenBlas and not used directly in the present program.



Let's go.
First you need to create a "local" folder in your home directory. Then inside of this folder, create two directories "lib" and "include".
If at some point you don't understand the steps or the organisation of the folders, check the end of this installation tutorial where you have the finale structure of the local folder.
You need to install two libraries.
I. Installing OpenBlas : 

go into the "local" folder and clone the Openblas github repository

	$ git clone https://github.com/xianyi/OpenBLAS.git

then go into the directory 

	$ cd OpenBLAS
	
and make the program

	$ make

Then you have to locate the file called something like 

	libopenblas_barcelonap-r0.3.10.dev.a 

and move it to the "lib" directory. Then rename it **libopenblas.a**.
So you should have a file at ${HOME}/local/lib called libopenblas.a (where ${HOME} is your home directory).

II. Installing Json_spirit

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

You then need to modify the json_spirit_value.h file. Replace the lines : 
	
	#if defined(_MSC_VER)
	  // comment out the value types you don't need to reduce build times and intermediate file sizes
	  #define JSON_SPIRIT_WVALUE_ENABLED
	  #define JSON_SPIRIT_MVALUE_ENABLED
	  #define JSON_SPIRIT_WMVALUE_ENABLED
	#endif

by 

	#define JSON_SPIRIT_MVALUE_ENABLED

The installation part is finished. Now you still need to make the program. In order to do that, go to the Readme files inside the two folders (ImpuritySolver, Autocoherence). Don't forget to do this much as you wouldn't be able to run the program without it.
This is the minimal structure of your ${HOME}/local folder you need.

${HOME}/local ==>	lib 	==> libopenblas.a
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
%========= USAGE ==========%

You can either use the programs using what is written in the Autocoherence/README.md and ImpuritySolver/README.md files or use the python library that I wrote in order to launch simulations. 
We will describe here the second part. 
All the scripts I write about here are located in the localScripts folder at the root of the repository. In order to use those, copy the content of this directory on your personal - or work - computer. Please make sure to copy all the content of the directory because the structure of the code is important for it to work.
Then you will need to create a ssh key in order to connect to the ComputeCanada computers. It is fairly easy don't worry.
First on your local compute use the command : 
$ ssh-keygen -b 4096 -t rsa
These lines should appear : 

	Generating public/private rsa key pair.
	Enter file in which to save the key (/home/username/.ssh/id_rsa):

Press enter if it is ok. Then :

	Enter passphrase (empty for no passphrase):

You can use a password to protect you connection if you wish : 

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
	|         .. +oo+*|
	|          ... o..|
	+-----------------+

There you created a pair of rsa keys. 
The final step is to copy the public key on the compute canada computer. To do so create a .ssh directory in your home and a authorized_keys file in thid directory

	$ mkdir ~/.ssh
	$ touch ~/.ssh/authorized_keys

Then you need to copy the content of the id_rsa.pub local file that you just created using the ssh-keygen command inside this authorized_keys file (it is one line). You should be able to connect to your supercomputer space without a password now. If not so, please check https://docs.computecanada.ca/wiki/Using_SSH_keys_in_Linux for help. This step is mandatory in order to be able to connect from the jupyter notebooks.

You have multiple files in  order to do multiple things.
Before doing anything, please read those lines and follows those simple steps

* You NEED to change the first part of the distant_consts.py file in order to reflect your credentials. 

* The variable distant_main_dir should be the main directory of your simulation space (in this directory one should find the structure of the repository - or zip file)

* You will also need to change the account that submits the jobs. In order to do that, you will need to change it at multiple locations in the "scripts" and "scripts/BACKUP_PATH" folder (approximately all the .sh files)

* Don't forget to change the number of nodes and of processor per nodes in scripts/BACKUP_PATH/parameter_run.sh to suit the super-computer you use and your needs.

We are all set.
Here is a quick description of all the files that you can use. More detailed descriptions are given inside the individual jupyter notebooks.
First the two programs that you will use the most : 

* Fast distant allows one to monitor all the current simulations from one place. It allows to almost do any operation to a simulation.

* New simulation allows one to start a new simulation from a Hyb seed

Then some programs that may help you but that you will use less often : 

* Convergence distant allows you to monitor the convergence folder by folder
* Convergence local allows you to check the behavior of your old simulations
* Resume Simulation allows you to resume your old simulation

Finally, I give my program that allows one to study easily the results that come out of the program.
This file is called Order-Parameter and its decription is also located at the beginning of the file. It is very convenient to plot 2D and 3D graphs.



