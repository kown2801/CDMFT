#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.8.0/bin/python3.8
import os
import json
import shutil
import sys
import subprocess

def create_dir_safely(dir):
	try:  
	    os.mkdir(dir)
	except OSError:  
		if not os.path.exists(dir):
			print("Creation of the directory " + dir + " failed\n")
			print("Sorry, we need this directory to proceed, stopping... \n")
			exit()

def find_available_data_dir_name(basename):
	current_name = basename + ""
	i = 1
	while os.path.isdir(current_name):
		current_name = basename + "_" + str(i)
		i+=1
	return current_name

def generate_simulation(all_args):
	this_directory_from_file_dir = "../../scripts"
	os.chdir(os.path.dirname(os.path.realpath(__file__))) #The working directory is now the root of the repository
	#Create the name of the folder to host this new simulation
	data_dir_name = "ep" + all_args["ep"] + "_beta" + all_args["beta"] + "_mu" + all_args["mu"] + "_U" + all_args["U"]
	if "tpd" in all_args:
		data_dir_name += "_tpd" + all_args["tpd"]
	if "tppp" in all_args:
		data_dir_name += "_tppp" + all_args["tppp"]

	backup_path = "./BACKUP_START"
	path_to_main_dir = "../"
	datas_dir =  os.path.join(path_to_main_dir,"ComputedData")
	create_dir_safely(datas_dir)
	files_path = os.path.join(datas_dir,data_dir_name)
	#We want to create a directory that does not exist for now (no overlap)
	files_path = find_available_data_dir_name(files_path)
	create_dir_safely(files_path)

	#We create the directories needed for simulations
	create_dir_safely(os.path.join(files_path,"DATA"))
	create_dir_safely(os.path.join(files_path,"IN"))
	create_dir_safely(os.path.join(files_path,"OUT"))

	#We put the Hyb file in
	shutil.copy(os.path.join(backup_path,"Hyb1.json"),os.path.join(files_path,"IN/Hyb1.json"))
	
	#From here we create the run.sh file to submit to slurm
	sh_origin = open(os.path.join(backup_path,"parameter_run.sh"),"r")
	lines = sh_origin.readlines()
	sh_origin.close()
	sh_destination = open(os.path.join(files_path,"run.sh"),"w")
	sh_destination.write('#!/bin/bash\n')
	sh_destination.write("#SBATCH --time=" + all_args["computing_time"] + "-00:00:00\n")
	del all_args["computing_time"]
	sh_destination.write('#SBATCH --job-name="' + all_args["ep"] + all_args["U"] + all_args["mu"] + '"\n')
	for l in lines:
		sh_destination.write(l)
	sh_destination.write("python " + os.path.join(this_directory_from_file_dir,"launch.py") + " " + files_path + " " + str(all_args["iterations"]) + " 1")
	del all_args["iterations"]
	sh_destination.close()
	#End creation run.sh

	#From here we create the json parameter file
	f = open(os.path.join(backup_path,"params0.meas.json"),"r")
	params_json = json.loads(f.read())
	f.close()
	params_json = params_json["Parameters"]
	for i in all_args:
		if i not in params_json:
			raise Exception(i + " can't be included in the params file if it is not declared in the one in scripts/BACKUP_START")
		params_json[i] = type(params_json[i])(all_args[i])
	params_json["HYB"] = "Hyb1.json"
	f = open(os.path.join(files_path,"IN/params1.json"),"w")
	if "LINK" in params_json:
		shutil.copy(os.path.join(backup_path,"Link.json"),os.path.join(files_path,"IN/Link.json"))
	else:
		shutil.copy(os.path.join(backup_path,"LinkA.json"),os.path.join(files_path,"IN/LinkA.json"))
		shutil.copy(os.path.join(backup_path,"LinkN.json"),os.path.join(files_path,"IN/LinkN.json"))
	f.write(json.dumps(params_json, indent=4,sort_keys=True))
	f.close()
	#End creation param1.json

	#Begin start simulation
	os.chdir(files_path)
	subprocess.run(["sbatch","run.sh"])
	#End start simulation


if __name__ == "__main__":
	all_args = dict()
	i = 1
	while len(sys.argv) > i+1:
		all_args[sys.argv[i]] = sys.argv[i+1]
		i+=2
	print(all_args)
	if "ep" in all_args and "beta" in all_args and "mu" in all_args and "U" in all_args and "computing_time" in all_args and "iterations" in all_args:
		generate_simulation(all_args)
	else:
		print("Usage : new_simulation.py parameter_name1 parameter_value1 parameter_name2 parameter_value2 ...\n At least ep, beta, mu, U, computing_time and iterations are required")
