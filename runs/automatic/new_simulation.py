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
			launch.write_to_log_file("Creation of the directory " + dir + " failed\n")
			launch.write_to_log_file("Sorry, we need this directory to proceed, stopping... \n")
			exit()
	to_log("Successfully created the directory " + dir + " \n")

def to_log(text):
	f = open("logfile","a")
	f.write(text + "\n")
	f.close()

def find_available_data_dir_name(basename):
	current_name = basename + ""
	i = 1
	while os.path.isdir(current_name):
		current_name = basename + "_" + str(i)
		i+=1
	return current_name


to_log("iu")
def generate_simulation(ep,beta,mu,U,tpd,tppp,MEASUREMENT_TIME,Computing_days,iterations,params_type):
	path_to_main_dir = "../.."
	data_dir_name = "ep" + str(ep) + "_beta" + str(beta) + "_mu" + str(mu) + "_U" + str(U)
	if tpd != None:
		data_dir_name += "_tpd" + str(tpd)
	if tppp != None:
		data_dir_name += "_tppp" + str(tppp)
	backup_path = "./BACKUP_START"
	current_path = os.path.join(os.getcwd(), path_to_main_dir)
	datas_dir =  os.path.join(current_path,"DataDifferentparameters")
	create_dir_safely(datas_dir)
	files_path = os.path.join(datas_dir,data_dir_name)
	#We want to create a directory that does not exist for now (no overlap)
	files_path = find_available_data_dir_name(files_path)
	create_dir_safely(files_path)
	create_dir_safely(os.path.join(files_path,"DATA"))
	create_dir_safely(os.path.join(files_path,"IN"))
	create_dir_safely(os.path.join(files_path,"OUT"))
	shutil.copy(os.path.join(backup_path,"Hyb1.json"),os.path.join(files_path,"IN/Hyb1.json"))
	
	#From here we create the run.sh file to submit to slurm
	sh_origin = open(os.path.join(backup_path,"parameter_run.sh"),"r")
	lines = sh_origin.readlines()
	sh_origin.close()
	sh_destination = open(os.path.join(files_path,"run.sh"),"w")
	sh_destination.write('#!/bin/bash\n')
	sh_destination.write("#SBATCH --time=" + Computing_days + "-00:00:00\n")
	sh_destination.write('#SBATCH --job-name="' + str(ep) + str(beta) + str(mu) + '"\n')
	for l in lines:
		sh_destination.write(l)
	sh_destination.write("../launch.py " + current_path + " " + files_path + " " + str(iterations) + " 1")
	sh_destination.close()
	#End ceration run.sh

	#From here we create the json parameter file
	f = open(os.path.join(backup_path,"params0.meas - " + params_type + ".json"),"r")
	params_json = json.loads(f.read())
	f.close()
	params_json = params_json["Parameters"]
	params_json["U"] = U
	params_json["ep"] = ep
	params_json["beta"] = beta
	params_json["mu"] = mu
	if tpd != None:
		params_json["tpd"] = tpd
	if tppp != None:
		params_json["tppp"] = tppp
	params_json["MEASUREMENT_TIME"] = MEASUREMENT_TIME
	params_json["HYB"] = "Hyb1.json"
	f = open(os.path.join(files_path,"IN/params1.json"),"w")
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
	if(len(sys.argv) >= 9):
		ep = float(sys.argv[1])
		beta = float(sys.argv[2])
		mu = float(sys.argv[3])
		U = float(sys.argv[4])
		MEASUREMENT_TIME = int(sys.argv[5])
		Computing_days = sys.argv[6]
		iterations = int(sys.argv[7])
		params_type = sys.argv[8]
		tpd = tppp = None
		if(len(sys.argv) >= 10):
			tpd = float(sys.argv[9])
		if(len(sys.argv) >= 11):
			tppp = float(sys.argv[10])
		generate_simulation(ep,beta,mu,U,tpd,tppp,MEASUREMENT_TIME,Computing_days,iterations,params_type)
	else:
		to_log("Usage : new_simulation.py ep beta mu U MEASUREMENT_TIME Computing_days iterations [tpd tppp]")