#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.8.0/bin/python3.8
import os
import json
import shutil
import sys
import subprocess
import re
import glob

def create_dir_safely(dir):
	try:  
	    os.mkdir(dir)
	except OSError:  
		if not os.path.exists(dir):
			launch.write_to_log_file("Creation of the directory " + dir + " failed\n")
			launch.write_to_log_file("Sorry, we need this directory to proceed, stopping... \n")
			exit()
	to_log("Successfully created the directory " + dir + " \n")

def get_float(text):
	return re.findall("\d+\.\d+", text)[0]
def to_log(text):
	f = open("logfile","a")
	f.write(text + "\n")
	f.close()

def resume_simulation(data_dir_name,iteration_min,iteration_max,Computing_days):
	path_to_main_dir = "../.."
	backup_path = "./BACKUP_START"
	current_path = os.path.join(os.getcwd(), path_to_main_dir)
	datas_dir =  os.path.join(current_path,"DataDifferentparameters")
	
	files_path = os.path.join(datas_dir,data_dir_name)
	if iteration_min == "-1":
		all_iterations = []
		for file in glob.glob(os.path.join(files_path,"IN/params") + "*"):
			number = re.findall("params([0-9]*)\.json", file, flags=0)
			if number:
				all_iterations.append(int(number[0]))
		iteration_min = max(all_iterations)
	#From here we create the run.sh file to submit to slurm
	sh_origin = open(os.path.join(backup_path,"parameter_run.sh"),"r")
	lines = sh_origin.readlines()
	sh_origin.close()
	sh_destination = open(os.path.join(files_path,"run.sh"),"w")
	sh_destination.write('#!/bin/bash\n')
	sh_destination.write("#SBATCH --time=" + Computing_days + "-00:00:00\n")
	parameters = data_dir_name.split("_")
	ep = get_float(parameters[0])
	beta = get_float(parameters[1])
	mu = get_float(parameters[2])
	sh_destination.write('#SBATCH --job-name=""' + ep + beta + mu + '""\n')
	for l in lines:
		sh_destination.write(l)
	sh_destination.write("../launch.py " + current_path + " " + files_path + " " + str(iteration_max) + " " + str(iteration_min))
	sh_destination.close()
	#End ceration run.sh

	#Begin start simulation
	os.chdir(files_path)
	subprocess.run(["sbatch","run.sh"])
	#End start simulation


if __name__ == "__main__":
	if(len(sys.argv) >= 5):
		folder = sys.argv[1]
		iteration_min = sys.argv[2]
		iteration_max = sys.argv[3]
		Computing_days = sys.argv[4]
		resume_simulation(folder,iteration_min,iteration_max,Computing_days)
	else:
		to_log("Usage : resume_simulation.py folder iteration_min iteration_max Computing_days")