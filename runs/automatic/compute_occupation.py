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

def get_last_computed_meas(folder):
	i = 1
	try:
		while 1:
			f = open(os.path.join(folder,"OUT/params" + str(i) + ".meas.json"),"r")
			i+=1
			f.close() 
	except Exception as e:
		pass
	print(str(i-1) + " is supposed to be the highest iteration")
	return i-1

def one_varying_parameter(dir_name,iteration_min,iteration_max,Computing_time):
	path_to_main_dir = "../.."
	backup_path = "./BACKUP_START"
	supra_path = os.path.join(os.getcwd(), path_to_main_dir)
	files_path = supra_path + "/DataDifferentparameters/" + dir_name
	if iteration_max == None:
		iteration_min = "1"
		iteration_max = get_last_computed_meas(files_path)
		Computing_time = str(iteration_max//33 + 1)
		iteration_max = str(iteration_max)

	#From here we create the run.sh file to submit to slurm
	sh_origin = open(os.path.join(backup_path,"occupation_run.sh"),"r")
	lines = sh_origin.readlines()
	sh_origin.close()
	sh_destination = open(os.path.join(files_path,"run_occupation.sh"),"w")
	sh_destination.write('#!/bin/bash\n')
	sh_destination.write("#SBATCH --time=0" + Computing_time + ":00:00\n")
	sh_destination.write('#SBATCH --job-name="dmft_occupation"\n')
	for l in lines:
		sh_destination.write(l)
	sh_destination.write("../launch_occupation.py " + supra_path + " " + files_path + " " + str(iteration_min) + " " + str(iteration_max))
	sh_destination.close()
	#End creation run.sh
	#Begin start simulation
	os.chdir(files_path)
	subprocess.run(["sbatch","run_occupation.sh"])
	#End start simulation

if __name__ == "__main__":
	if(len(sys.argv) >= 5):
		folder = sys.argv[1]
		iteration_min = sys.argv[2]
		iteration_max = sys.argv[3]
		Computing_time = sys.argv[4]
	elif(len(sys.argv) >= 2):
		folder = sys.argv[1]
		iteration_min = None
		iteration_max = None
		Computing_time = None
	else:
		print("Usage : compute_occupation.py folder iteration_min iteration_max")
		exit(0)
	one_varying_parameter(folder,iteration_min,iteration_max,Computing_time)