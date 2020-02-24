#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.5.4/bin/python3.5
import launch
import os
import json
import shutil
import sys
def multiple_energies(energy_range,iterations):
	current_path = os.getcwd() 
	try:  
	    os.mkdir(current_path + "/DataDifferentEnergies")
	except OSError:  
		launch.write_to_log_file("Creation of the directory " + current_path + "/DataDifferentEnergies failed\n")
		if not os.path.exists(current_path + "/DataDifferentEnergies"):
			launch.write_to_log_file("Sorry, we need this directory to proceed, stopping... \n")
			exit()
	else:  
	    launch.write_to_log_file("Successfully created the directory " + current_path + "/DataDifferentEnergies \n")
	for i,energy in enumerate(energy_range):
		try:  
			files_path = current_path + "/DataDifferentEnergies/Energy" + str(i)
			os.mkdir(files_path)
			os.mkdir(files_path + "/DATA")
			os.mkdir(files_path + "/IN")
			os.mkdir(files_path + "/OUT")
		except OSError:  
			launch.write_to_log_file("Creation of the directories in " + current_path + "/DataDifferentEnergies/Energy" + str(i) + " failed\n")
			if not os.path.exists(current_path + "/DataDifferentEnergies/Energy" + str(i)):
				launch.write_to_log_file("Sorry, we need this directory to proceed, stopping... \n")
				exit()
		else:  
			launch.write_to_log_file("Successfully created the directory  " + current_path + "/DataDifferentEnergies/Energy" + str(i) + "\n")
		f = open(current_path + "/BACKUP_START/params0.meas.json","r")
		params_json = json.loads(f.read())
		f.close()
		params_json["Parameters"]["ep"] = energy
		f = open(files_path + "/OUT/params0.meas.json","w")
		shutil.copy(current_path + "/BACKUP_START/LinkA.json",files_path + "/IN/LinkA.json")
		shutil.copy(current_path + "/BACKUP_START/LinkN.json",files_path + "/IN/LinkN.json")
		f.write(json.dumps(params_json, indent=4,sort_keys=True))
		f.close()
		launch.run(iterations,current_path,files_path)
		#Then launch the simulation with the right directories



if __name__ == "__main__":
	if(len(sys.argv) >= 5):
		start = int(sys.argv[1])
		stop = int(sys.argv[2])
		step = int(sys.argv[3])
		iteration = int(sys.argv[4])
		energy_range = range(start,stop,step)
	else:
		launch.write_to_log_file("Usage : different_energies.py start stop step iterations_number")
		launch.write_to_log_file("But that's okay for this time")
		energy_range = range(0,10,2)
		iteration = 10
	multiple_energies(energy_range,iteration)