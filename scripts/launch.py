#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.5.4/bin/python3.5

from subprocess import call
import sys
import datetime
import os


def run(iterations_max,files_dir,iteration_start = 0):

	def to_log(string):
		f = open(os.path.join(files_dir,"logfile"),"a")
		f.write(string)
		f.close()
	#files_dir is here the path of the folder to launch from the main dir (example ComputedData/ep9.0_beta60.0_mu12.41_U12.0_tpd1.4_tppp1.0)
	path_to_main_dir = "../"
	os.chdir(os.path.join(os.path.dirname(os.path.realpath(__file__)))) #The working directory is now the scripts directory
	iterations=iteration_start
	autocoherence_dir = os.path.join(path_to_main_dir,"Autocoherence")
	solver_dir = os.path.join(path_to_main_dir,"ImpuritySolver")
	input_dir = os.path.join(os.path.join(files_dir,"IN"),"")
	output_dir = os.path.join(os.path.join(files_dir,"OUT"),"")

	data_dir = os.path.join(os.path.join(files_dir,"DATA"),"")
	to_log("Let's go for " + str(iterations_max) + " iterations" + "\n")
	if iterations == 0:
		#In case you need to start from scratch (not working for now)
		call([os.path.join(autocoherence_dir,"CDMFT"), output_dir, input_dir, data_dir, "params", str(iterations)])
		iterations+=1
	while iterations <= iterations_max:
		#In the general case
		to_log("begin iteration " + str(iterations)  + " at: " + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) + "\n")
		if not os.path.exists(input_dir + "params" + str(iterations) + ".json"): #Check if the files exists to stop the program (it doesn't stop otherwise)
			to_log("There was an error, an input parameters file was not produced ("+ input_dir + "params" + str(iterations) + ".json" + ")")
			exit(0)
		to_log("Calling the Impurity Solver")
		call(["srun", os.path.join(solver_dir,"IS"),input_dir, output_dir,"params" + str(iterations)])
		to_log("Calling the Autocoherence")
		call([os.path.join(autocoherence_dir,"CDMFT"), output_dir, input_dir, data_dir, "params", str(iterations)])
		to_log("end iteration " + str(iterations)  + " at: " + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) + "\n")

		#We call the plaquette computations on another sbatch to not keep a lot of processors idle
		call(["sbatch", "run_occupation.sh", output_dir, data_dir, "params", str(iterations)])
		#We compute the order parameter
		call(["./actions.py", "-a","order_parameter", "-f",os.path.basename(files_dir)])
		iterations +=1
	return 0


if __name__ == "__main__":
	if(len(sys.argv) >= 4):
		files_dir = sys.argv[1]
		iterations_max = int(sys.argv[2])
		iteration_start = int(sys.argv[3])
		run(iterations_max,files_dir,iteration_start)
	else:
		print("Usage : launch.py files_dir iterations_max iterations_start")
