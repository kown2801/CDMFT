#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.5.4/bin/python3.5

from subprocess import call
import sys
import datetime
import os
import time

max_number_of_tries = 15

def run(iterations_max,files_dir,iteration_start = 0):

	def to_log(string):
		f = open(os.path.join(files_dir,"logfile"),"a")
		f.write(string)
		f.close()
	#files_dir is here the path of the folder to launch from the main dir (example ComputedData/ep9.0_beta60.0_mu12.41_U12.0_tpd1.4_tppp1.0)
	path_to_main_dir = "../"
	os.chdir(os.path.join(os.path.dirname(os.path.realpath(__file__)))) #The working directory is now the scripts directory

	iterations=iteration_start
	autocoherence_dir = os.path.join(path_to_main_dir,"SelfConsistency")
	solver_dir = os.path.join(path_to_main_dir,"ImpuritySolver")
	input_dir = os.path.join(os.path.join(files_dir,"IN"),"")
	output_dir = os.path.join(os.path.join(files_dir,"OUT"),"")

	data_dir = os.path.join(os.path.join(files_dir,"DATA"),"")
	if iterations_max != -1:
		to_log("Let's go for " + str(iterations_max) + " iterations" + "\n")
	else:
		to_log("Let's go for some iterations" + "\n")

	if iterations == 0:
		#In case you need to start from scratch (not working for now)
		call([os.path.join(autocoherence_dir,"CDMFT"), output_dir, input_dir, data_dir, "params", str(iterations)])
		iterations+=1	
	def do_one_iteration(iteration_number):
		#In the general case
		to_log("begin iteration " + str(iteration_number)  + " at: " + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) + "\n")
		if not os.path.exists(input_dir + "params" + str(iteration_number) + ".json"): #Check if the file exists to prevent errors and retry the previous iteration
			to_log("There was an error, an input parameters file was not produced ("+ input_dir + "params" + str(iteration_number) + ".json" + "), Retrying the iteration before. Try number " + str(number_of_tries))
			return -1
		to_log("Calling the Impurity Solver")
		call(["srun", os.path.join(solver_dir,"IS"),input_dir, output_dir,"params" + str(iterations)])
		#If the srun call failed, we should retry this iteration.
		if not os.path.exists(os.path.join(output_dir,"params" + str(iterations) + ".meas.json")):
			return 0
		#We need to verify that there are no nan numbers (we replace them with 0s)
		out_file = open(os.path.join(output_dir,"params" + str(iterations) + ".meas.json"))
		data = out_file.read()
		out_file.close()
		if "nan" in data:
			to_log("A nan was present in the output file")
		data = data.replace("nan","0")
		out_file  = open(os.path.join(output_dir,"params" + str(iterations) + ".meas.json"),"w")
		out_file.write(data)
		out_file.close()
		#We can now go on with the cycle
		to_log("Calling the self-consistency")
		call([os.path.join(autocoherence_dir,"CDMFT"), output_dir, input_dir, data_dir, "params", str(iterations)])
		to_log("end iteration " + str(iterations)  + " at: " + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) + "\n")

		#We call the plaquette computations on another sbatch to not keep a lot of processors idle
		call(["sbatch", "run_occupation.sh", output_dir, input_dir, data_dir, "params", str(iterations)])
		#We compute the order parameter
		call(["./actions.py", "-a","order_parameter", "-f",os.path.basename(files_dir)])
		return 1
	number_of_tries = 0
	while (iterations_max == -1 or iterations <= iterations_max) and number_of_tries <= max_number_of_tries:
		return_code = do_one_iteration(iterations)
		iterations+=return_code
		if return_code != 1:
			number_of_tries += 1
			time.sleep(60)

	return 0


if __name__ == "__main__":
	if(len(sys.argv) >= 4):
		files_dir = sys.argv[1]
		iterations_max = int(sys.argv[2])
		iteration_start = int(sys.argv[3])
		run(iterations_max,files_dir,iteration_start)
	else:
		print("Usage : launch.py files_dir iterations_max iterations_start")
