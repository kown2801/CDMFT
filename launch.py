#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.5.4/bin/python3.5

from subprocess import call
import sys
import datetime


def run(iterations_max,main_dir,files_dir):
	iterations=0

	autocoherence_dir = main_dir + "/Autocoherence/"
	solver_dir = main_dir + "/ImpuritySolver/"
	input_dir = files_dir + "/IN/"
	output_dir = files_dir + "/OUT/"
	data_dir = files_dir + "/DATA/"
	write_to_log_file("Let's go for " + str(iterations_max) + " iterations" + "\n")
	if iterations == 0:
		call([autocoherence_dir + "CDMFT", output_dir, input_dir, data_dir, "params", str(iterations)])
		iterations+=1
	while iterations <= iterations_max:
		write_to_log_file("begin iteration " + str(iterations)  + " at: " + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) + "\n")
		call(["srun", solver_dir + "IS",input_dir, output_dir,"params" + str(iterations)])
		call([autocoherence_dir + "CDMFT", output_dir, input_dir, data_dir, "params", str(iterations)])
		call([autocoherence_dir + "GFULL", output_dir, data_dir, "params", str(iterations)])
		write_to_log_file("end iteration " + str(iterations)  + " at: " + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) + "\n")
		iterations +=1
	return 0

def write_to_log_file(string):
	f = open("logfile","a")
	f.write(string)
	f.close()



if __name__ == "__main__":
	if(len(sys.argv) >= 2):
		iteration = int(sys.argv[1])
	else:
		write_to_log_file("Usage : launchy.py iterations_number")
		write_to_log_file("But that's okay for this time")
		iteration = 10
	run(iteration,"/home/kowalski/scratch/FullCDMFT/Supra","/home/kowalski/scratch/FullCDMFT/Supra")
