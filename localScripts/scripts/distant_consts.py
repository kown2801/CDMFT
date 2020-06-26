import os
################ Variables you CAN and SHOULD change ##################
username = "kowalski"
domain = "cedar.computecanada.ca"
distant_main_dir = "/scratch/kowalski/FullCDMFT/git-version/CDMFT/"


########## You normally don't need to touch those variables #########
ssh_address = username + "@" + domain
data_dir_from_main = "ComputedData"
scripts_dir_from_main = "scripts"
distant_data_dir = os.path.join(distant_main_dir,data_dir_from_main)
distant_scripts_dir = os.path.join(distant_main_dir,scripts_dir_from_main)
launch_numpy_actions = "cd " + distant_scripts_dir + ";" + "./actions.sh" + " "
launch_sbatch_actions = "cd " + distant_scripts_dir + ";sbatch " + "actions.sh" + " "