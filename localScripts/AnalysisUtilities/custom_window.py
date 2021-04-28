import pandas
import tkinter
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
import os
import re
from subprocess import call
from threading import Thread
from AnalysisUtilities import distant_consts as CONSTS

all_data_folder = CONSTS.distant_data_dir
scripts_dir = CONSTS.distant_scripts_dir

def print_and_mean(data,ax,mean_over):
    n = len(data)
    if mean_over > n:
        print("Not enough data to mean over but proceeding anyway with a lower data count")
        mean_over = n
    X = np.linspace(mean_over,n-1,n-mean_over)
    panda_data = pandas.DataFrame(data=data)
    ax.plot(data,label="Raw Data")
    ax.plot(X,panda_data.rolling(mean_over).mean().values.flatten()[mean_over:],label="Rolling mean")

#This class is a simple graphical interface that allows you to handle one simulation folder. The jupyter notebook Fast Distant allows you to do that for all folders in one shot.

class Process_window:
    def __init__(self,folder,mean_over,c):
        self.folder = folder
        self.mean_over = mean_over
        self.connection = c
        self.not_all_occupation = False
        self.root = tkinter.Tk()
        self.root.wm_title(folder["name"])
        self.fig = Figure(figsize=(8, 9), dpi=100)
        self.ax1 = self.fig.add_subplot(311)
        self.ax2 = self.fig.add_subplot(312)
        self.ax3 = self.fig.add_subplot(313)
        self.plot_data()
        self.create_buttons()
        canvas = FigureCanvasTkAgg(self.fig, master=self.root)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, self.root)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        tkinter.mainloop()
        
    def plot_data(self):
        #Plotting the downloaded data into the graphical interface
        graph,sign_graph,N_graph,pn_graph = \
            self.folder["graph"],self.folder["sign_graph"],self.folder["N_graph"],self.folder["pn_graph"],
        try:
            if len(graph.shape) > 1:
                for i in range(graph.shape[1]):
                    #We need to print multiple data to observe convergence
                    if graph[0,i]!= 0 and graph[0,0]:
                        #If the starting values are not 0, we can simply normalize
                        print_and_mean(graph[:,i]/graph[0,i]*graph[0,0],self.ax1,self.mean_over)
                        #Else we just want to print lines that coincide on the first point to have somewhat nice plots
                    elif graph[0,i] != 0:
                        print_and_mean(graph[:,i] - graph[0,i],self.ax1,self.mean_over)
                    else:
                        print_and_mean(graph[:,i] + graph[0,0],self.ax1,self.mean_over)
            else:
                print_and_mean(graph,self.ax1,self.mean_over)
        
            self.ax1.set_title("Data analysis")
            print_and_mean(sign_graph,self.ax2,self.mean_over)
            self.ax2.set_title("Sign")
            self.ax3.set_title("Total occupation")
            if len(N_graph) != len(pn_graph):
                if len(N_graph) > len(pn_graph):
                    print("L'occupation n'a pas été calculé en entier :", len(N_graph), len(pn_graph))
                    self.not_all_occupation = True
                    N_graph = N_graph[:len(pn_graph)]
                else:
                    self.occs()
                    messagebox.showwarning("Error in the Occupations", "For this simulation folder (" + self.folder["name"] + "), there was too much oxygen occupation data. The eronous data will be diplayed here for reference and a correction job has just been launched")
                    pn_graph = pn_graph[:len(N_graph)]
            print_and_mean(2*N_graph + 4*pn_graph,self.ax3,self.mean_over)
            self.ax1.legend()
            self.fig.subplots_adjust(top=0.95)
        except Exception as e:
            print("Not enough data in " + self.folder["name"] + " : " + str(len(N_graph)))
            print("Error message : " + str(e))
            return
    
    
    def create_one_button(self,text,function,frame):
        #Creates one button on the interface
        button = tkinter.Button(master=frame, text=text)
        button.config(command=lambda : self.button_clicked(button,function))
        button.pack(side=tkinter.LEFT)
        return button
    
    def create_buttons(self):
        #Create all the interface buttons, inputs and on click functions
        self.buttonframe = tkinter.Frame(self.root)
        button_cancel = self.create_one_button("Cancel Simulation",self.cancel,self.buttonframe)
        button_copy = self.create_one_button("Copy to computer",self.copy,self.buttonframe)
        button_occs = self.create_one_button("Run all occupations",self.occs,self.buttonframe)
        if self.not_all_occupation:
            button_occs.configure(bg="orange")
        button_delete = self.create_one_button("Delete last iteration",self.delete,self.buttonframe)
        button = tkinter.Button(master=self.buttonframe, text="Next",command=self._quit)
        button.pack(side=tkinter.LEFT)
        
        self.resumeframe = tkinter.Frame(self.root)
        self.iterations_input = tkinter.Entry(master=self.resumeframe)
        self.iterations_input.pack(side=tkinter.LEFT)
        self.iterations_input.insert(tkinter.END, str(len(self.folder["graph"]) + 35))
        self.time_input = tkinter.Entry(master=self.resumeframe)
        self.time_input.pack(side=tkinter.LEFT)
        self.time_input.insert(tkinter.END,"1")
        self.button_resume = self.create_one_button("Resume Simulation",self.resume_according_inputs,self.resumeframe)
        self.create_one_button("+1 day",self.resume_one_day,self.resumeframe)
        self.create_one_button("+2 days",self.resume_two_days,self.resumeframe)
        
        self.resumeframe.pack(side=tkinter.BOTTOM)
        self.buttonframe.pack(side=tkinter.BOTTOM)
        
    
    def _quit(self):
        #Stops the graphical interface 
        self.root.quit()     # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate
      
    def button_clicked(self,button,action):
        #Handles the clicked buttons and changes the button color in case the operation fails
        try:
            action()
            button.configure(bg = "green")
            button.configure(activebackground='#90ee90')
        except Exception as e:
            button.configure(bg = "red")
            button.configure(activebackground='#ffcccb')
            print(str(e))
    
    def cancel(self):
        #Stopping the job for the folder
        all_files = self.connection.run("cd " + os.path.join(all_data_folder,self.folder["name"]) + ";ls",hide=True).stdout.split()
        all_job_numbers = []
        for file in all_files:
            number = re.findall("slurm\-([0-9]*)\.out", file, flags=0)
            if number:
                all_job_numbers.append(int(number[0]))
        try:
            current_job = max(all_job_numbers)
            result = self.connection.run("scancel " + str(current_job),hide=True)
            print(result.stdout)  
        except Exception as e:
            print(e)
            print("It is very likely that there is no more slurm files in the folder")

    def copy(self):
        #Copying the file on the computer and removing the now useless files. It fails if the N and pn files don't have same dimensions
        try:
            self.cancel()
        except:
            pass
        self.connection.run(CONSTS.launch_numpy_actions + " -a prepare_copy -f" + self.folder["name"],hide=True)
        if self.folder["N_graph"].shape != self.folder["pn_graph"].shape:
            raise Exception("N and pn graphs should have the same shape when copying. Either run all occupations or check the simulation folder to correct this problem")
        def go_host_key():
            self.folder["name"]
            call(["scp","-rC",CONSTS.ssh_address + ":" + os.path.join(all_data_folder,self.folder["name"]) ,os.path.join(os.getcwd(),"AllData/transfered")])
            f = open("AllData/transfer.done","a")
            f.write(self.folder["name"] + "\n")
            f.close()
            local_transfered_dir = os.path.join("AllData/transfered/",self.folder["name"])
            if os.path.isdir(local_transfered_dir) and\
            os.path.isdir(os.path.join(local_transfered_dir,"DATA")) and\
            os.path.isdir(os.path.join(local_transfered_dir,"OUT")) and\
            os.path.isdir(os.path.join(local_transfered_dir,"IN")):
                self.connection.run("rm -r " + os.path.join(all_data_folder,self.folder["name"]))
        Thread(target=go_host_key).start()
        
    def occs(self):
        #Sends the signal to compute all occupations for the current folder
        self.connection.run(CONSTS.launch_sbatch_actions + " -a occupations -f " + self.folder["name"],hide=True)
    def delete(self):
        #Sends the signal to delete the last iteration
        self.connection.run(CONSTS.launch_numpy_actions + " -a delete_last -f " + self.folder["name"],hide=False)
        
        
    def resume_custom(self):#Sends the signal to resume the simulation
        try:
            self.cancel()
        except:
            pass
        result = self.connection.run("cd " + scripts_dir + ";./resume_simulation.sh "+ \
                       self.folder["name"] + " -1 " + str(self.more_iterations) + " " + str(self.more_time),hide=True)
        print(result.stdout)
    
    def resume_according_inputs(self):
        self.more_iterations = self.iterations_input.get()
        self.more_time = self.time_input.get()
        self.resume_custom()    
        
    def resume_one_day(self):
        self.more_iterations = -1
        self.more_time = 1
        self.resume_custom()
        
    def resume_two_days(self):
        self.more_iterations = -1
        self.more_time = 2
        self.resume_custom()
        
        