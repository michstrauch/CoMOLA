# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#
#   Name:       filehandler.py
#   Purpose:    This module is responsible for file operations 
#
#   Author:     Carola Paetzold, Christian Schweitzer, Michael Strauch
#   Contact:    michael.strauch@ufz.de
#
#               Helmholtz Centre for Environmental Research - UFZ
#               Department Computational Landscape Ecology - CLE
#               Permoserstrasse 15
#               D-04318 Leipzig, Germany
#               http://www.ufz.de
#
#   Created:    Mo Apr 14 2014
#
#   Copyright:  (c) Carola Paetzold / Christian Schweitzer 2014
#
#   Licence:    This program is free software:
#               you can redistribute it and/or modify it under the terms
#               of the GNU General Public License as published by the
#               Free Software Foundation, either version 3 of the License,
#               or (at your option) any later version. This program is
#               distributed in the hope that it will be useful, but
#               WITHOUT ANY WARRANTY; without even the implied warranty
#               of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#               See the GNU General Public License for more details.
#               You should have received a copy of the GNU General
#               Public License along with this program.
#               If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------
import os
import time 
import pylab
import csv
import subprocess
from inspyred import ec
import shutil
import numpy as np
import sys

import requirements as req

wrkDir = os.path.abspath('.')
timestamp_file = time.strftime("%d-%m-%Y_%H-%M-%S_")
if os.path.exists(os.path.join(wrkDir, "output", "help_file.txt")):
    fobj = open(os.path.join(wrkDir, "output", "help_file.txt"), "r") 
    for line in fobj:
        timestamp_file = line         
    fobj.close()
    
# Array for worst fitness values if individuals are filtered of plausibility
worst_fitness = np.array([], dtype=np.float64)

#------------------------------------------------------------------------------
#   Log-file handling
#------------------------------------------------------------------------------
def InitLogFile():   
    """ Create an optimization_log file for the log information of the optimization."""
      
    fileName = timestamp_file + "optimization_log.txt"   
    logFile = open(os.path.join(wrkDir, "output", fileName), "w")
    logFile.close()
    # clear the helping folder for logging of parallel running processes
    if os.path.exists(os.path.join(wrkDir, "output", "child_processes")):
        for filename in os.listdir(os.path.join(wrkDir, "output", "child_processes")):
            os.remove(os.path.join(wrkDir, "output", "child_processes", filename))

# lines for printing the checked new_cand from logical_variator
def WriteCandidateList(candidate_list):
    """ Write checked candidates in a file.

        input:
            candidate_list is a list with checked candidates and informations
                           of the new individual generation with logical_variator
                           (maphandler.py)
    """

    fileName = timestamp_file + "candidate_list.txt"  
    if os.path.isfile(os.path.join(wrkDir, "output", fileName)):
        logFile = open(os.path.join(wrkDir, "output", fileName), "a")
    else:
        logFile = open(os.path.join(wrkDir, "output", fileName), "w") 
    for element in candidate_list:
        msg = "%s\n" %element 
        logFile.writelines(msg)
    msg = "\nnext\n" 
    logFile.writelines(msg)
    logFile.close()    

def WriteLogMsg(msg,ind_number=0):
    """ Write the message with special formatting in the optimization_log or ind_number file.

        input:
            msg contains the message for the log file
            ind_number is the individual number of the current population
    """
    
    if not os.path.exists(os.path.join(wrkDir, "output", "child_processes")):
       os.mkdir(os.path.join(wrkDir, "output", "child_processes"))    
    timestamp = time.strftime("%d-%m-%Y %H:%M:%S | ")
    if ind_number != 0:
        fileName = "ind_number" + str(ind_number) + ".txt"
        if os.path.isfile(os.path.join(wrkDir, "output", "child_processes", fileName)):
            logFile = open(os.path.join(wrkDir, "output", "child_processes", fileName), "a")
        else:
            logFile = open(os.path.join(wrkDir, "output", "child_processes", fileName), "w")
    else:
        fileName = timestamp_file + "optimization_log.txt"
        logFile = open(os.path.join(wrkDir, "output", fileName), "a") 
    msg = timestamp + msg + "\n"
    logFile.writelines(msg)
    logFile.close()
    print (msg)

def join_ind_number_log():
    """ Join the logs of the child processes in the optimization_log file."""
    
    fileName = timestamp_file + "optimization_log.txt"
    logFile = open(os.path.join(wrkDir, "output", fileName), "a")    
    dirList = os.listdir(os.path.join(wrkDir, "output", "child_processes"))
    dirList.sort()
    for file in dirList:
        fobj = open(os.path.join(wrkDir, "output", "child_processes", file), "r") 
        for line in fobj:
            logFile.writelines(line)         
        fobj.close()
        os.remove(os.path.join(wrkDir, "output", "child_processes", file))
    logFile.close()

    
def inspyred_logging(): 
    """ Create an log file for the information of the inspyred logger."""
     
    # inspyred logging for troubleshooting
    # http://inspyred.github.io/troubleshooting.html#use-and-consult-the-logs
    import logging
    logger = logging.getLogger('inspyred.ec')
    logger.setLevel(logging.DEBUG)
    #file_handler = logging.FileHandler('inspyred.log', mode='w')
    fileName = timestamp_file + 'inspyred.log'
    file_handler = logging.FileHandler(os.path.join(wrkDir, "output", fileName), mode='w')
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

#------------------------------------------------------------------------------
#   create output folder for log files
#------------------------------------------------------------------------------   
def create_output_folder():
    """ Create an output folder for saving the log files."""
    
    #timestamp = time.strftime("%d-%m-%Y_%H-%M-%S_")
    # create an output folder for the log-files
    try:
        os.makedirs(os.path.join(wrkDir,"output"))
    # output folder exists
    except:
        pass    

#------------------------------------------------------------------------------
#   save plot as pdf
#------------------------------------------------------------------------------   
def savePlot_pdf(algorithm):
    """ Save plot as pdf file.

        input:
            optimization algorithm
    """
    
    #timestamp = time.strftime("%d-%m-%Y_%H-%M-%S_")
    if algorithm == "NSGA2":
        fileName = timestamp_file + "NSGA2_result_plot.pdf"
    elif algorithm == "GA":
        fileName = timestamp_file + "GA_result_plot.pdf"
    pylab.savefig(os.path.join(wrkDir, "output", fileName),format='pdf')
    
#------------------------------------------------------------------------------
#   save plot as png
#------------------------------------------------------------------------------   
def savePlot_png(algorithm):
    """ Save plot as png file.

        input:
            optimization algorithm
    """
    
    #timestamp = time.strftime("%d-%m-%Y_%H-%M-%S_")
    if algorithm == "NSGA2":
        fileName = timestamp_file + "NSGA2_result_plot.png"
    elif algorithm == "GA":
        fileName = timestamp_file + "GA_result_plot.png"

    # for constrained_tournament_selection: create a second plot with feasible solutions
    if os.path.exists(os.path.join(wrkDir, "output", fileName)):
        if algorithm == "NSGA2":
            fileName = timestamp_file + "NSGA2_result_plot_feasible.png"
        elif algorithm == "GA":
            fileName = timestamp_file + "GA_result_plot_feasible.png"

    pylab.savefig(os.path.join(wrkDir, "output", fileName),format='png')

#------------------------------------------------------------------------------
#   plot with ec.analysis.generation_plot of the statistics_file.csv
#------------------------------------------------------------------------------ 
def plot_statistics_file():
    """ Plot a graphic with ec.analysis.generation_plot of the statistics_file.csv.""" 
    
    # Plot results - True/False = to use errorbars or not
    fileName = timestamp_file + "statistics_file.csv"
    ec.analysis.generation_plot(os.path.join(wrkDir, "output", fileName), True)
    
#------------------------------------------------------------------------------
#   Log-file handling for inspyred
#------------------------------------------------------------------------------
def init_inspyred_logfiles():
    """ Create log files for the inspyred internal log file options.
        
        The inspyred package has an internal log file possibility.
        To make use of this, file objects have to be passes to the 
        algorithms.
    """
    
    # Save the results from each generation to 'statistics_file'  
    fileName = timestamp_file + "statistics_file.csv"
    stats_file = open(os.path.join(wrkDir, "output", fileName), "w")   
    # plot only run without header 
    #header = "generation number, population size, worst, best,"
    #header += "median, average, standard deviation \n"    
    #stats_file.write(header)
         
    # Save all individual results from each run to 'individuals_file'
    fileName = timestamp_file + "individuals_file.csv"        
    individ_file = open(os.path.join(wrkDir, "output", fileName), "w")
    header = "generation number, individual number,"
    header += "fitness, string representation of candidate \n"
    individ_file.write(header)

    return stats_file, individ_file
 
#------------------------------------------------------------------------------
#   Documentation of the input data
#------------------------------------------------------------------------------
def save_input_data():
    """ Create an text file for saving the input data.
        
        input data: config.ini
    """
    
    #timestamp = time.strftime("%d-%m-%Y_%H-%M-%S_")
    fileName = timestamp_file + "input_data.txt"
    # import config.py
    fobj_config = open("config.ini", "r") 
    # open or create file for saving data
    fobj_docu = open(os.path.join(wrkDir, "output", fileName), "w")  
    for line in fobj_config: 
        fobj_docu.write(line + "\n") 
    fobj_config.close()
    fobj_docu.close()

#------------------------------------------------------------------------------
#   Documentation of the model console outputs in one output file
#------------------------------------------------------------------------------
def summarize_console_outputs(number_individuals, number_generation,individuals, external_models, not_accepted_ind):
    """ Summarize the console outputs of the models after 
        each population evaluation in one output file.
        
        input data: console.txt files
                    number of individuals
                    generation number
                    population list
                    individuals without model runs
    """
    
    fileName = timestamp_file + "model_outputs.txt"
    i=0
    fobj_docu = open(os.path.join(wrkDir, "output", fileName), "a+")  
    fobj_docu.seek(0, os.SEEK_END)
    fobj_docu.write("\n") 
    fobj_docu.write("generation %s" %number_generation + "\n\n") 
    while i < number_individuals:
        if i+1 not in not_accepted_ind:
            fobj_docu.write("population %s" %individuals[i] + "\n\n") 
            if i == 0:
                model_folder = 'models'
            else:
                model_folder = 'models_%s' %i
            subfolder_files = os.listdir(os.path.join(wrkDir, model_folder))
            subfolder_files.sort()
            for k in subfolder_files:
                # check if k is a folder
                if os.path.isdir(os.path.join(wrkDir, model_folder, k)):
                    for m in external_models:
                        # subfolder == possible model folder
                        if k == m:
                            try:
                                fobj_docu.write("model path %s" %os.path.join(wrkDir, model_folder, k) + "\n") 
                                fobj_txt = open(os.path.join(wrkDir, model_folder, k, 'console.txt'), "r")
                                for line in fobj_txt: 
                                    fobj_docu.write(line + "\n") 
                                fobj_txt.close()
                            except:
                                pass
        i += 1
    fobj_docu.close() 

#------------------------------------------------------------------------------  
#   Collect all the fitness values of one generations
#------------------------------------------------------------------------------
def collect_fitness_values(opt_algorithm, number_individuals, fitness, external_models, output_files, not_accepted_ind, file_worst_fitness):
    """Read the fitness values from external models, write them into the log file
       and append them to the fitness list.
 
       input data: optimization algorithm
                   number of individuals
                   list of fitness values
                   external models
                   output files  
                   individuals without model runs  
                   file with worst fitness values               
     """
  
    # Array for worst fitness values if individuals are filtered of plausibility
    global worst_fitness
    if worst_fitness.size == 0 and file_worst_fitness != 'None':
        worst_fitness = np.genfromtxt(os.path.join(wrkDir, 'input', file_worst_fitness), dtype=float)
    
    WriteLogMsg("not_accepted_ind %r" %not_accepted_ind)
     
    i=0
    count_worst_fitness = 0
    count_real_fitness = 0
    while i < number_individuals:
        fitness_model = []
        # individual break plausibility rules
        if i+1 in not_accepted_ind:
            # return worst fitness values for the excluded individual
            # one worst fitness value
            count_worst_fitness = 0
            if worst_fitness.size == 1:
                fitness_model.append(worst_fitness.item(0))
                count_worst_fitness += 1
            # more than one worst fitness value
            else:
                for n in worst_fitness:
                    fitness_model.append(n)
                    count_worst_fitness += 1
        else:
            if i == 0:
                model_folder = 'models'
            else:
                model_folder = 'models_%s' %i
            subfolder_files = os.listdir(os.path.join(wrkDir, model_folder))
            subfolder_files.sort()
            count_real_fitness = 0
            for k in subfolder_files:
                # check if k is a folder
                if os.path.isdir(os.path.join(wrkDir, model_folder, k)):
                    j=0
                    for m in external_models:
                        # subfolder == possible model folder
                        if k == m:
                            values = read_fitness_value(os.path.join(wrkDir, model_folder, k, output_files[j]))
                            for n in values:
                                fitness_model.append(n)
                                count_real_fitness += 1
                            #WriteLogMsg("Fitness value from external model %s: %s" % (os.path.join(wrkDir, model_folder, k),fitness_model))
    
                        else:
                            j += 1
                            
        if (count_real_fitness != 0 and count_worst_fitness != 0) and count_real_fitness != count_worst_fitness:
            msg = "Error: Number of worst fitness values (%s) is not equal with the number of the model fitness values (%s). Please check the worst fitness values." %(count_worst_fitness, count_real_fitness)
            WriteLogMsg(msg)
            raise SystemError("Error: Number of worst fitness values (%s) is not equal with the number of the model fitness values (%s). Please check the worst fitness values."%(count_worst_fitness, count_real_fitness))
            req.close_window
                            
        #WriteLogMsg("Fitness values: %s, %s" %(type(fitness_model),fitness_model))
        if opt_algorithm == "GA":
            fitness.append(fitness_model[0])                
        else:                                        
            # if NSGA 2 then multiple parameters are expected
            fitness.append(ec.emo.Pareto(fitness_model)) 
        #print("fitness-werte bei model-Ordner wechsel %s" %fitness)
        i += 1
        
    return fitness

#------------------------------------------------------------------------------  
#   Change parameters in parameter file
#------------------------------------------------------------------------------
def change_parameter_values(genom, ind_number):    
    """Update genom in csv file in the models folders chosen by optimization algorithm.

       input data:
            genom holds the individual
            ind_number - describe in which model folder the map should be saved
                         or if best values are saved in the output folder than 
                         the number attributed the individual number to the map
		
    """
   
    if (ind_number == 1):
        folder_name = 'models'
    else:
        folder_name = 'models_%s' %(ind_number-1)
    
    msg = "Update genome in models folder %s" %(ind_number)
    WriteLogMsg(msg)
    
    for folder in os.listdir(os.path.join(wrkDir,folder_name)):
        #msg = "Update genome %r in %s" %(genom,os.path.join(wrkDir,folder_name,folder,'genom.csv')) 
        if os.path.isdir(os.path.join(wrkDir,folder_name,folder)):
            in_file = open(os.path.join(wrkDir,folder_name,folder,'genom.csv'), "w")
            in_file.write("genom\n")
            for i in genom:
                line = "%d\n" % i
                in_file.write(line)  
            # Close file
            in_file.close()
        else:
            msg = "Error in executing external model %s. Directory does not exist." % os.path.join(wrkDir,folder_name,folder)
            WriteLogMsg(msg,ind_number)
            raise SystemError("Error in executing external model. Directory does not exist.")
            req.close_window

#------------------------------------------------------------------------------  
#   Documentation of the best solutions
#------------------------------------------------------------------------------
def save_best_solutions(best_solutions, dimension):    
    """Save the best solutions of the optimization process in an csv file.

       input data:
            best_solutions holds the best individuals and their fitness values
            dimension is the number of implemented models
		
    """
    
    fileName = timestamp_file + "best_solutions.csv"

    # create two files for constraint_tournament_selection 
    if os.path.exists(os.path.join(wrkDir, "output", fileName)):
        fileName = timestamp_file + "best_feasible_solutions.csv"

    in_file = open(os.path.join(wrkDir, "output", fileName), "w")
    in_file.write("individual,fitness values\n")
 
    for element in best_solutions:
        if dimension == 1:
            in_file.write("%s,[%f]\n" % (element.candidate, element.fitness))
        elif dimension == 2:
            in_file.write("%s,[%f,%f]\n" % (element.candidate, element.fitness[0], element.fitness[1])) 
        elif dimension == 3:
            in_file.write("%s,[%f,%f,%f]\n" % (element.candidate, element.fitness[0], element.fitness[1], 
                element.fitness[2])) 
        else:
            in_file.write("%s,[%f,%f,%f,%f]\n" % (element.candidate, element.fitness[0], element.fitness[1], 
                element.fitness[2], element.fitness[3])) 
    # Close file
    in_file.close()
        
#------------------------------------------------------------------------------  
#   Write map in ascii file
#------------------------------------------------------------------------------        
def WriteMap(header_all, map_info, modelfolder=False, ind_number=0, feasible=True):
    """Write the header and map information in an ascii file.
    
       input data: 
            header holds the original header information
            map_info holds the map information
            modelfolder - if True than write the ascii-map in the modelfolder
                          else write it in the output-folder
            ind_number - describe in which model folder the map should be saved
                         or if best values are saved in the output folder than 
                         the number attributed the individual number to the map
            information for constraint_tournament_selection if individual is feasible
    """  
    
    # update maps in model folders    
    if ind_number != 0 and modelfolder == True:
        if ind_number == 1:
            folder_name = 'models'
        else:
            folder_name = 'models_%s' %(ind_number-1)
 
        for folder in os.listdir(os.path.join(wrkDir,folder_name)):
            in_file = open(os.path.join(wrkDir,folder_name,folder,'map.asc'), 'w')
            for element in header_all:
                if type(element) == bytes:
                    in_file.write(element.decode(encoding='UTF-8')) 
                else:
                    in_file.write(element)
            for row in map_info:
                for item in row:
                    in_file.write("%s " % item)  
                in_file.write("\n")
            # Close file
            in_file.close()
    else:
        # write best maps in output folder
        if ind_number != 0 and feasible == True: 
            fileName = timestamp_file + "best_ascii_map%s.asc" %ind_number
        # for constraint_tournament_selection: write best maps in output folder 
        #   with information if individual is feasible
        elif ind_number != 0 and feasible == False:
            fileName = timestamp_file + "best_ascii_map%s_infeasible.asc" %ind_number
        # write ID-patch-map in output folder
        else:
            fileName = timestamp_file + "patch_ID_map.asc"
        
        #WriteLogMsg("Update map in %s" %os.path.join(wrkDir, "output", fileName))
        in_file = open(os.path.join(wrkDir, "output", fileName), "w")
        for element in header_all:
            if type(element) == bytes:
                in_file.write(element.decode(encoding='UTF-8')) 
            else:
                in_file.write(element)
        for row in map_info:
            for item in row:
                in_file.write("%s " % item)  
            in_file.write("\n") 
        # Close file
        in_file.close()

#-------------------------------------------------------------------------------------  
#   Read fitness values from external models and append them on global fitness vector
#-------------------------------------------------------------------------------------
def read_fitness_value(file):
    """Read the fitness value from external model and append it on global fitness vector.

       input:
           output file path of a model 
    """
            
    # Open file and read value     
    # in_file is an instance of lists with strings
    fitness = []
    in_file = csv.reader(open(file, "r"))
    
    # transfer the list contents from in_file in the fitness list (one list with floats)
    for row in in_file: 
        fitness.append(float(row[0])) 
    
    # fitness holds the fitness values of the external models
    return fitness
#-------------------------------------------------------------------------------------  
#   Run external model
#-------------------------------------------------------------------------------------
def run_model(file_path, file_path_R, file_path_python, RPy2, number):
    """Run the external model from file_path.

       input:
           file_path is the script file path of a model
           file_path_R is the path for the R file
           file_path_python is the path for the python file
           RPy2 information is RPy2 is available
           number is the individual number of the current population
    """

    # check if model script, R path and python path exist
    # if not then break the program with an error message
    if not os.path.isfile(file_path):
        msg = "Error in executing external model %s. File does not exist." % file_path
        WriteLogMsg(msg, number)
        raise SystemError("Error in executing external model. File does not exist. %s" % file_path)
        req.close_window
    elif not os.path.isfile(file_path_R):
        msg = "Error in executing external model %s. File does not exist." % file_path_R
        WriteLogMsg(msg, number)
        raise SystemError("Error in executing external model. File does not exist. %s" % file_path_R)
        req.close_window
    elif not os.path.isfile(file_path_python):
        msg = "Error in executing external model %s. File does not exist." % file_path_python
        WriteLogMsg(msg, number)
        raise SystemError("Error in executing external model. File does not exist. %s" % file_path_python)
        req.close_window
    
    # Split the filename, ext holds the file extension
    name, ext = os.path.splitext(os.path.basename(file_path))
    folder, file = os.path.split(file_path)
    #print("name, extension: %s, %s" %(name,ext))
    
    if ext == ".R" and RPy2 == "True":
        
        # delete lines for setting a new work space in the R file
        infile = open(file_path,'r')
        newopen = open(os.path.join(os.path.dirname(file_path), 'help_script.txt'), 'w')
        for line in infile :
            if ('setwd' not in line) or ('setwd' in line and '#' in line ):
                    newopen.write(line)
        newopen.close()
        infile.close()
        os.remove(file_path)
        os.rename(os.path.join(os.path.dirname(file_path), 'help_script.txt'), file_path)
        # append setwd for setting the R workspace on the right path at the beginning of the file
        setwd_path = r'%s'%os.path.dirname(r'%s'%file_path)
        setwd_path = setwd_path.replace('\\', '/')
        with open(file_path, "r+") as f: s = f.read(); f.seek(0); f.write('setwd("%s")\n'%setwd_path + s)
        f.close
        
        # Import for run R model via RPy2
        import rpy2.robjects as robjects
        begin=time.time() 
        r=robjects.r
        # redirect console outputs to a logFile
        fileName = 'console' + '.txt'
        logFile = os.path.join(folder, fileName)  
        r.sink(logFile, append=False)
        # open R script and setwd on the current path of the script
        r.source(file_path, chdir = True)
        end=time.time()
        msg = "The model %s ran for %d seconds." %(file_path, end-begin)
        WriteLogMsg(msg, number)
        r.sink()
          
    elif ext == ".R" and RPy2 == "False":
        
        fileName = 'console' + '.txt'
        logFile = os.path.join(folder, fileName)
        logFile = logFile.replace('\\', '/') 
        
        # delete lines for setting a new work space in the R file
        infile = open(file_path,'r')
        newopen = open(os.path.join(os.path.dirname(file_path), 'help_script.txt'), 'w')
        for line in infile :
            if (('setwd' not in line) or ('setwd' in line and '#' in line )) and ('sink' not in line):
                    newopen.write(line)
        newopen.close()
        infile.close()
        os.remove(file_path)
        os.rename(os.path.join(os.path.dirname(file_path), 'help_script.txt'), file_path)
        # append setwd for setting the R workspace on the right path at the beginning of the file
        setwd_path = r'%s'%os.path.dirname(r'%s'%file_path)
        setwd_path = setwd_path.replace('\\', '/')
        # without redirect the output of the R file in the logFile
        #with open(file_path, "r+") as f: s = f.read(); f.seek(0); f.write('setwd("%s")\n'%setwd_path + s)
        # redirect the output of the R file in the logFile
        with open(file_path, "r+") as f: s = f.read(); f.seek(0); f.write('setwd("%s")\nsink("%s", append=FALSE)\n'%(setwd_path,logFile) + s +'\nsink()\n')
        f.close 
    
        # Execute external model
        try:

            # this line create an output file with all infile and output lines
            #cmd = '"%s" --vanilla <%s >%s'%(file_path_R,file_path,os.path.join(folder, fileName))
          
            # or only the output lines will be documented
            # via append sink() in the R file 
            cmd = [file_path_R, 'CMD', 'BATCH', file_path]  
            begin=time.time()
            subprocess.check_call(cmd)
            end=time.time()
            msg = "The model %s ran for %d seconds." %(file_path, end-begin)
            WriteLogMsg(msg, number)  
            
            # check fitness values of the models
            #fitness = []
            #read_fitness_value(os.path.join(os.path.dirname(file_path),'R_output.csv'), fitness)
            #print ("fitnesswerte von %s" %os.path.join(os.path.dirname(file_path),'R_output.csv'))
            #for i in fitness:
                #print i

        except:
            exctype, value = sys.exc_info()[:2]
            msg = "Error in executing external model %s. Error number %s and error message: %s." % (file_path,exctype,value)
            WriteLogMsg(msg, number)
            raise SystemError("Error in executing external model.")
            req.close_window

    elif ext == ".py":
        cmd = [file_path_python, file_path]
       
        # Execute external model
        try:
            fileName = 'console' + '.txt'
            logFile = open(os.path.join(folder, fileName), 'w')
            begin=time.time()    
            subprocess.check_call(cmd, stdout=logFile)
            end=time.time()
            msg = "The model %s ran for %d seconds." %(file_path, end-begin)
            WriteLogMsg(msg, number)
            logFile.close()

        except:
            exctype, value = sys.exc_info()[:2]
            msg = "Error in executing external model %s. Error number %s and error message: %s." % (file_path,exctype,value)
            WriteLogMsg(msg, number)
            raise SystemError("Error in executing external model.")
            req.close_window
        
    else:
        msg = "It is no call for this file extension implemented %s." % file_path
        WriteLogMsg(msg, number)
        raise SystemError("It is no call for this file extension implemented.")
        req.close_window
    
#-------------------------------------------------------------------------------------  
#   Prepare string value in an executable attribute
#-------------------------------------------------------------------------------------
def preparing_attribute(type, attribute):
    """Prepare string value in an executable attribute.

       input:
           type is the inspyred operator type for example variator, terminator etc.
           attribute is the config.ini setting of the operator
    """
            
    # split attribute     
    list = attribute.split(',')
    # list with one argument
    if len(list) == 1:
        arg = "ec.%ss.%s" % (type, list[0])
    else:
        arg = "["
        for element in list:
            arg += "ec.%ss.%s" % (type, element)
            if element != list[-1]:
                arg += ", "
        arg += "]"
            
    return arg

#-------------------------------------------------------------------------------------  
#   Write filestamp_file of the main process in a file, so the multiprocessing 
#   processes can write in the same log file
#-------------------------------------------------------------------------------------
def save_timestamp(timestamp_file_parents):
    """Write the timestamp of the timestamp_file from main process in a file.

       input:
           timestamp for the created files in the output folder
    """
            
    fileName = "help_file.txt"   
    helpFile = open(os.path.join(wrkDir, "output", fileName), "w") 
    msg = timestamp_file
    helpFile.writelines(msg)
    helpFile.close()
    
#-------------------------------------------------------------------------------------  
#   Check if enough copies of the models folder exists for the multiprocessing.
#   If not than create the missing folders.
#-------------------------------------------------------------------------------------
def copy_models(number_genom):
    """Check if enough copies of the models folder exists for the multiprocessing 
       to run the models for every genom parallel. If not than create the missing folders.

       input:
           number_genom is the population size
    """
       
    for i in range(1,number_genom):
        # check if the copy i of the models folder exist
        folder_name = 'models_%s' %i
        # if not than copy the models folder with the new name
        if not os.path.exists(os.path.join(wrkDir,folder_name)):
            msg = "Create the helping folder ..."
            WriteLogMsg(msg) 
            shutil.copytree(os.path.join(wrkDir,'models'),os.path.join(wrkDir,folder_name))
            msg = "Done."
            WriteLogMsg(msg) 
    
#-------------------------------------------------------------------------------------  
#   Delete the helping models folder
#-------------------------------------------------------------------------------------
def delete_models():
    """Delete the helping models folder before termination of the optimization tool."""
    
    msg = "Delete the helping folders ..."
    WriteLogMsg(msg)    
    i = 1
    folder_name = 'models_%s' %i
    while os.path.exists(os.path.join(wrkDir,folder_name)):
        shutil.rmtree(os.path.join(wrkDir,folder_name))
        i += 1
        folder_name = 'models_%s' %i
    msg = "Done."
    WriteLogMsg(msg) 

#------------------------------------------------------------------------------  
#   Test functions
#------------------------------------------------------------------------------   
#if __name__ == "__main__":            
    #test functions
    
#------------------------------------------------------------------------------
#
#   EOF
#
#------------------------------------------------------------------------------
