# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#
#   Name:       optiAlgorithm.py
#   Purpose:    This module provides the optimization algorithms.
#
#   Author:     Carola Paetzold, Christian Schweitzer, Michael Strauch
#   Contact:    michael.strauch@ufz.de

#               Helmholtz Centre for Environmental Research - UFZ
#               Department Computational Landscape Ecology - CLE
#               Permoserstrasse 15
#               D-04318 Leipzig, Germany
#               http://www.ufz.de
#
#   Created:    Mo Apr 14 2014
#
#   Copyright:  (c) Carola Paetzold / Christian Schweitzer / Michael Strauch 2018
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
#   Imports
#------------------------------------------------------------------------------
import random
import csv
import glob
import ast

import os
import sys
import time
import multiprocessing
try:
    # Python 2
    import Queue
except ImportError:
    # Python 2 and 3
    import queue as Queue

from inspyred import ec

# import helper functions
import filehandler as fh
from filehandler import WriteLogMsg
import config as cfg
from maphandler import generate_genom
from maphandler import generate_parameter
from maphandler import individual_filter
from maphandler import transform_individual_ascii_map
from maphandler import get_from_maphandler
from requirements import close_window
from __init__ import options

wrkDir = os.path.abspath('.')

#------------------------------------------------------------------------------  
#   Configuration / global variables
#------------------------------------------------------------------------------
     
# set file path for R
file_path_R = cfg.modelConfig.file_path_R

# file with HRUs
file_HRU = cfg.mapConfig.file_HRU

# external model folders
# minimal 1 model, maximal 4 models
model1_folder = cfg.modelConfig.model1_folder
try:
    model2_folder = cfg.modelConfig.model2_folder
    model3_folder = cfg.modelConfig.model3_folder
    model4_folder = cfg.modelConfig.model4_folder
except AttributeError:
    pass

file_model1 = cfg.modelConfig.file_model1
try:
    file_model2 = cfg.modelConfig.file_model2
    file_model3 = cfg.modelConfig.file_model3
    file_model4 = cfg.modelConfig.file_model4
except AttributeError:
    pass

# optimization algorithm
opt_algorithm = cfg.modelConfig.opt_algorithm

# # maximum number of possible land use options
max_range = cfg.modelConfig.max_range

# number of current generation
global nmbr_generation
nmbr_generation = 0
global custom_individuals
custom_individuals =[]

# array for the start individual
start_individual = []

#------------------------------------------------------------------------------  
#   Models: second level process handling for multiprocessing
#------------------------------------------------------------------------------
def process_handling(queue):
    """Process handling for multiprocessing in the second level."""
    
    # get tasks of the holding stack and call start function to execute the model
    while True:  

        try:
            # wait in maximum 1 second if no element is found in the queue
            # relevant for Linux 
            argument = queue.get(True,1)
            number = str(argument[0]) + str(argument[1])
            model = argument[2]
            msg = "Start model %s" %model
            WriteLogMsg(msg,number)
            fh.run_model(model, file_path_R, cfg.modelConfig.file_path_python, cfg.modelConfig.RPy2_available, number)
        
        #except:
        # Python 2 and 3
        except Exception as e:
            if type(e) != Queue.Empty:
                WriteLogMsg("Error: %s, %s" %(str(type(e)),str(e)))
                #WriteLogMsg("Unexpected error: %s" %sys.exc_info()[0])
                break
            else:
                break

#------------------------------------------------------------------------------  
#   Execute external models
#------------------------------------------------------------------------------
def execute_models(ind_number):
    """Execute external models."""
    
    msg = "Run external models ..."
    WriteLogMsg(msg,ind_number)
    
    if (ind_number == 1):
        folder_name = 'models'
    else:
        folder_name = 'models_%s' %(ind_number-1)

    # Run models depending on the optimization algorithm                    
    if opt_algorithm == "GA":
        # Run model 1
        fh.run_model(os.path.join(wrkDir, folder_name, model1_folder,file_model1), file_path_R, cfg.modelConfig.file_path_python, cfg.modelConfig.RPy2_available,ind_number)
           
    elif opt_algorithm == "NSGA2":                                        
        # for NSGAII, in maximum 4 objectives are allowed 
        # run the models in parallel on separate CPU-cores

        # count models and add tasks for the multiprocessing processes to the queue
        work_queue = multiprocessing.Queue()
        number_models = 1
        argument = [ind_number,number_models,os.path.join(wrkDir, folder_name, model1_folder,file_model1)]
        work_queue.put(argument)
        try:  
            argument = [ind_number,number_models+1,os.path.join(wrkDir, folder_name, model2_folder,file_model2)]
            work_queue.put(argument)
            number_models += 1 
            argument = [ind_number,number_models+1,os.path.join(wrkDir, folder_name, model3_folder,file_model3)]
            work_queue.put(argument) 
            number_models += 1   
            argument = [ind_number,number_models+1,os.path.join(wrkDir, folder_name, model4_folder,file_model4)]
            work_queue.put(argument) 
            number_models += 1
        except:
            pass        

        # define maximum number of processes to run in parallel
        if options.nthreads == "max cpu cores":
            nthreads = multiprocessing.cpu_count()
        else:
            nthreads = int(options.nthreads)

        # holds the multiprocessing processes
        jobs = []
        i = 0

        # create the multiprocessing processes
        while i < min(number_models, nthreads):
            p = multiprocessing.Process(target=process_handling,args=(work_queue,))
            jobs.append(p)
            p.start()
            i += 1

        for j in jobs:
            # wait until all multiprocessing processes are finished
            j.join()
    else:
        msg = "The selected optimization algorithm is not implemented."
        WriteLogMsg(msg,ind_number)
            
    msg = "Done."
    WriteLogMsg(msg,ind_number)

#------------------------------------------------------------------------------  
#   Individuals: first level process handling for multiprocessing
#------------------------------------------------------------------------------
def genome_process_handling(queue_arg):
    """Process handling for multiprocessing in the first level."""
    
    queue = queue_arg[0]
    map_info = queue_arg[1]
    patchID_map_info = queue_arg[2]
    header_all_info = queue_arg[3]
    
    while True:  
        try:
            # get the individual of the holding stack and start the evaluation
            # wait in maximum 1 second if no element is found in the queue
            # relevant for Linux
            argument = queue.get(True,1)
            ind_number = argument[0]
            individual = argument[1]
            msg = "Evaluation of individual %d" %(ind_number)   
            WriteLogMsg(msg,ind_number) 
            # change individual in genom.csv file
            fh.change_parameter_values(individual,ind_number)
            if (file_HRU == 'None' and cfg.modelConfig.map == 'True') or (file_HRU != 'None' and cfg.modelConfig.map == 'True' and cfg.mapConfig.file_ID_map != 'None'):
                # save or change individual as map in map.asc file
                transform_individual_ascii_map(individual,True,ind_number, map_info, patchID_map_info, header_all_info)
            # start external model         
            execute_models(ind_number)
        #except:
        # Python 3 and 2
        except Exception as e:
            if type(e) != Queue.Empty:
                WriteLogMsg("Error: %s, %s" %(str(type(e)),str(e)))
                #WriteLogMsg("Unexpected error: %s" %sys.exc_info()[0])
                break
            else:
                break

#------------------------------------------------------------------------------  
#   Evaluate individuals
#------------------------------------------------------------------------------
def evaluate(candidates, args):
    """Evaluate individuals."""   

    individuals = candidates

    # array for not accepted (infeasible) individuals
    not_accepted_ind = []

    # increment the generation number
    global nmbr_generation
    nmbr_generation += 1

    if len(individuals[0]) < 101:
        msg = "Population for generation %d: " % nmbr_generation
        WriteLogMsg(msg)
    
    i = 1    
    genome_queue = multiprocessing.Queue()

    # list with infeasible individuals if constrained_tournament_selection is selected
    if 'constrained_tournament_selection' in cfg.ea.selector:
        infeasible_ind = []

    # log the new population set
    for param in individuals:        
        if len(param) < 101:
            msg = "%d, %r" % (i, param)
            WriteLogMsg(msg)        
        
        # check if individuals are subject to special_termination
        # (genome consists of zeros)
        if all(item == 0 for item in param):
            not_accepted_ind.append(i)

        # check if individuals are feasible and constrained_tournament_selection is not selected
        # or constrained_tournament_selection is selected -> run models for all individuals
        elif (('constrained_tournament_selection' not in cfg.ea.selector) and individual_filter(param) == True) or ('constrained_tournament_selection' in cfg.ea.selector):
            # add tasks for the multiprocessing processes to the queue
            argument = [i,param]
            genome_queue.put(argument)

            # mark infeasible individuals for constrained_tournament_selection
            if 'constrained_tournament_selection' in cfg.ea.selector and individual_filter(param) == False:
                infeasible_ind.append(i)

        else:
            not_accepted_ind.append(i) 
        i += 1

    # transfer also the variables for map creation to the subprocesses
    map_info, patchID_map_info, header_all_info = get_from_maphandler() 
    queue_arg=[genome_queue, map_info, patchID_map_info, header_all_info]

    # check/create helping models folder for multiprocessing
    fh.copy_models(i-1)

    # a list with results for each individual
    fitness = []

    # count models
    number_models = 1
    try: 
        file_model2
        number_models += 1
        file_model3
        number_models += 1
        file_model4
        number_models += 1
    except:
        pass

    # define maximum number of processes to run in parallel
    if options.nthreads == "max cpu cores":
        nthreads = multiprocessing.cpu_count()
    else:
        nthreads = int(options.nthreads)

    # hold the multiprocessing processes
    jobs = []
    k = 0

    # every multiprocess for an individual generates later maximum number_models multiprocesses
    # if you have 2 cores and 4 models than you should generate only one multiprocess 
    # in the first level because you need the 2 cores for the multiprocessing of the models
    processes = 1
    while (max(nthreads, number_models) >= (processes * number_models)) and processes < i:
        processes += 1

    # create the multiprocessing processes
    while k < (processes-1):
        p = multiprocessing.Process(target=genome_process_handling,args=(queue_arg,))
        jobs.append(p)
        p.start() 
        k += 1

    for j in jobs:
        # wait until all multiprocessing processes are finished
        j.join()

    # Collect the fitness values of all individuals from one generation and return a list of them
    output_files = []
    external_models = []

    try:
        output_files.append(cfg.modelConfig.file_output1)
        output_files.append(cfg.modelConfig.file_output2)
        output_files.append(cfg.modelConfig.file_output3)
        output_files.append(cfg.modelConfig.file_output4)

    except AttributeError:
        pass

    try:
        external_models.append(cfg.modelConfig.model1_folder)
        external_models.append(cfg.modelConfig.model2_folder)
        external_models.append(cfg.modelConfig.model3_folder)
        external_models.append(cfg.modelConfig.model4_folder)

    except AttributeError:
        pass

    # add the logging informations from the child processes in the optimization_log file
    fh.join_ind_number_log()

    # add the model outputs of one generation to the special output file     
    fh.summarize_console_outputs(i-1,nmbr_generation,individuals, external_models, not_accepted_ind)
    
    # collect the fitness values of all individuals and models
    fitness = fh.collect_fitness_values(opt_algorithm, i-1, fitness, external_models, output_files, not_accepted_ind, cfg.mapConfig.file_worst_fitness)

    # for constrained_tournament_selection: print numbers of infeasible individuals  
    if 'constrained_tournament_selection' in cfg.ea.selector:
        WriteLogMsg("infeasible_ind: %s" % infeasible_ind)

    msg = "Fitness values are: %r \n" % fitness
    WriteLogMsg(msg)

    return fitness

#------------------------------------------------------------------------------
#   Genetic Algorithm function
#------------------------------------------------------------------------------
def GA():
    """Starts the optimization with the GA algorithm."""   

    begin = time.time()
    # initialize random generator with system time
    rand = random.Random()
    rand.seed()

    # original start individual of the input data
    global start_individual

    # generate the original start individual from the input data
    # return it including the non static land use indices
    start_individual, nonstatic_elements = generate_genom(max_range, file_HRU,cfg.mapConfig.file_ASCII_map, 
                                      cfg.mapConfig.file_transformation, cfg.mapConfig.file_ID_map, 
                                      cfg.mapConfig.four_neighbours)

    if len(start_individual) == 0:
        msg = "Error: The generated start individual has no elements."
        WriteLogMsg(msg) 
        raise SystemError("Error: The generated start individual has no elements.")
        close_window

    # determine that 'Bounder' conditions of candidates are equal to  
    # the integer values of the non static land use indices
    bounder_discrete = nonstatic_elements

    # initialize inspyred log files  
    stats_file,individ_file = fh.init_inspyred_logfiles()  

    # initialize and run GA
    ea = ec.GA(rand)
    # public attributes
    # GA is predefined with rank_selection
    if cfg.ea.selector != 'rank_selection':
        exec ("%s%s" % ('ea.selector = ', fh.preparing_attribute('selector',cfg.ea.selector)))
        msg = 'Selector of the optimization algorithm changed to: %s' % cfg.ea.selector
        WriteLogMsg(msg)
    # GA is predefined with generational_replacement
    if cfg.ea.replacer != 'generational_replacement':
        exec ("%s%s" % ('ea.replacer = ', fh.preparing_attribute('replacer',cfg.ea.replacer)))
        msg = 'Replacer of the optimization algorithm changed to: %s' % cfg.ea.replacer
        WriteLogMsg(msg)
    # specify how the new candidates should be varied
    # GA is predefined with n_point_crossover,bit_flip_mutation as variators
    if cfg.ea.variator != 'n_point_crossover,bit_flip_mutation' and cfg.ea.variator != 'bit_flip_mutation,n_point_crossover':
        exec ("%s%s" % ('ea.variator = ', fh.preparing_attribute('variator',cfg.ea.variator)))         
        msg = 'Variator of the optimization algorithm changed to: %s' % cfg.ea.variator
        WriteLogMsg(msg)       
    # GA is predefined with num_selected = pop_size
    if cfg.ea.num_selected != cfg.ea.pop_size:
        msg = 'Num_selected of the optimization algorithm changed to: %s' % cfg.ea.num_selected
        WriteLogMsg(msg) 
    exec ("%s%s" % ('ea.migrator = ', fh.preparing_attribute('migrator',cfg.ea.migrator)))
    exec ("%s%s" % ('ea.archiver = ', fh.preparing_attribute('archiver',cfg.ea.archiver)))
    if cfg.ea.archiver != 'best_archiver':
        msg = 'Archiver of the optimization algorithm changed to: %s' % cfg.ea.archiver
        WriteLogMsg(msg)
    exec ("%s%s" % ('ea.observer = ', fh.preparing_attribute('observer',cfg.ea.observer)))         
    # specify when the optimization should terminate
    exec ("%s%s" % ('ea.terminator = ', fh.preparing_attribute('terminator',cfg.ea.terminator)))

    # run optimization, when finished final_pop holds the results
    final_pop = ea.evolve(generator = generate_parameter, 
                    # evaluate is the function to start external models
                    # return results for the optimization algorithm
                    evaluator = evaluate, 
                    # define population size
                    pop_size = cfg.ea.pop_size,
                    # maximize or minimize the problem
                    maximize = cfg.ea.maximize, 
                    # bound the parameters to an interval
                    # choose integer values between 1 and max_range in this case
                    bounder = ec.DiscreteBounder(bounder_discrete), 
                    # minimum population diversity allowed (when using diversity_termination default 0.001)
                    min_diversity = cfg.ea.min_diversity,
                    # maximum number of evaluations (default pop_size) 
                    max_evaluations = cfg.ea.max_evaluations,
                    # maximum number of generations  
                    max_generations = cfg.ea.max_generations, 
                    # number of elites to consider (default 0)                    
                    num_elites = cfg.ea.num_elites,
                    # number of individuals to be selected (default NSGA2 pop_size)
                    num_selected = cfg.ea.num_selected,
                    # tournament size (default NSGA2 2) 
                    tournament_size = cfg.ea.tournament_size,
                    # the rate at which crossover is performed (default 1.0)
                    crossover_rate = cfg.ea.crossover_rate,
                    # mutation rate
                    mutation_rate = cfg.ea.mutation_rate,
                    # number of crossover points used (default 1)
                    num_crossover_points = cfg.ea.num_crossover_points,  
                    # a positive integer representing the number of 
                    # closest solutions to consider as a “crowd” (default 2)
                    crowding_distance = cfg.ea.crowding_distance,  
                    # statistic file
                    statistics_file = stats_file,
                    # individuals file
                    individuals_file = individ_file)                     

    # read out the best individuals
    final_arc = ea.archive

    # for constrained_tournament_selection: 
    # create a copy of final_arc only with feasible individuals (for csv file with best feasible solutions)
    if 'constrained_tournament_selection' in cfg.ea.selector:
        final_arc_feasible = [] 

    end = time.time()
    WriteLogMsg("The optimization process needed %d seconds." %(end-begin))

    msg = 'Best Solutions: \n'
    WriteLogMsg(msg)
    
    # save the map as ascii file in output folder

    f_count=1
    for f in final_arc:
        # for constrained_tournament_selection: with information if individual is infeasible
        # and copy feasible solutions in final_arc_feasible
        if 'constrained_tournament_selection' in cfg.ea.selector:
            if individual_filter(f.candidate) == False:
                WriteLogMsg("(infeasible) %s" % f)
                # save the map as ascii file in output folder
                if file_HRU == 'None' or (file_HRU != 'None' and cfg.mapConfig.file_ID_map != 'None'):
                    transform_individual_ascii_map(f.candidate,False,f_count,None,None,None,False)
            else:
                WriteLogMsg("%s" % f)
                # save the map as ascii file in output folder
                if file_HRU == 'None' or (file_HRU != 'None' and cfg.mapConfig.file_ID_map != 'None'):
                    transform_individual_ascii_map(f.candidate,False,f_count)
                final_arc_feasible.append(f)
        else:
            WriteLogMsg("%s" % f)
            # save the map as ascii file in output folder
            if file_HRU == 'None' or (file_HRU != 'None' and cfg.mapConfig.file_ID_map != 'None'):
                transform_individual_ascii_map(f.candidate,False,f_count)
        f_count += 1

    if cfg.ea.maximize == 'True':
        if 'constrained_tournament_selection' in cfg.ea.selector and individual_filter(f.candidate) == False:
            WriteLogMsg("\nFinal infeasible individual: %s, [%f]" % (max(final_pop).candidate,max(final_pop).fitness))
        else:
            WriteLogMsg("\nFinal individual: %s, [%f]" % (max(final_pop).candidate,max(final_pop).fitness))
    else:
        if 'constrained_tournament_selection' in cfg.ea.selector and individual_filter(f.candidate) == False:
            WriteLogMsg("\nFinal infeasible individual: %s, [%f]" % (min(final_pop).candidate,min(final_pop).fitness))
        else:
            WriteLogMsg("\nFinal individual: %s, [%f]" % (min(final_pop).candidate,min(final_pop).fitness))

    # save the map as ascii file in output folder

    # log the best solutions in a csv file
    fh.save_best_solutions(final_arc,1) 

    # for constrained_tournament_selection: log the best feasible solutions in a csv file
    if 'constrained_tournament_selection' in cfg.ea.selector:
        fh.save_best_solutions(final_arc_feasible,1)

#------------------------------------------------------------------------------
#   Nondominated Sorting Genetic Algorithm (NSGA-II)
#------------------------------------------------------------------------------
def NSGA2():
    """Starts the optimization with the NSGA-II algorithm."""   

    begin = time.time()
    # initialize random generator with system time
    rand = random.Random()
    rand.seed()

    global start_individual  
    global nmbr_generation
    global custom_individuals
    
    archive_is_available = False
    previous_archive = []


    if cfg.ea.start_from_previous_gen == True:
        generation_files = glob.glob("output/*_individuals_file.csv")  # Adjust path if necessary
        if generation_files:
            ea = ec.emo.NSGA2(rand)
            from inspyred.ec import Individual
            from inspyred.ec.emo import Pareto
            # Load the previously saved Pareto archive if it exists
            pareto_archive_file = "output/pareto_archive.csv"
            if os.path.exists(pareto_archive_file):
                with open(pareto_archive_file, "r") as file:
                    reader = csv.reader(file)
                    next(reader)  # Skip header
                    ea.archive = []
                    #ea.archive = [ec.emo.Pareto(eval(row[1]), eval(row[2])) for row in reader]
                    for row in reader:
                        candidate = eval(row[1])  # Assuming it's stored in a serializable format
                        fitness_vals = eval(row[2])
                        ind = Individual(candidate=candidate)
                        ind.fitness = Pareto(fitness_vals)
                        ea.archive.append(ind)
                print(f"Resumed Pareto front with {len(ea.archive)} individuals.")
                archive_is_available = True
                previous_archive = ea.archive
            else:
                print("No previous Pareto front archive found. Starting fresh.")

            # Sort files lexicographically to find the latest generation file
            generation_files.sort()
            last_saved_generation_file = generation_files[-1]  # Get the most recent file
            with open(last_saved_generation_file, mode='r') as file:
                reader = csv.reader(file)
                header = next(reader)  # Skip the header
                data = list(reader)
                nested_list = []


                if data:
                    # Parse the last generation number
                    nmbr_generation = int(data[-1][0])  # Last row's generation number
                    if cfg.ea.max_generations <= nmbr_generation:
                        msg = (f"Configured maximum generations ({cfg.ea.max_generations}) "
                               f"is less than or equal to the last saved generation ({nmbr_generation}). "
                               f"No further generation is performed")
                        WriteLogMsg(msg)
                        return 

                    # Extract the populations from lst generation
                    last_five_rows = data[-cfg.ea.pop_size:]
                    for row in last_five_rows:
                        # Extract columns from the 6th index onward
                        extracted_row = row[6:]  # Index 6 corresponds to column 7
                        # Convert values to integers and remove brackets if present
                        #cleaned_row = [int(value.strip("[]")) for value in extracted_row]
                        cleaned_row = [int(value.strip("[] ").strip()) for value in extracted_row]
                        # Append the cleaned row to the nested list
                        nested_list.append(cleaned_row)
                        custom_individuals.append(cleaned_row)

                    start_individual= nested_list
                    # Generate nonstatic elements
                    nonstatic_elements = generate_genom(
                        max_range, file_HRU, cfg.mapConfig.file_ASCII_map,
                        cfg.mapConfig.file_transformation, cfg.mapConfig.file_ID_map,
                        cfg.mapConfig.four_neighbours, return_only_nonstatic=True)
                else:
                    raise ValueError("The most recent file exists but contains no data.")

            msg = f"Resuming from generation {nmbr_generation}: {last_saved_generation_file}"
            WriteLogMsg(msg)
    else:
        # generate the original start individual from the input data
        # return it including the non static land use indices
        start_individual, nonstatic_elements = generate_genom(max_range, file_HRU,cfg.mapConfig.file_ASCII_map, 
                                      cfg.mapConfig.file_transformation, cfg.mapConfig.file_ID_map, 
                                      cfg.mapConfig.four_neighbours, return_only_nonstatic=False)
        print(f"Starting  with the individuals: {start_individual}")

        ea = ec.emo.NSGA2(rand)

    if len(start_individual) == 0:
        msg = "Error: The generated start individual has no elements."
        WriteLogMsg(msg) 
        raise SystemError("Error: The generated start individual has no elements.")
        close_window

    # determine that 'Bounder' conditions of candidates are equal to 
    # the integer values of the non static land use indices
    bounder_discrete = nonstatic_elements

    # initialize inspyred log files
    stats_file,individ_file = fh.init_inspyred_logfiles()

    # initialize and run NSGA2
    #ea = ec.emo.NSGA2(rand)
    # public attributes
    # NSGA2 is predefined with tournament_selection
    if cfg.ea.selector != 'tournament_selection':
        exec ("%s%s" % ('ea.selector = ', fh.preparing_attribute('selector',cfg.ea.selector)))
        msg = 'Selector of the optimization algorithm changed to: %s' % cfg.ea.selector
        WriteLogMsg(msg)
    # NSGA2 is predefined with nsga_replacement
    if cfg.ea.replacer != 'nsga_replacement':
        exec ("%s%s" % ('ea.replacer = ', fh.preparing_attribute('replacer',cfg.ea.replacer)))
        msg = 'Replacer of the optimization algorithm changed to: %s' % cfg.ea.replacer
        WriteLogMsg(msg)
    # NSGA2 is predefined with best_archiver
    if cfg.ea.archiver != 'best_archiver':
        exec ("%s%s" % ('ea.archiver = ', fh.preparing_attribute('archiver',cfg.ea.archiver)))
        msg = 'Archiver of the optimization algorithm changed to: %s' % cfg.ea.archiver
        WriteLogMsg(msg)
    exec ("%s%s" % ('ea.migrator = ', fh.preparing_attribute('migrator',cfg.ea.migrator)))
    # file observer prints after each generation the best, worst, mean etc. values into the statistic and individual file
    exec ("%s%s" % ('ea.observer = ', fh.preparing_attribute('observer',cfg.ea.observer)))      
    # specify how the new candidates should be varied
    exec ("%s%s" % ('ea.variator = ', fh.preparing_attribute('variator',cfg.ea.variator)))          
    # specify when the optimization should terminate
    exec ("%s%s" % ('ea.terminator = ', fh.preparing_attribute('terminator',cfg.ea.terminator)))   
    # NSGA2 is predefined with num_selected = pop_size
    if cfg.ea.num_selected != cfg.ea.pop_size:
        msg = 'Num_selected of the optimization algorithm changed to: %s' % cfg.ea.num_selected
        WriteLogMsg(msg) 
    # NSGA2 is predefined with tournament_size = 2
    if cfg.ea.tournament_size != 2:
        msg = 'Tournament_size of the optimization algorithm changed to: %s' % cfg.ea.tournament_size
        WriteLogMsg(msg)
  
    from inspyred.ec import observers  # Import existing observers
    from inspyred.ec.observers import pareto_archive_observer  # Import the new observer
    # Ensure all necessary observers are included
    ea.observer = [observers.file_observer, observers.best_observer, pareto_archive_observer]

    if not ea.archive:
        print("⚠ Warning: Pareto archive is empty before evolution starts.")
    else:
        print(f"✅ Archive loaded successfully with {len(ea.archive)} individuals.")
    
    #remaining_generations = cfg.ea.max_generations - nmbr_generation
    # run optimization, when finished final_pop holds the results
    final_pop = ea.evolve(generator = generate_parameter, 
                    # evaluate is the function to start external models
                    # return results for the optimization algorithm
                    evaluator = evaluate, 
                    # define population size
                    pop_size = cfg.ea.pop_size, 
                    # maximize or Minimize the problem (default True)                  
                    maximize = cfg.ea.maximize,   
                    # bound the parameters to an interval 
                    # DiscreteBounder: choose integer values between 1 and max_range
                    bounder = ec.DiscreteBounder(bounder_discrete),
                    # minimum population diversity allowed (when using diversity_termination default 0.001)
                    min_diversity = cfg.ea.min_diversity,
                    # maximum number of generations  
                    max_generations = cfg.ea.max_generations,
                    # maximum number of evaluations (default pop_size) 
                    max_evaluations = cfg.ea.max_evaluations,
                    # number of elites to consider (default 0)                    
                    num_elites = cfg.ea.num_elites,
                    # number of individuals to be selected (default NSGA2 pop_size)
                    num_selected = cfg.ea.num_selected,
                    # tournament size (default NSGA2 2) 
                    tournament_size = cfg.ea.tournament_size,
                    # rate at which crossover is performed (default 1.0)
                    crossover_rate = cfg.ea.crossover_rate,
                    # rate at which mutation is performed (default 0.1)
                    mutation_rate = cfg.ea.mutation_rate,
                    # number of crossover points used (default 1)
                    num_crossover_points = cfg.ea.num_crossover_points,  
                    # a positive integer representing the number of 
                    # closest solutions to consider as a “crowd” (default 2)
                    crowding_distance = cfg.ea.crowding_distance, 
                    # statistic file
                    statistics_file = stats_file,
                    # individuals file
                    individuals_file = individ_file, 
                    is_available = archive_is_available, 
                    archive = previous_archive, 
                    custom_individual= custom_individuals, num_generation = nmbr_generation)                     

    final_arc = ea.archive

    # for constrained_tournament_selection: 
    # create a copy of final_arc only with feasible individuals (for csv file with best feasible solutions)
    if 'constrained_tournament_selection' in cfg.ea.selector:
        final_arc_feasible = []    

    end = time.time()
    msg = "The optimization process needed %d seconds." %(end-begin)
    fh.WriteLogMsg(msg)

    msg = 'Best Solutions: \n'
    WriteLogMsg(msg)

    f_count=1
    for f in final_arc:
        # for constrained_tournament_selection: with information if individual is infeasible
        # and copy feasible solutions in final_arc_feasible
        if 'constrained_tournament_selection' in cfg.ea.selector:
            if individual_filter(f.candidate) == False:
                WriteLogMsg("(infeasible) %s" % f)
                # save the map as ascii file in output folder
                if file_HRU == 'None' or (file_HRU != 'None' and cfg.mapConfig.file_ID_map != 'None'):
                    transform_individual_ascii_map(f.candidate,False,f_count,None,None,None,False)
            else:
                WriteLogMsg("%s" % f)
                # save the map as ascii file in output folder
                if file_HRU == 'None' or (file_HRU != 'None' and cfg.mapConfig.file_ID_map != 'None'):
                    transform_individual_ascii_map(f.candidate,False,f_count)
                final_arc_feasible.append(f)

        else:
            msg = "%s" % f
            #msg = "%f" % f
            WriteLogMsg(msg)
            # save the map as ascii file in output folder
            if file_HRU == 'None' or (file_HRU != 'None' and cfg.mapConfig.file_ID_map != 'None'):
                transform_individual_ascii_map(f.candidate,False,f_count)

        if f_count == 1:
            len_fitness = len(f.fitness)
        
        f_count += 1

    # plot the best solution in a 2, 3 or 4 dimensional plot 
    # 2 dimensional plot
    if (cfg.ea.plot_results == True and len_fitness == 2):
        # log the best solutions in a csv file
        fh.save_best_solutions(final_arc,2)  
        import pylab
        x = []
        y = []
        for f in final_arc:
            x.append(f.fitness[0])
            y.append(f.fitness[1])
        pylab.scatter(x, y, color='r')
        fh.savePlot_png(opt_algorithm)
        pylab.show()

        # for constrained_tournament_selection: create a second plot with feasible solutions
        if 'constrained_tournament_selection' in cfg.ea.selector:
            # log the best feasible solutions in a csv file
            fh.save_best_solutions(final_arc_feasible,2)  
            x = []
            y = []
            for f in final_arc_feasible:
                x.append(f.fitness[0])
                y.append(f.fitness[1])
            pylab.scatter(x, y, color='r')
            fh.savePlot_png(opt_algorithm)
            pylab.show()


    # 3 and 4 dimensional plots
    if (cfg.ea.plot_results == True and (len_fitness == 3 or len_fitness == 4)):
        import warnings 

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt 

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            x = []
            y = []
            if len_fitness == 3 or len_fitness == 4:
                z = []
            if len_fitness == 4:
                c = []
            for f in final_arc:
                x.append(f.fitness[0])
                y.append(f.fitness[1])
                if len(f.fitness) == 3 or len(f.fitness) == 4:
                    z.append(f.fitness[2])
                if len(f.fitness) == 4:
                    c.append(f.fitness[3])

            if len_fitness == 3:
                # log the best solutions in a csv file
                fh.save_best_solutions(final_arc,3)
                ax.scatter(x, y, z, c='r')
            if len_fitness == 4:
                # log the best solutions in a csv file
                fh.save_best_solutions(final_arc,4)
                ax.scatter(x, y, z, c=c, cmap=plt.hot())
            fh.savePlot_png(opt_algorithm)
            plt.show() 

            # for constrained_tournament_selection: create a second plot with feasible solutions
            if 'constrained_tournament_selection' in cfg.ea.selector:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')

                x = []
                y = []
                if len_fitness == 3 or len_fitness == 4:
                    z = []
                if len_fitness == 4:
                    c = []
                for f in final_arc_feasible:
                    x.append(f.fitness[0])
                    y.append(f.fitness[1])
                    if len(f.fitness) == 3 or len(f.fitness) == 4:
                        z.append(f.fitness[2])
                    if len(f.fitness) == 4:
                        c.append(f.fitness[3])

                if len_fitness == 3:
                    # log the best solutions in a csv file
                    fh.save_best_solutions(final_arc_feasible,3)
                    ax.scatter(x, y, z, c='r')
                if len_fitness == 4:
                    # log the best solutions in a csv file
                    fh.save_best_solutions(final_arc_feasible,4)
                    ax.scatter(x, y, z, c=c, cmap=plt.hot())
                fh.savePlot_png(opt_algorithm)
                plt.show()
    
    # print results without plotting (if plot_results was set to false)
    if cfg.ea.plot_results == False:
        if len_fitness == 2:
            # log the best solutions in a csv file
            fh.save_best_solutions(final_arc,2)
            if 'constrained_tournament_selection' in cfg.ea.selector:
                # log the best feasible solutions in a csv file
                fh.save_best_solutions(final_arc_feasible,2)
        if len_fitness == 3:
            # log the best solutions in a csv file
            fh.save_best_solutions(final_arc,3)
            if 'constrained_tournament_selection' in cfg.ea.selector:
                # log the best feasible solutions in a csv file
                fh.save_best_solutions(final_arc_feasible,3)
        if len_fitness == 4:
            # log the best solutions in a csv file
            fh.save_best_solutions(final_arc,4)
            if 'constrained_tournament_selection' in cfg.ea.selector:
                # log the best feasible solutions in a csv file
                fh.save_best_solutions(final_arc_feasible,4)

#------------------------------------------------------------------------------  
#   Test functions
#------------------------------------------------------------------------------   
if __name__ == "__main__":    

    start_individual, nonstatic_elements = generate_genom(max_range, file_HRU,cfg.mapConfig.file_ASCII_map, 
                                  cfg.mapConfig.file_transformation, cfg.mapConfig.file_ID_map, 
                                  cfg.mapConfig.four_neighbours, return_only_nonstatic = False)
    """print("pop Laenge: %s" %len(pop))
    print("pop type: %s" %type(pop))
    print("Listenelemente:")
    for i in pop:
        print i
    print("elemente mit index:")
    for i in range(0, len(pop)):
        print pop[i]
    print("start_individual in NSGA2: "%pop)"""

  
#------------------------------------------------------------------------------
#
#   EOF
#
#------------------------------------------------------------------------------
