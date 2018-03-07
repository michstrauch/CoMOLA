# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#
#   Name:       config.py
#   Purpose:    This script includes all configuration parameters
#               for an optimization run
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
#   Copyright:  (c) Carola Paetzold / Christian Schweitzer / Michael Strauch 2017
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
from distutils.util import strtobool
try:
    import ConfigParser
except ImportError:
    # for Python 3
    import configparser as ConfigParser
import requirements as req

# get absolute path of this file
wrkDir = os.path.abspath(".")

#-------------------------------------------------------------------------------------  
#   Read ini file
#-------------------------------------------------------------------------------------
def read_ini_file():
    """Read the ini file."""
    
    # RawConfigParser(Default): default value dictionary are set for all sections      
    config = ConfigParser.RawConfigParser()
    # dict_default holds the default values of the configuration parameters
    dict_default_model = {}
    # default model settings
    # custom settings under [config_model] in ini file
    dict_default_model.update({'file_path_r' : "C:\\Program Files\\R\\R-3.0.1\\bin\\x64\\R.exe"})
    dict_default_model.update({'file_path_python' : "C:\\Python27\\Python.exe"})
    dict_default_model.update({'opt_algorithm' : 'NSGA2'})
    dict_default_model.update({'rpy2_available' : 'False'})
    dict_default_model.update({'map' : 'False'})
    dict_default_model.update({'del_help_folders' : 'False'})
    dict_default_model.update({'update_files1' : 'None'})
    dict_default_model.update({'update_files2' : 'None'})
    dict_default_model.update({'update_files3' : 'None'})
    dict_default_model.update({'update_files4' : 'None'})
    dict_default_alg = {}
    # default optimization settings
    # custom settings under [config_optimization_algorithm] in ini file
    dict_default_alg.update({'pop_size' : '100'})
    dict_default_alg.update({'max_generations' : '1'})
    dict_default_alg.update({'max_repair_trials' : '0'})
    dict_default_alg.update({'write_tabu_memory' : 'False'})
    dict_default_alg.update({'plot_results' : 'False'})
    dict_default_alg.update({'maximize' : 'True'})
    dict_default_alg.update({'selector' : 'default_selection'})
    dict_default_alg.update({'variator' : 'default_variation'})
    dict_default_alg.update({'replacer' : 'default_replacement'})
    dict_default_alg.update({'migrator' : 'default_migration'})
    dict_default_alg.update({'archiver' : 'default_archiver'})
    dict_default_alg.update({'observer' : 'file_observer'})
    dict_default_alg.update({'terminator' : 'default_termination'})
    dict_default_alg.update({'crossover_rate' : '1.0'})
    dict_default_alg.update({'priority' : 'True'})
    dict_default_alg.update({'feasible_first_pop' : 'False'})
    dict_default_alg.update({'extreme_seeds' : 'False'})
    dict_default_alg.update({'num_crossover_points' : '1'})
    dict_default_alg.update({'mutation_rate' : '0.1'})
    dict_default_alg.update({'num_elites' : '0'})
    dict_default_alg.update({'min_diversity' : '0.001'})
    dict_default_alg.update({'tournament_size' : '2'})
    dict_default_alg.update({'num_selected' : '1'})
    dict_default_alg.update({'crowding_distance' : '2'})
    dict_default_alg.update({'penalty_function' : '0'})
    # default map analysis settings
    # custom settings under [config_map_analysis] in ini file
    dict_default_map = {}
    dict_default_map.update({'file_patch_map' : 'None'})
    dict_default_map.update({'file_landuse_map' : 'None'})
    dict_default_map.update({'file_hru' : 'None'})
    dict_default_map.update({'file_transition' : 'None'})
    dict_default_map.update({'file_area' : 'None'})
    dict_default_map.update({'four_neighbours' : 'False'})
    dict_default_map.update({'file_worst_fitness' : 'None'})
    # dict_model, dict_alg and dict_map hold all configuration parameters and their values as dictionaries 
    dict_model = {}
    dict_alg = {}
    dict_map = {}
    
    try:
        # read ini and transfer section list into dictionaries
        config.read(os.path.join(wrkDir,"config.ini"))
        config_model = config.items('config_model')
        config_opti_alg = config.items('config_optimization_algorithm')
        config_map = config.items('config_map_analysis')
    
        for element in config_model:
            dict_model[element[0]] = element[1]
        for element in config_opti_alg:
            dict_alg[element[0]] = element[1]
        for element in config_map:
            dict_map[element[0]] = element[1]
            
        # default exceptions
        try:
            dict_default_alg.update({'max_evaluations' : dict_alg['pop_size']})
        except KeyError: # if no pop_size is given in ini file
            dict_alg.update({'pop_size' : dict_default_alg['pop_size']}) 
            dict_default_alg.update({'max_evaluations' : dict_alg['pop_size']})  
        if dict_model['opt_algorithm'] == 'NSGA2':
            dict_default_alg['selector'] = 'tournament_selection'
            dict_default_alg['replacer'] = 'nsga_replacement'
            dict_default_alg['archiver'] = 'best_archiver'
            dict_default_alg['num_selected'] = dict_alg['pop_size']
        elif dict_model['opt_algorithm'] == 'GA':
            dict_default_alg['selector'] = 'rank_selection'
            dict_default_alg['replacer'] = 'generational_replacement'
            dict_default_alg['variator'] = 'n_point_crossover,bit_flip_mutation'
            dict_default_alg['num_selected'] = dict_alg['pop_size']
            dict_default_alg['archiver'] = 'best_archiver'
        try:
            if dict_alg['selector'] == 'truncation_selection':
                dict_default_alg['num_selected'] = dict_alg['pop_size']
        except KeyError: # troubleshooting if no value for selector is given in the ini file
            pass
            
        # check if default values should be added to the dictionaries
        for element in dict_default_model.keys():
            if element not in dict_model:
                dict_model.update({element : dict_default_model[element]})
        for element in dict_default_alg.keys():
            if element not in dict_alg:
                dict_alg.update({element : dict_default_alg[element]})
        for element in dict_default_map.keys():
            if element not in dict_map:
                dict_map.update({element : dict_default_map[element]})                
                
    except ConfigParser.NoSectionError:       
        print("An error occurred when the program tries to read the ini file. Exists a config.ini file in the main folder? Exist a config_model and a config_optimization_algorithm section in the config.ini?")
        req.close_window() 
       
    return dict_alg, dict_model, dict_map

# read the ini file
# ini_list holds dict_alg, dict_model and dict_map from read_ini_file
ini_list = read_ini_file()
dict_alg = ini_list[0]
dict_model = ini_list[1] 
dict_map = ini_list[2] 

#------------------------------------------------------------------------------
#   Config model
#------------------------------------------------------------------------------
class ModelConfig:
    
    def __init__(self):
        """Configure the model settings.
            
           Variables:                   Description:
           ----------                   ------------
           file_path_R                | file path for R  
           file_path_Python           | file path for Python          
           modelx_folder              | folder name of model x (1 <= x <= 4)
           file_modelx                | file name of the model x script
           file_outputx               | file name of the output file from model x
           file_outputx               | file names from model x which should be updated 
                                      | in the helping folders at the start of the tool 
           max_range                  | maximum number of possible land use options
           opt_algorithm (string)     | definition of the optimization algorithm,
                                        available choices are GA or NSGA2   
           RPy_available (string)     | if RPy2 is available then True, False otherwise 
           map                        | if True then transfer individuals as ascii maps into the model folders,
                                      | else save the individuals as string of integers in a csv file 
           del_help_folders           | if True than delete and create the helping folders at the start of the process,
                                      | else update only the changed files in the existing helping folders 
                                      
        """
        # set current working directory
        self.file_ini = wrkDir

        # set up file path for R (not necessary for RPy2)
        self.file_path_R = dict_model['file_path_r']
        
        # set up file path for Python
        self.file_path_python = dict_model['file_path_python']
        
        # should the helping folders be deleted and created for each optimization?
        self.del_help_folders = dict_model['del_help_folders']
        
        # folder of the first model - path will be generated dynamically  
        self.model1_folder = dict_model['model1_folder']
           
        # file of the first model script
        self.file_model1 = dict_model['file_model1']
 
        # file with resulting fitness values
        self.file_output1 = dict_model['file_output1']
        
        # file with the files which should be updated in the helping folders
        self.update_files1 = dict_model['update_files1']
        self.update_files2 = dict_model['update_files2']
        self.update_files3 = dict_model['update_files3']
        self.update_files4 = dict_model['update_files4']
        
        try:
            # folder of the second model 
            self.model2_folder = dict_model['model2_folder']
            
            # file of the second model script
            self.file_model2 = dict_model['file_model2']
            
            # file with resulting fitness values
            self.file_output2 = dict_model['file_output2']
            
            # folder of the third model
            self.model3_folder = dict_model['model3_folder']
            
            # file of the third model script
            self.file_model3 = dict_model['file_model3']
            
            # file with resulting fitness values
            self.file_output3 = dict_model['file_output3']
            
            # folder of the fourth model 
            self.model4_folder = dict_model['model4_folder']
            
            # file of the fourth model script
            self.file_model4 = dict_model['file_model4']
            
            # file with resulting fitness values
            self.file_output4 = dict_model['file_output4']
        except:
            pass
        
        # maximum number of possible land use options
        self.max_range = int(dict_model['max_range'])

        # optimization algorithm
        self.opt_algorithm = dict_model['opt_algorithm']
        
        # RPy2 is available
        self.RPy2_available = dict_model['rpy2_available']
        
        # individual transfer to model folder as string (False) or as ascii map (True)
        self.map = dict_model['map']

modelConfig = ModelConfig()

#------------------------------------------------------------------------------
#   Optimization algorithm configuration
#------------------------------------------------------------------------------
class EaConfig:
    """Parameter settings for the evolutionary algorithm"""

    def __init__(self):
        
        self.pop_size = int(dict_alg['pop_size'])
        self.maximize = dict_alg['maximize']        
        self.selector = dict_alg['selector']
        self.variator = dict_alg['variator']
        self.replacer = dict_alg['replacer']
        self.migrator = dict_alg['migrator']
        self.archiver = dict_alg['archiver']
        self.observer = dict_alg['observer']
        self.terminator = dict_alg['terminator']        
        self.crossover_rate = float(dict_alg['crossover_rate'])
        self.priority = dict_alg['priority']
        self.feasible_first_pop = dict_alg['feasible_first_pop']
        self.extreme_seeds = strtobool(dict_alg['extreme_seeds'])
        self.num_crossover_points = int(dict_alg['num_crossover_points'])        
        self.mutation_rate = float(dict_alg['mutation_rate'])
        self.num_elites = int(dict_alg['num_elites'])
        self.min_diversity = float(dict_alg['min_diversity'])
        self.num_selected = int(dict_alg['num_selected'])
        self.crowding_distance = int(dict_alg['crowding_distance'])
        self.tournament_size = int(dict_alg['tournament_size'])
        self.penalty_function = int(dict_alg['penalty_function'])
        self.max_generations = int(dict_alg['max_generations'])
        self.max_evaluations = int(dict_alg['max_evaluations'])
        self.max_repair_trials = int(dict_alg['max_repair_trials'])
        self.write_tabu_memory = strtobool(dict_alg['write_tabu_memory'])
        self.plot_results = strtobool(dict_alg['plot_results'])
        
ea = EaConfig()

#------------------------------------------------------------------------------
#   Map analysis configuration
#------------------------------------------------------------------------------
class MapConfig:
    """Parameter settings for the map analysis"""

    def __init__(self):
      
        # file with transition matrix defining which land use class can change into which other
        self.file_transformation = dict_map['file_transition']
        # file with total area (min-max) constraints
        self.file_difference = dict_map['file_area']
        # optional file with HRUs (start with HRUs as basis for the genom instead of a map)
        self.file_HRU = dict_map['file_hru']
        # file for the original ASCII map
        self.file_ASCII_map = dict_map['file_landuse_map']
        # file for the patch ID map
        self.file_ID_map = dict_map['file_patch_map']
        # should 4 or 8 neighbouring cells be considered for patch generation
        self.four_neighbours = dict_map['four_neighbours']
        # file for the default worst fitness values
        self.file_worst_fitness = dict_map['file_worst_fitness']
        
mapConfig = MapConfig()

#------------------------------------------------------------------------------
#
#   EOF
#
#------------------------------------------------------------------------------
