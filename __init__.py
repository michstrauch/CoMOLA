# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#
#   Name:       __init__.py
#   Purpose:    This script includes the main function of the optimization tool.  
#               It uses the inspyred package (http://inspyred.github.io/).
#               
#   Author:     Carola Paetzold, Michael Strauch
#   Contact:    michael.strauch@ufz.de
#
#               Helmholtz Centre for Environmental Research - UFZ
#               Department Computational Landscape Ecology - CLE
#               Permoserstrasse 15
#               D-04318 Leipzig, Germany
#               http://www.ufz.de
#
#   Created:    Mar 19 2014
#
#   Copyright:  (c) Carola Paetzold, Michael Strauch 2018
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
import config as cfg
import requirements as req
import time
import shutil
import sys
import distutils.core
import argparse
from argparse import ArgumentParser

# ------------------------------------------------------------------------------
# Maximum number of threads as command line argument
# ------------------------------------------------------------------------------

# defaults
default_nthreads = "max cpu cores"

# parsing arguments

parser = ArgumentParser(
    description = 'CoMOLA - Constrained Multi-objective Optimization of Land use Allocation'
)

parser.add_argument(
    "-t", "--threads",
    help = "number of threads to use, defaults to " + str(default_nthreads),
    dest = "nthreads",
    default = default_nthreads,
    metavar = 4
)

options = parser.parse_args()

# Get absolute path of this file
wrkDir = os.path.abspath(".")

#-------------------------------------------------------------------------------------  
#   Update/delete the helping model folders
#-------------------------------------------------------------------------------------
def update_help_folders(fh):
    """Delete the helping model folders if cfg.modelConfig.del_help_folders = 'True'
       else update files/folders from cfg.modelConfig.update_filesx 
       if this variable is not 'None' (default).

    """
    
    if cfg.modelConfig.del_help_folders == 'True':
        # delete helping models folder
        fh.delete_models()
    else:
        # update files in helping model folders
        fh.WriteLogMsg("Update files in helping folders ...")  
        i = 1
        folder_name = 'models_%s' %i
        while os.path.exists(os.path.join(wrkDir,folder_name)):
            for number in range(1,5):
                try:
                    # if update_files is activated, check if files have changed; if so, update these files in the helping folders
                    if number == 1:
                        sub_folder =  cfg.modelConfig.model1_folder
                        update_files = cfg.modelConfig.update_files1
                    elif number == 2:
                        sub_folder =  cfg.modelConfig.model2_folder
                        update_files = cfg.modelConfig.update_files2
                    elif number == 3:
                        sub_folder =  cfg.modelConfig.model3_folder
                        update_files = cfg.modelConfig.update_files3
                    else:
                        sub_folder =  cfg.modelConfig.model4_folder
                        update_files = cfg.modelConfig.update_files4
                             
                    if update_files != 'None':
                        for file in update_files.split(','): 
                            if os.path.isfile(os.path.join(wrkDir, folder_name, sub_folder, file)):
                                if os.path.getmtime(os.path.join(wrkDir, folder_name, sub_folder, file)) \
                                != os.path.getmtime(os.path.join(wrkDir, 'models', sub_folder, file)): 
                                    shutil.copy(os.path.join(wrkDir, 'models', sub_folder, file), \
                                                os.path.join(wrkDir, folder_name, sub_folder, file))
                            elif os.path.isdir(os.path.join(wrkDir, folder_name, sub_folder, file)):
                                if os.path.getmtime(os.path.join(wrkDir, folder_name, sub_folder, file)) \
                                != os.path.getmtime(os.path.join(wrkDir, 'models', sub_folder, file)): 
                                    distutils.dir_util.copy_tree(os.path.join(wrkDir, 'models', sub_folder, file), \
                                                os.path.join(wrkDir, folder_name, sub_folder, file))
                            # os.path.isfile(os.path.join(wrkDir, folder_name, sub_folder, file)) does not exist
                            else:
                                if os.path.isdir(os.path.join(wrkDir, 'models', sub_folder, file)):
                                    shutil.copytree(os.path.join(wrkDir, 'models', sub_folder, file), \
                                                os.path.join(wrkDir, folder_name, sub_folder, file))
                                elif os.path.isfile(os.path.join(wrkDir, 'models', sub_folder, file)):
                                    shutil.copy(os.path.join(wrkDir, 'models', sub_folder, file), \
                                                os.path.join(wrkDir, folder_name, sub_folder, file))
                                else:
                                    fh.WriteLogMsg("Path %s can`t be updated, because the original path does not exist. Please check config.ini." %os.path.join(wrkDir, 'models', sub_folder, file)) 
                                    raise SystemError("Path %s can`t be updated, because the original path does not exist. Please check config.ini." %os.path.join(wrkDir, 'models', sub_folder, file))
                                    req.close_window  
                except:
                    #print "Unexpected error:", sys.exc_info()[0]
                    #raise  
                    break    
            i += 1
            folder_name = 'models_%s' %i
        fh.WriteLogMsg("Update files is done.")  
    
#------------------------------------------------------------------------------  
#   Entry point - the main function
#------------------------------------------------------------------------------
def main():
    """Starts the optimization process."""
    
    # get start time 
    begin=time.time()
    
    print("number of threads: {}".format(options.nthreads))
    
    # if a help file of a previous run exists, delete it 
    if os.path.exists(os.path.join(wrkDir, "output", "help_file.txt")):
        os.remove(os.path.join(wrkDir, "output", "help_file.txt"))
    
    # check the Python version and requirements of the tool
    req.check_requirements()
    
    # if no exit activated until now, execute the next steps     
    import optiAlgorithm
    import filehandler as fh
    
    # create an output folder for the log-files
    fh.create_output_folder()
    
    # initialize log file
    fh.InitLogFile()
    
    # save the input data
    fh.save_input_data()
    
    # inspyred logging for troubleshooting
    fh.inspyred_logging()
    
    # create a help file
    fh.save_timestamp(fh.timestamp_file)
    
    # check if helping folders should be created or if files should be updated
    update_help_folders(fh)
    
    # if repair_mutation is selected then special_termination should be selected too
    if 'repair_mutation' in cfg.ea.variator:
        if 'special_termination' not in cfg.ea.terminator:
            termination = cfg.ea.terminator
            termination = '%s,special_termination' % termination
            cfg.ea.terminator = termination
            fh.WriteLogMsg("Special_termination is added because of the selection of the repair_mutation.")
    
    # check if file is available with worst fitness values
    if cfg.mapConfig.file_worst_fitness == 'None' and (cfg.mapConfig.file_transformation != 'None' or cfg.mapConfig.file_difference!= 'None'):
        fh.WriteLogMsg("No input file available with worst fitness values. Please check config.ini.") 
        raise SystemError("Error: No input file available with worst fitness values.")
        req.close_window
       
    # run optimization
    algorithm_selected = cfg.modelConfig.opt_algorithm
    msg = "Selected algorithm: %s" % algorithm_selected 
    fh.WriteLogMsg(msg)    
    
    if algorithm_selected == "GA":    
        optiAlgorithm.GA() 
        fh.plot_statistics_file()    
    elif algorithm_selected  == "NSGA2":
        optiAlgorithm.NSGA2()
        
    # delete the help file
    os.remove(os.path.join(wrkDir, "output", "help_file.txt"))
    
    # end time 
    end=time.time()
    msg = "Total runtime of the tool: %d seconds." %(end-begin)
    fh.WriteLogMsg(msg)

#------------------------------------------------------------------------------  
#   Call main function
#------------------------------------------------------------------------------   
if __name__ == "__main__":        
    main()
    
###############################################################################
#
#   EOF
#
###############################################################################
