# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#
#   Name:       requirements.py
#   Purpose:    This module is responsible for checking requirements of the tool.
#
#   Author:     Carola  Paetzold
#   Contact:    carola.paetzold@ufz.de
#
#               Helmholtz Centre for Environmental Research - UFZ
#               Department Computational Landscape Ecology - CLE
#               Permoserstrasse 15
#               D-04318 Leipzig, Germany
#               http://www.ufz.de
#
#   Created:    Wed May 21 2014
#
#   Copyright:  (c) Carola  Paetzold 2014
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
import sys
import importlib
import config as cfg

#------------------------------------------------------------------------------  
#   Check the Python Version and requirements of the tool
#------------------------------------------------------------------------------
def check_requirements():
    # path to the used Python version
    python_folder = sys.executable
    accepted_versions = [(3,9),(3,10),(3,11)] 
    # check for accepted python versions
    if not sys.version_info[:2] in accepted_versions:
        print ("Warning: The tool is testet for Python 3.11. You started the tool with Python version from %s." %python_folder)  
    # check if the libraries matplotlib and RPy2 exist
    matplotlib_exist = importlib.util.find_spec("matplotlib")
    rpy2_exist = importlib.util.find_spec("rpy2")
    if not matplotlib_exist:
        if not rpy2_exist:
            print("Error, you need the matplotlib library but it is not installed in Python version from %s. If you would like to use RPy2 than you should install RPy2 too." %python_folder)
        else:
            print("Error, you need the matplotlib library but it is not installed in Python version from %s." %python_folder)
        close_window()     
    #  matplotlib exist, check if RPy2 is used and installed
    elif cfg.modelConfig.RPy2_available == "True" and not rpy2_exist:
        #if not rpy2_exist:
        print("Error, you want use RPy2 but it is not installed in Python version from %s." %python_folder)
        close_window() 
#------------------------------------------------------------------------------  
#   Close the program after a key entry
#------------------------------------------------------------------------------        
def close_window():
    """Prevent that the window is closed after the error message appears."""
    
    # wait for an input before close the window
    # the user can read the error messages
    # raw_input does not exist in Python 3
    #raw_input('Please press enter to close the program.')
    # Python 3
    input('Please press enter to close the program.')
    exit(1) 
  
#------------------------------------------------------------------------------
#
#   EOF
#
#------------------------------------------------------------------------------
