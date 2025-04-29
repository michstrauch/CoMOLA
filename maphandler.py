# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#
#	Name: l	  maphandler.py
#	Purpose:	This module includes all functions for the usage of maps.
#
#	Author:	 Carola Paetzold, Sven Lautenbach, Michael Strauch
#	Contact:	michael.strauch@ufz.de
#
#				Helmholtz Centre for Environmental Research - UFZ
#				Department Computational Landscape Ecology - CLE
#				Permoserstrasse 15
#				D-04318 Leipzig, Germany
#				http://www.ufz.de
#
#	Created:	We Oct 29 2014
#
#	Copyright:	(c) Carola Paetzold / Sven Lautenbach / Michael Strauch 2017
#
#	Licence:	This program is free software:
#				you can redistribute it and/or modify it under the terms
#				of the GNU General Public License as published by the
#				Free Software Foundation, either version 3 of the License,
#				or (at your option) any later version. This program is
#				distributed in the hope that it will be useful, but
#				WITHOUT ANY WARRANTY; without even the implied warranty
#				of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#				See the GNU General Public License for more details.
#				You should have received a copy of the GNU General
#				Public License along with this program.
#				If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------
import os
import sys
import csv
import time
import random
import numpy as np
from math import pow

import requirements as req
from filehandler import WriteLogMsg
from filehandler import WriteMap
from filehandler import WriteCandidateList

import config as cfg

wrkDir = os.path.abspath('.')
# first_ind is True for the first call of the generate_parameter function
first_ind = True
# array for the start individual
start_individual = []
# array for original ASCII map
map = []
# array for the patch ID map
patchID_map = []
# array for the header information from ASCII map
header = []
# list with header information from ASCII map - one string per line
header_all = []
# array for land use transition constraints
trans_matrix = []
# array for total area constraints of each land use class
min_max_diff = []
# array for static land use types
static_elements = []
# array for non static land use types
nonstatic_elements = []
# array for the area percentage of the patches
map_proportion = []
# dictionary for possible land use options per element
possible_elements = {}
# array with segments of impossible candidates
impossible_cand = []
# Boolean variable for starting the special_termination
end_optimization = False
# for constraint tournament selection
# dictionary as memory for the analysis results of the candidate violation
violation_memory = {}
# dictionary for area percentage of the static land use categories
static_area = {}
# dictionary with patch indices and total area per land use class
start_land_cover = {}
# dictionary with information per land use class if extreme seeds were created before
extreme_seeds_dict = {}
global custom_index
custom_index = 0

#-------------------------------------------------------------------------------------	
#	Start the termination of the optimization algorithm
#-------------------------------------------------------------------------------------
def logical_termination():
	"""If repair_mutation is selected and the last possible candidate is set as an individual 
		then the variable end = True and the termination can start after the determination of 
		the fitness values of the last population.

	"""
	# Boolean variable for starting the special_termination
	global end_optimization
	return end_optimization
	
#-------------------------------------------------------------------------------------	
#	Read an ASCII-map
#-------------------------------------------------------------------------------------
def read_ascii_map(file):
	"""Read an ascii file.
		Return the map as matrix and header data as array.
	
		input data: ascii file
	"""
	
	msg = "Read ascii file ..."
	WriteLogMsg(msg) 

	# read header information in the header array
	header_all = open(file,'rb').readlines()[0:6]
	header = []
	for line in header_all:
		line_split = line.split()
		header.append(line_split[1])
	# read the map in the map matrix
	map = np.genfromtxt(file, dtype=int, skip_header=6, filling_values='-1')
	
	print("map \n%s" %map)
	return header, header_all,	map

#-------------------------------------------------------------------------------------	
#	Determine the coordinates of neighboring cells of cell (col, row)
#-------------------------------------------------------------------------------------
def getNbh(col, row, ncols, nrows, four_neighbours):
	"""Determine the neighboring cells of the cell (col,row) and
		return the coordinates as arrays separated in nbhs_col and nbhs_row.
		The combination of the elements gives the coordinates of the neighbouring cells.
		
		input:
			col and row are coordinates of the reviewed element
			ncols, nrows are numbers of rows and columns in the map
			four_neighbours if True than neighboring cells are scanned else 8
	"""
	
	# assuming that a cell in the center has 8 neighbouring cells
	if four_neighbours == 'False':
		# cell is no edge cell
		if col > 0 and row > 0 and row < nrows -1 and col < ncols -1:
			nbhs_col = [x + col for x in[-1, -1, -1,  0, 0,	 1, 1, 1]]
			nbhs_row = [x + row for x in[-1,  0,  1, -1, 1, -1, 0, 1]]
		# cell is a left edge element but no corner element
		elif col == 0 and row > 0 and row < nrows -1:
			nbhs_col= [x + col for x in[0, 1, 1, 0, 1]]
			nbhs_row= [x + row for x in[-1, -1, 0, 1, 1]]	
		# cell is a right edge element but no corner element
		elif col == ncols -1 and row > 0 and row < nrows -1:
			nbhs_col= [x + col for x in[-1, -1, -1,	 0, 0]]
			nbhs_row= [x + row for x in[-1,	 0,	 1, -1, 1]]
		# cell is an upper edge element but no corner element
		elif row == 0 and col > 0 and col < ncols -1:
			nbhs_col= [x + col for x in[-1, -1,	 0, 1, 1 ]]
			nbhs_row= [x + row for x in[ 0,	 1, 1, 0, 1 ]]
		# cell is a bottom edge element but no corner element	
		elif row == nrows -1 and col > 0 and col < ncols -1:
			nbhs_col= [x + col for x in[-1, -1,	 0,	 1, 1 ]]
			nbhs_row= [x + row for x in[ -1, 0, -1, -1, 0 ]] 
		# cell is in the left upper corner
		elif col == 0 and row == 0:
			nbhs_col= [x + col for x in[ 0, 1, 1]]
			nbhs_row= [x + row for x in[ 1, 0, 1]]
		# cell is in the left bottom corner
		elif col == 0 and row == nrows -1:
			nbhs_col= [x + col for x in[ 0,	 1,	 1]]
			nbhs_row= [x + row for x in[ -1, 0, -1]] 
		# cell is in the right upper corner
		elif col == ncols -1 and row == 0:
			nbhs_col= [x + col for x in[ -1, -1, 0]]
			nbhs_row= [x + row for x in[  0,  1, 1]]
		# cell is in the right bottom corner
		else:
			nbhs_col= [x + col for x in[ -1, -1, 0 ]]
			nbhs_row= [x + row for x in[ -1,  0, -1]] 
			
	# assuming that a cell in the center has neighbouring cells
	elif four_neighbours == 'True':
		# cell is no edge cell
		if col > 0 and row > 0 and row < nrows -1 and col < ncols -1:
			nbhs_col = [x + col for x in[-1,  0, 0, 1]]
			nbhs_row = [x + row for x in[ 0, -1, 1, 0]]
		# cell is a left edge element but no corner element
		elif col == 0 and row > 0 and row < nrows -1:
			nbhs_col= [x + col for x in[0, 1, 0]]
			nbhs_row= [x + row for x in[-1, 0, 1]]	
		# cell is a right edge element but no corner element
		elif col == ncols -1 and row > 0 and row < nrows -1:
			nbhs_col= [x + col for x in[-1,	 0, 0]]
			nbhs_row= [x + row for x in[ 0, 1, -1]]		
		# cell is an upper edge element but no corner element
		elif row == 0 and col > 0 and col < ncols -1:
			nbhs_col= [x + col for x in[-1, 0, 1]]
			nbhs_row= [x + row for x in[ 0, 1, 0]]
		# cell is an bottom edge element but no corner element	
		elif row == nrows -1 and col > 0 and col < ncols -1:
			nbhs_col= [x + col for x in[-1, 0,	1]]
			nbhs_row= [x + row for x in[ 0, -1, 0]] 
		# cell is in the left upper corner
		elif col == 0 and row == 0:
			nbhs_col= [x + col for x in[ 0, 1]]
			nbhs_row= [x + row for x in[ 1, 0]]
		# cell is in the left bottom corner
		elif col == 0 and row == nrows -1:
			nbhs_col= [x + col for x in[ 0,	 1]]
			nbhs_row= [x + row for x in[ -1, 0]] 
		# cell is in the right upper corner
		elif col == ncols -1 and row == 0:
			nbhs_col= [x + col for x in[ -1, 0]]
			nbhs_row= [x + row for x in[  0, 1]]
		# cell is in the right bottom corner
		else:
			nbhs_col= [x + col for x in[ -1, 0 ]]
			nbhs_row= [x + row for x in[  0, -1]]

	else:
		msg = "Error: ini input for four_neighbours is not correct. Please check."
		WriteLogMsg(msg) 
		raise SystemError("Error: ini input for four_neighbours is not correct")
		req.close_window

	return [nbhs_row, nbhs_col]

#-------------------------------------------------------------------------------------	
#	Determination of patch elements
#-------------------------------------------------------------------------------------
def determine_patch_elements(row, col, map, patch_map, patch_ID, cls, four_neighbours):
	"""This recursive function scans all patch elements 
		and returns the coordinates of these elements.
		
		input:
			col and row are coordinates of the parent element
			map is the original ascii map
			patch_map is a map with patch_IDs for each patch element
			patch_ID is the ID of the new patch
			cls is the land use index of the patch
			four_neighbours if True than 4 neighboring cells are scanned else 8
	"""
	# determine coordinates of neighboring cells
	new_nbhs_row, new_nbhs_col	= getNbh(col, row, map.shape[1], map.shape[0], four_neighbours)
	# stack for patch elements whose neighboring cells should be determined
	nbhs_row = []
	nbhs_col = []
	for i in range(len(new_nbhs_row)):
		# add new neighboring cells to nbhs_row/col if new cells belong to cls and are not jet marked as patch element
		# the cell is no patch element if it has another land use id
		if map[new_nbhs_row[i], new_nbhs_col[i]] == cls and patch_map[new_nbhs_row[i], new_nbhs_col[i]] == 0:
			nbhs_row.append(new_nbhs_row[i])	 
			nbhs_col.append(new_nbhs_col[i])  
	while len(nbhs_row) > 0:
		# cells could be double in nbhs_row/col
		if patch_map[nbhs_row[0], nbhs_col[0]] == 0:
			# mark all patch elements in patch_map with patch_ID
			patch_map[nbhs_row[0], nbhs_col[0]] = patch_ID											
			# get coordinates of neighboring cells of this cell
			new_nbhs_row, new_nbhs_col	= getNbh(nbhs_col[0], nbhs_row[0], map.shape[1], map.shape[0], four_neighbours)
			for i in range(len(new_nbhs_row)):
				# add new neighboring cells to nbhs_row/col if new cells belong to cls and are not jet marked as patch element
				if map[new_nbhs_row[i], new_nbhs_col[i]] == cls and patch_map[new_nbhs_row[i], new_nbhs_col[i]] == 0:
					nbhs_row.append(new_nbhs_row[i])	 
					nbhs_col.append(new_nbhs_col[i])
		# delete this checked neighboring cell of the array 
		del nbhs_row[0]
		del nbhs_col[0]

	return patch_map

#-------------------------------------------------------------------------------------	
#	Determine patch elements of a patch ID map and check equality of land use index
#-------------------------------------------------------------------------------------
def determine_IDmap_patch_elements(row, col, patch_map, map, neighbors, cls, landuse, error_dic, four_neighbours):
	"""This recursive function scans all patch elements of the patch ID map,
		check if all patch elements have the same land use index
		and return the coordinates of the patch elements.
		
		input:
			col and row are coordinates of the parent element
			patch_map is the given patch ID map
			map is the original ascii map
			neighbors is a matrix for marking the scanned cells
			cls is the patch index
			landuse is the index of the first scanned patch element
			four_neighbours if True than 4 neighboring cells are scanned else 8
	"""
	# determine coordinates of neighboring cells
	new_nbhs_row, new_nbhs_col	= getNbh(col, row, map.shape[1], map.shape[0], four_neighbours)
	# stack for patch elements whose neighboring cells should be determined
	nbhs_row = []
	nbhs_col = []
	for i in range(len(new_nbhs_row)):
		# add new neighboring cells to nbhs_row/col if new cells belong to cls in patch_map and are not jet marked as scanned
		if patch_map[new_nbhs_row[i], new_nbhs_col[i]] == cls and neighbors[new_nbhs_row[i], new_nbhs_col[i]] == 0:
			nbhs_row.append(new_nbhs_row[i])
			nbhs_col.append(new_nbhs_col[i])
	# print ("nbhs_row, nbhs_col von (%s,%s): %s, %s" %(row, col, nbhs_row, nbhs_col)) 
	while len(nbhs_row) > 0:
		# if cell was not scanned before
		if neighbors[nbhs_row[0], nbhs_col[0]] == 0:
			# mark this cell in neighbors with True (scanned)
			neighbors[nbhs_row[0], nbhs_col[0]] = 1
			# and check if land use type for all patch elements is equal
			if map[nbhs_row[0], nbhs_col[0]] != landuse:
				error_dic.update({"(%s,%s)" %(nbhs_row[0], nbhs_col[0]) : "more than one land use index for one patch" })
			# determine coordinates of neighboring cells from this cell
			new_nbhs_row, new_nbhs_col	= getNbh(nbhs_col[0], nbhs_row[0], map.shape[1], map.shape[0], four_neighbours)
			for i in range(len(new_nbhs_row)):
				# add new neighboring cells to nbhs_row/col if new cells belong to cls in patch_map and not jet marked as scanned
				if patch_map[new_nbhs_row[i], new_nbhs_col[i]] == cls and neighbors[new_nbhs_row[i], new_nbhs_col[i]] == 0:
					nbhs_row.append(new_nbhs_row[i])
					nbhs_col.append(new_nbhs_col[i])		
		# delete this checked neighboring cell of the array 
		del nbhs_row[0]
		del nbhs_col[0]

	return neighbors, error_dic

#-------------------------------------------------------------------------------------	
#	Cluster the cells of the map into patches
#-------------------------------------------------------------------------------------
def create_patch_ID_map(map, NODATA_value, static_elements, four_neighbours):
	"""This function clusters the cells of the original map into patches
		and returns a patch ID map as a 2 dimensional array and the start individual as vector.
	
		input: 
			map is the original ascii map
			NODATA_value is the NODATA_value of the original map
			static_elements are the land use indices excluded from the optimization
			four_neighbours if True than 4 neighboring cells are scanned else 8
	"""
	
	msg = "Cluster cells of the original map into patches ..."
	WriteLogMsg(msg) 
	begin = time.time()
	
	patches= np.zeros([map.shape[0], map.shape[1]], int)
	ids = 0
	NoData = int(NODATA_value)
	genom = []
	# loop over all cells
	for row in range(0, map.shape[0]):
		if (row + 0.0) % (round(map.shape[0] / 10.0)) == 0 and (time.time() - begin) > 2 :
			progress = ((row+0.0) / map.shape[0]) * 100 
			WriteLogMsg("progress for scanning the rows in percent: %s" %progress)
		for col in range(0, map.shape[1]):
			# patchID = 0 used for static_elements
			# map element was not scanned before as patch element and is not a static element or the NODATA_value
			if patches[row,col]==0 and static_elements.count(map[row, col])==0 and map[row,col]!=NoData:
				cls = map[row, col]
				# increment scanned patch ID
				ids += 1
				# marke this cell as scanned patch element 
				patches[row, col] = ids
				determine_patch_elements(row,col, map, patches, ids, cls, four_neighbours) 
				# add the map cell value to the individual vector
				genom.append(cls)
	
	end = time.time()
	WriteLogMsg("Runtime for the generation of the patch_ID_map: %d seconds." %(end-begin))						
	return patches, genom	

#-------------------------------------------------------------------------------------	
#	Read the patch ID map and create the start individual
#-------------------------------------------------------------------------------------
def read_patch_ID_map(file, map, NODATA_value, static_elements, four_neighbours):
	"""This function reads a given patch ID map, checks its plausibility
		and returns the patch ID map as a 2 dimensional array and the start individual as vector.
	
		input: 
			file with the patch ID map (ascii format)
			map is the original ascii map
			NODATA_value is the NODATA_value of the original map
			static_elements are the land use indices excluded from the optimization
			four_neighbours if True than 4 neighboring cells are scanned else 8
	"""
	
	msg = "Read the patch ID map ..."
	WriteLogMsg(msg) 
	
	# transform map into a matrix
	try:
		patches = np.genfromtxt(file, dtype=int, skip_header=6, filling_values='-1')
	except ValueError as e:
		msg = "Error: Number of values in one row or column is not correct. Please check.\n %s" %e
		WriteLogMsg(msg) 
		raise SystemError("Error: Number of values in one row or column is not correct. Please check.\n %s" %e)
		req.close_window
	# check if number of rows and columns are equal to the original map
	if (map.shape[0] != patches.shape[0]) or (map.shape[1] != patches.shape[1]):
		msg = "Error: Number of rows or columns of the original map (rows: %s, columns: %s) and the patch ID map (rows: %s, columns: %s) are not equal. Please check." %(map.shape[0],map.shape[1],patches.shape[0],patches.shape[1])
		WriteLogMsg(msg) 
		raise SystemError("Error: Number of rows or columns of the original map and the patch ID map are not equal.")
		req.close_window
	# check that land use indices of the patch elements are equal 
	# and that the land use index is not a static element
	# if all checks are okay then add the land use index to the genome of the start individual
	max_patchID = patches.max()
	genom = np.zeros(max_patchID, int)
	help_map = np.zeros([map.shape[0], map.shape[1]], bool) # default is 0/False
	error_dic = {}
	for row in range(0, help_map.shape[0]):
		for col in range(0, help_map.shape[1]):
			# map element was not scanned before
			if help_map[row,col]==0:
				# first case: patch element = 0 
				# check if the index in the original map is a static element or NODATA_value
				if patches[row,col]==0:
					if static_elements.count(map[row, col])==0 and map[row,col]!=int(NODATA_value):
						error_dic.update({"(%s,%s)" %(row, col) : "zero but no static element and not a NODATA_value" })
				# patch element != 0
				else:
					# check that cell is not a static element or NODATA_value
					if static_elements.count(map[row, col])!=0 or map[row,col]==int(NODATA_value):
						error_dic.update({"(%s,%s)" %(row, col) : "non-zero but static element or NODATA_value" })
					# if file_HRU == None then check that patch ID was not registered in the genome before
					elif cfg.mapConfig.file_HRU == 'None' and genom[patches[row,col]-1] != 0:
						error_dic.update({"(%s,%s)" %(row, col) : "non-zero but patch ID is used for more than one patch" })
					# if file_HRU != None then check if land use of the two patches are equal
					elif cfg.mapConfig.file_HRU != 'None' and genom[patches[row,col]-1] != 0 and genom[patches[row,col]-1] != map[row,col]:
						error_dic.update({"(%s,%s)" %(row, col) : "non-zero but patch ID is used for more than one land use" })
					# plausibility checks are okay, then
					else:
						# add the land use index to the genome
						if genom[patches[row,col]-1] == 0:
							genom[patches[row,col]-1] = map[row,col]
						# scan all patch elements
						# check if land use is equal
						# and mark them in the help_map					
						determine_IDmap_patch_elements(row, col, patches, map, help_map, patches[row,col], genom[patches[row,col]-1], error_dic, four_neighbours)
				# mark the scanned cell in help_map 
				# so the cell is not scanned again by the determine_IDmap_patch_elements
				help_map[row,col] = 1 
	# if at least one plausibility check was not succesfull
	# print the error messages and stop the program
	if len(error_dic) != 0:
		msg = "Error: Some cells of the patch ID map break the plausibility rules. Please check."
		WriteLogMsg(msg)
		for key in error_dic:
			msg = "%s %s" %(key,error_dic[key])
			WriteLogMsg(msg) 
		raise SystemError("Error: Some cells of the patch ID map break the plausibility rules.")
		req.close_window	
							  
	return patches, genom.tolist() 

#-------------------------------------------------------------------------------------	
#	Check that original map and transition matrix are correct
#-------------------------------------------------------------------------------------
def check_matrices(map, trans_matrix):
	"""This function checks the trans_matrix for the following conditions to be true:
		- all diagonal elements are ones
		- all elements except the first column and row are only ones and zeros
		- first row and column elements are equal
		- all indices of the map are included in the transition matrix 
		
		This function checks the map for the following conditions to be true::
		- no elements are set by the default -1
		
		If one condition is not fulfill the program ends with an error message.
		
		input:
			map is the original ascii map
			trans_matrix holds the land use transition constraints
	""" 
	if len(trans_matrix) != 0:
		# check if all diagonal elements are ones
		for row in range(1,trans_matrix.shape[0]):
			if trans_matrix[row][row] != 1:
				msg = "Error: Not all diagonal elements of the transition matrix are ones. Please check."
				WriteLogMsg(msg) 
				raise SystemError("Error: Not all diagonal elements of the transition matrix are ones.")
				req.close_window
		# check if all elements except the first column and row are only ones and zeros
		if not (len(np.unique(trans_matrix[1:,1:]))==2 and 1 in np.unique(trans_matrix[1:,1:]) and 0 in np.unique(trans_matrix[1:,1:])):
			msg = "Error: Not all elements of the transition matrix (excluding first column and row) are ones and zeros. Please check."
			WriteLogMsg(msg) 
			raise SystemError("Error: Not all elements of the transition matrix (excluding first column and row) are ones and zeros.")
			req.close_window  
		# check if first row and column elements are equal
		for row in range(1,trans_matrix.shape[0]):
			if trans_matrix[row][0] != trans_matrix[0][row]:
				msg = "Error: Elements of first row and column are not equal. Please check."
				WriteLogMsg(msg) 
				raise SystemError("Error: Elements of first row and column are not equal.")
				req.close_window		
		# no elements are set by the default -1
		if (-1 in np.unique(trans_matrix)):
				msg = "Error: Elements of the transition matrix are -1 (default value). Please check."
				WriteLogMsg(msg) 
				raise SystemError("Error: Elements of the transition matrix are -1 (default value).")
				req.close_window 

		# check if all indices of the map are included in the transition matrix
		for element in np.unique(map):
			if element not in np.unique(trans_matrix):
				msg = "Error: Map element %s is not included in the transition matrix. Please check." % element
				WriteLogMsg(msg) 
				raise SystemError("Error: Map element %s is not included in the transition matrix." %element)
				req.close_window  
	
	# no elements are set by the default -1
	if (-1 in np.unique(map)):			
		msg = "Error: Elements of the original map are -1 (default value). Please check."
		WriteLogMsg(msg) 
		raise SystemError("Error: Elements of the original map are -1 (default value).")
		req.close_window 
	# no elements of the original map are zero (no data < -1 and land use id > 0)
	# 0 reserved for patch id map (patches which are excluded from the optimization)
	if (0 in np.unique(map)):			
		msg = "Error: Elements of the original map are 0. No data should be < -1 and  land use classes > 0. Please check."
		WriteLogMsg(msg) 
		raise SystemError("Error: Elements of the original map are 0. No data should be < -1 and  land use classes > 0. Please check.")
		req.close_window 
				
#-------------------------------------------------------------------------------------	
#	Determine all classes which are excluded from optimization 
#-------------------------------------------------------------------------------------
def determine_static_classes(trans_matrix, max_range):
	"""This function determines all classes which are excluded from optimization (static elements) 
		and returns arrays with the indices of static and non static elements.
		
		input:
			trans_matrix holding the land use transition constraints
			max_range is the maximum number of possible land use options
	""" 
	
	msg = "Determine static elements of the transition matrix ..."
	WriteLogMsg(msg)
	
	# identify all classes where column and row elements are zero (apart from diagonal elements)
	# that means that all patches of these classes cannot be converted
	static_elements = []
	nonstatic_elements = []
	# filter columns which fulfill the condition
	ones = 0
	# row-wise check for ones
	for row in range(1,trans_matrix.shape[0]):
		for col in range(1,trans_matrix.shape[1]):
			if trans_matrix[row][col] == 1:
				ones += 1
		# mark the candidate as static or non static element (row is checked)
		# if ones for row = 1 or max_range < land use index of trans_matrix
		if ones==1 or trans_matrix[row][0] > max_range:
			static_elements.append(trans_matrix[row][0])
		else:
			nonstatic_elements.append(trans_matrix[row][0])
		ones = 0
			
	# column-wise check for ones
	ones = 0
	index = 0
	if len(static_elements) != 0:
		for col in range(1,trans_matrix.shape[1]):
			if index < len(static_elements) and static_elements[index] <= max_range:
				if trans_matrix[0][col] == static_elements[index]:
					for row in range(1,trans_matrix.shape[0]):
						if trans_matrix[row][col] == 1:
							ones += 1
					if ones!=1:
						# index remains as it is for the next round 
						# because of length reduction from static_element
						nonstatic_elements.append(static_elements[index])
						del static_elements[index]
					else:
						index += 1
					ones = 0
				if len(static_elements) == 0:
					break
	
	return static_elements, nonstatic_elements
	
#------------------------------------------------------------------------------
#	Read the start individual and area proportion of each HRU from the HRU input file
#------------------------------------------------------------------------------
def read_HRUs(file_HRU):
	"""The function counts the HRUs listed in the HRU input file 
		and returns the result.

		input:
			filename of the HRU list		
	"""
	
	msg = "Read HRU file ..."
	WriteLogMsg(msg) 
	
	# array for area percentage and name of the HRUs
	global map_proportion
	
	# set default value for areas excluded from optimization
	map_proportion.append(float(0))

	genom = []
	
	reader = csv.reader(open(os.path.join(wrkDir, 'input', file_HRU), "r")) 
	i = 0
	sum_area = 0
	for row in reader: 
		if i != 0:
			map_proportion.append(float(row[2]))
			sum_area += float(row[2])
			genom.append(int(row[0]))
		i += 1
	if round(sum_area,1) <= 100.0:
		map_proportion[0] = 100 - sum_area
	else:
		msg = "Error: The input data for Area_rel are %s instead of 100 percent. Please check." %sum_area
		WriteLogMsg(msg) 
		raise SystemError("Error: The input data for Area_rel are %s instead of 100 percent. Please check." %sum_area)
		req.close_window
	
	# save start individual in global variable
	global start_individual 
	start_individual = genom
			
	return genom

#-------------------------------------------------------------------------------------	
#	Generate the genome of the start individual
#-------------------------------------------------------------------------------------
def generate_genom(max_range, file_HRU, map_file, trans_file, patchIDmap_file, four_neighbours, return_only_nonstatic):
	"""The function generates and returns the start individual and the non static land use indices
		based on an HRU file or ASCII map.
		It is called from generate_parameter in optiAlgorithm.py 
		and calls functions of maphandler.py for the individual generation.
	
		input data: 
				max_range is the maximum number of possible land use options for the individual
				file_HRU file with the original HRU list 
				map_file file with the original ascii map
				trans_file file with the land use transition constraints
				patchIDmap_file file with the patch ID map (ascii format)
				four_neighbours if True than 4 neighboring cells are scanned else 8
	"""

	global map
	global header
	global header_all
	global trans_matrix
	global static_elements
	global nonstatic_elements
	global patchID_map
	global map_proportion
	# dictionary for possible land use options per element
	global possible_elements

	# HRU and given patch ID map need an original map
	if file_HRU != 'None' and (map_file != 'None' or patchIDmap_file != 'None'):
		if map_file == 'None' or patchIDmap_file == 'None':
			msg = "Error: If you want to use maps with HRUs than you need an original map and a patch ID map. Please check."
			WriteLogMsg(msg) 
			raise SystemError("Error: If you want to use maps with HRUs than you need an original map and a patch ID map. Please check.")
			req.close_window

	# if HRU file is given
	if file_HRU != 'None':		
		msg = "Generate start individual based on an HRU file ..."
		WriteLogMsg(msg)
		genom = read_HRUs(file_HRU)
		
	# if ASCII map is given
	if map_file != 'None': 
		if file_HRU == 'None':
			msg = "Generate start individual based on an ASCII map..."
			WriteLogMsg(msg)				
		# read the original ascii map
		header, header_all, map = read_ascii_map(os.path.join(wrkDir, 'input', map_file))

	# no ASCII original map or HRU list is given	
	if map_file == 'None' and file_HRU == 'None':
		msg = "Error: No original ASCII map or HRU list is given by the config.ini. Please check."
		WriteLogMsg(msg) 
		raise SystemError("Error: No original ASCII map or HRU list is given by the config.ini.")
		req.close_window
		
	# if a transition matrix is given
	if trans_file != 'None':
		# read the transition matrix for land use change
		msg = "Read transition matrix..."
		WriteLogMsg(msg) 
		# write information into a matrix
		trans_matrix = np.genfromtxt(os.path.join(wrkDir, 'input',trans_file), dtype=int, filling_values='-1')
		# check if max_range is plausible
		if max_range > trans_matrix.max():
			msg = "Error: Input of max_range is bigger than the maximum element of the transition matrix. Please check."
			WriteLogMsg(msg) 
			raise SystemError("Error: Input of max_range is bigger than the maximum element of the transition matrix.")
			req.close_window
		if map_file != 'None':
			if max_range < map.max():
				msg = "Error: Input of max_range is smaller than the maximum element (%s) in the original ASCII map. Please check." %map.max()
				WriteLogMsg(msg) 
				raise SystemError("Error: Input of max_range is smaller than the maximum element in the original ASCII map.")
				req.close_window
			  
		# check format requirements
		check_matrices(map, trans_matrix)

		# determine static land use elements
		static_elements, nonstatic_elements = determine_static_classes(trans_matrix, max_range)

				
		# possible land use options for each land use class according to transition matrix
		if len(possible_elements) == 0:
			for i in range(1,cfg.modelConfig.max_range+1):
				if i not in static_elements:
					# determine indices of transition matrix where in row with first element = i are ones
					indices = np.asarray(np.nonzero(trans_matrix[np.nonzero(np.unique(trans_matrix[:,:1]) == i)[0][0],:])[0])
					# delete the index for 1 as land use name
					if 0 in indices:
						indices = np.delete(indices,0)
					# determine the land use classes of the indices
					if len(indices) !=0 and trans_matrix[0][indices[0]] <= cfg.modelConfig.max_range:
						land_use_classes = np.array([trans_matrix[0][indices[0]]])
						for k in range(1,len(indices)):
							if trans_matrix[0][indices[k]] <= cfg.modelConfig.max_range:
								land_use_classes = np.append(land_use_classes,trans_matrix[0][indices[k]])
					# save the possible land use classes
					possible_elements.update({ i : land_use_classes})
					
			WriteLogMsg("possible_elements per land use: ")
			for key in possible_elements:
				WriteLogMsg("land use, elements: %s,%s" %(key,possible_elements[key]))
		
	# if no transition matrix is given -> land use change is completely random 
	else:
		# check format requirements
		if map_file != 'None': 
			check_matrices(map, trans_matrix)
		# determine non_static elements ( 1 <= x <= max_range)
		for i in range(1,max_range+1):
			nonstatic_elements.append(i)
		
	# header[5] is the NODATA_value of the original ascii map
	# read patchID_map if it is available
	if patchIDmap_file != 'None':
		# read the patch ID map, check plausibility and create the start individual
		if file_HRU == 'None':
			msg = "Create the start individual based on the ID map ..."
			WriteLogMsg(msg)
		patchID_map, genom = read_patch_ID_map(os.path.join(wrkDir, 'input', patchIDmap_file), map, header[5], static_elements, four_neighbours)
	# the patch ID map cannot be created based on the original map for HRUs (one HRU can consist of more than one patch)
	elif file_HRU == 'None':
		# create the patch ID map based on the original map and the start individual
		patchID_map, genom = create_patch_ID_map(map, header[5], static_elements, four_neighbours)
	# log the patch ID map
	if file_HRU == 'None' or (file_HRU != 'None' and patchIDmap_file != 'None'): 
		WriteMap(header_all, patchID_map)
		
	# analyze the proportional area of the patches in %
	# for HRUs the map_proportion is obtained in read_HRUs 
	if file_HRU == 'None':
		unique, counts = np.unique(patchID_map, return_counts = True)		
		number_cells = patchID_map.shape[0]*patchID_map.shape[1]
		if 0 not in unique:
			# no cells are excluded from optimization
			map_proportion.append(float(0))
		for item in counts: 
			map_proportion.append(float(item)/float(number_cells)*100)
		# check if sum of area proportions is 100%
		sum = 0
		for i in map_proportion:
			sum = sum + i
		if round(sum) != 100:
			msg = "Error: The sum of patch/HRU areas is %s instead of 100 percent. Please check." %sum
			WriteLogMsg(msg) 
			raise SystemError("Error: The sum of patch/HRU areas is %s instead of 100 percent." %sum)
			req.close_window
	
	# check all indices of the map are included in the transition matrix
	if trans_file != 'None':
		for element in np.unique(genom):
			if element not in np.unique(trans_matrix):
				msg = "Error:  Land use class %s is not included in the transition matrix. Please check." % element
				WriteLogMsg(msg) 
				raise SystemError("Error: Land use class %s is not included in the transition matrix." %element)
				req.close_window
	
	# log the start individual and the non static land use indices
	if len(genom) > 100:
		msg = "start individual of the input data with length %d" %len(genom)
		WriteLogMsg(msg)
	else:
		msg = "start individual of the input data with length %d: %r" %(len(genom),genom)
		WriteLogMsg(msg)
	msg = "DiscreteBounder values: %r" %nonstatic_elements
	WriteLogMsg(msg)	
	
	# save start individual in global variable
	global start_individual 
	start_individual = genom

	WriteLogMsg("map proportion of the patches: %r" %map_proportion)
	
	if return_only_nonstatic:
		return nonstatic_elements
	
	return genom, nonstatic_elements

#------------------------------------------------------------------------------	 
#	Check if new_cand is element of the tabu memory
#------------------------------------------------------------------------------
def check_impossible_candidates(new_cand):
	"""Return True if new_cand is element of the tabu memory
		else False.
		
		input data: 
			new_cand is a part of a new candidate
	"""
	
	# array with segments of impossible candidates (tabu memory)
	global impossible_cand
	# array for the start individual
	global start_individual

	return_value = False
	copy_new_cand = new_cand[:] 
	# check complete candidate
	if copy_new_cand in impossible_cand:
		return_value = True
	# check if first part is in impossible_cand:
	else:
		for item in impossible_cand:
			if len(item) < len(copy_new_cand):
				for index in range(len(item)):
					# item not first part of new_cand
					if item[index] != copy_new_cand[index]:
						break
					# item is first part of new_cand -> new_cand added to tabu memory
					elif index == (len(item)-1) and item[index] == copy_new_cand[index]:
						return_value = True

	return return_value

#------------------------------------------------------------------------------	 
#	Check if created individual is plausible
#------------------------------------------------------------------------------
def individual_filter(new_cand):
	"""Check if the created individual violates the constraints. 
		Return True if the individual is feasible else False.
		
		The function is structured for changeable objects for comparison but at 
		the moment only the original start individual is used for the comparison 
		and the decision if the new candidate is feasible or not.
		
		input data: 
			new_cand is new candidate for checking the feasibility
	"""
	
	return_value = True
	# array for land use transition constraints
	global trans_matrix
	# array for the area percentages of the patches
	global map_proportion
	# array for total area constraints
	global min_max_diff
	# array for the start individual
	global start_individual
	# array for static land use types
	global static_elements
	
	compare_individual = start_individual

	i = 0
	if cfg.ea.start_from_previous_gen == True and i < cfg.ea.pop_size:
		i =i +1
		return_value= True
	else:
		# for each element of the individual, check if transition is allowed
		if cfg.mapConfig.file_transformation != 'None':				
			for i in range(0,len(new_cand)):
				if trans_matrix[np.nonzero(np.unique(trans_matrix[:,:1]) == compare_individual[i])[0][0]][np.nonzero(trans_matrix[0] == new_cand[i])[0][0]] != 1:
					return_value = False
					break
	
	# check if total area constraints are satisfied
	if cfg.mapConfig.file_difference != 'None':
		for i in np.unique(new_cand):
			sum_new = 0 
			for m in np.nonzero(new_cand==i)[0]:
				# index m+1 because of the exclusion of patch ID 0 from optimization
				sum_new += map_proportion[m+1]
			if not (min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] <= sum_new <= min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]):
				return_value = False
				break
		
		# check if area rules for land use types in new_cand do not apply -> min should be zero
		for i in range(1,cfg.modelConfig.max_range+1):
			if (i not in static_elements) and (i not in np.unique(new_cand)):
				if min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] > 0:
					return_value = False
					break
		
	return return_value

#------------------------------------------------------------------------------	 
#	Determine possible land use classes for next new_cand element
#------------------------------------------------------------------------------
def determine_possible_elements(possible_land_use, new_cand):
	"""This function is used during repair mutation and returns the 
		possible land use classes for the next new_cand element taking 
		into account the tabu memory.
		
		input data: 
			possible_land_use is an array with possible land use classes 
								 derived from the transition matrix
			new_cand is a part of a new candidate
	"""
	
	# array with segments of impossible candidates (tabu memory)
	global impossible_cand

	copy_poss_el = np.copy(possible_land_use)
	copy_new_cand = new_cand[:] 
	# create all new candidate parts which can be created with one of the possible_elements and the new_cand
	# and check if the new one is excluded by the tabu memory
	for element in possible_land_use:
		if len(copy_new_cand) == len(new_cand):
			copy_new_cand.append(element)
		else:
			copy_new_cand[len(copy_new_cand)-1] = element
		if copy_new_cand in impossible_cand:
			copy_poss_el = np.delete(copy_poss_el,np.nonzero(copy_poss_el == copy_new_cand[len(copy_new_cand)-1])[0][0])
	
	return copy_poss_el

#------------------------------------------------------------------------------	 
#	Determine possible land use classes for next new_cand element
#------------------------------------------------------------------------------
def determine_possible_landuse(possible_land_use, index, new_cand):
	"""This function is used during extreme seed generation and returns 
		the possible land use classes for an index of the candidate taking 
		into account the tabu memory.
		
		input data: 
			possible_land_use is an array with possible land use classes 
								 derived from the transition matrix
			index is the index which should be checked
			new_cand is a part of a new candidate
	"""
	
	# array with segments of impossible candidates (tabu memory)
	global impossible_cand

	copy_poss_el = np.copy(possible_land_use)
	copy_new_cand = new_cand[:] 
	# create all new candidate parts which can be created with one of the possible_elements and the new_cand
	# and check if the new one is excluded by the tabu memory
	for element in possible_land_use:
		copy_new_cand[index] = element
		if copy_new_cand in impossible_cand:
			copy_poss_el = np.delete(copy_poss_el,np.nonzero(copy_poss_el == element)[0][0])
	
	return copy_poss_el


#------------------------------------------------------------------------------	 
#	Check tabu memory for redundant entries
#------------------------------------------------------------------------------
def check_redundancy(new_cand):
	"""Return tabu memory with new_cand as additional element 
		but without redundant entries.
		
		example: if [1,1,3,4] is excluded and now [1,1,3] too, than [1,1,3,4] is 
				 not relevant anymore 
		
		input data: 
			new_cand is a part of an excluded candidate
			impossible_cand is the tabu memory
	"""
	
	# dictionary for possible land use options per element
	global possible_elements
	# array with segments of impossible candidates (tabu memory)
	global impossible_cand
	# original start individual derived from input data
	global start_individual
	
	search_object = new_cand[:]
	sublist = []
	last_elements = []

	# check if new_cand is part of the impossible_cand
	# if True -> do nothing
	# if False -> check for redundant entries and add new_cand
	if not (check_impossible_candidates(search_object)):  
		# check all items with one element more than the part of the new candidate 
		# for redundancy 
		if len(search_object) != len(start_individual):
			for item in impossible_cand:
				check_var = True
				# check items with one more element than the new_cand
				if len(item) == (len(search_object)+1):
					# check if all elements except the last one are equal with the search_object
					for index in range(len(search_object)):
						if search_object[index] != item[index]:
							check_var = False
							break
					# save potential redundant items and the last elements in separate lists
					if check_var == True:
						sublist.append(item)
						last_elements.append(item[len(item)-1])
	
			# delete redundant items
			if cfg.mapConfig.file_transformation != 'None':
				if len(last_elements) == len(possible_elements[start_individual[len(search_object)]]):
					for item in sublist:
						impossible_cand.remove(item)
			else:
				if len(last_elements) == cfg.modelConfig.max_range:
					for item in sublist:
						impossible_cand.remove(item)

		# add the excluded part of candidate to the tabu memory
		impossible_cand.append(search_object)
		
	return impossible_cand

#------------------------------------------------------------------------------	 
#	Generate new candidate for optimization algorithm with filter function
#------------------------------------------------------------------------------
def filter_variator(candidate): 
	"""Generate new candidates with the filter method (filter_mutation) if the given candidate 
		was used before or the candidate is not plausible.
		The function generates randomly new individuals and checks them for feasibility.
		It returns the first individual which is conform with the constraints.
		
		input: 
			candidate is given by the optimization algorithm or the start individual 
						for the first population
	"""
	
	# array for total area constraints
	global min_max_diff
	# dictionary for possible land use options per element
	global possible_elements
	# original start individual derived from input data
	global start_individual
	
	begin = time.time()
	timeout = begin + 60*5
	max_loops = 1000000
	# new candidate list 
	new_cand = []
	# check if candidate is part of the tabu memory
	# if True -> create a new candidate
	# if False -> check feasibility -> if True then use the candidate else create a new candidate
	if check_impossible_candidates(candidate) == True or (check_impossible_candidates(candidate) == False and individual_filter(candidate) == False):
		# create individual according to land use transition and/or total area constraints
		if cfg.mapConfig.file_transformation != 'None' or  cfg.mapConfig.file_difference != 'None':
			WriteLogMsg("create individual in accordance with transition matrix and/or total area rules") 
			cand_accept = False
			count_loops = 0 

			while cand_accept == False:	 

				if count_loops % 100000 == 0:
					WriteLogMsg("count_loops: %s" %count_loops)
					WriteLogMsg("new_cand: %s" %new_cand) 
				if count_loops == 2:
					WriteLogMsg("count_loops=2: start_ind: %s" %candidate)
					WriteLogMsg("count_loops=2: new_cand: %s" %new_cand)

				del new_cand[:] 
				# for each element within the genome
				for i in range(len(start_individual)):
					if cfg.mapConfig.file_transformation != 'None':
						# generate a list with possible land use options according to transition matrix
						pos_el_new_cand = possible_elements[start_individual[len(new_cand)]]			
					else:
						# generate a list with possible land use options from whole range of options
						pos_el_new_cand = np.arange(1,cfg.modelConfig.max_range+1)
					# select one option from this list randomly
					clbrValue = random.choice(pos_el_new_cand.tolist())
					new_cand.append(clbrValue)
				# check feasibility of the new individual
				cand_accept = individual_filter(new_cand)
				if cand_accept == True:
					WriteLogMsg("%r generated individuals break the plausibility rules." % count_loops)
				count_loops += 1

		# create individuals purely randomly (neither transition nor area constraints are defined)
		else:
			WriteLogMsg("create individual purely randomly (neither transition nor area constraints are defined)") 
			for i in range(len(start_individual)):
				# generate a list with possible land use options from whole range of options
				clbrValue = random.randint(1, cfg.modelConfig.max_range)
				new_cand.append(clbrValue)
	# use the candidate given by the optimization algorithm
	else:
		new_cand = candidate[:]
	
	# add accepted new candidate to impossible candidates 
	check_redundancy(new_cand)
	end = time.time()
	WriteLogMsg("The generation of a new candidate needed %d seconds." %(end-begin))
	#WriteLogMsg("%r" % new_cand)
	return new_cand
		  
#------------------------------------------------------------------------------	 
#	Generate new candidate for optimization algorithm by rules based on the 
#	transformation and min/max information
#------------------------------------------------------------------------------
def logical_variator(candidate, first_generation='False'):	
	"""Generate feasible new candidates.
		Return the first individual which is conform with the defined constraints.
		
		input: 
			candidate is given by the optimization algorithm or the start individual 
						for the first population
			first_generation is only True for the first generation, later the new_cand 
						will be generated without priority 

	""" 

	# original start individual derived from input data
	global start_individual
	# array for transition constraints
	global trans_matrix
	# array for total area constraints
	global min_max_diff
	# dictionary for possible land use options per element
	global possible_elements
	# array for non-static land use types
	global nonstatic_elements
	# array for the area percentage of the patches
	global map_proportion
	# Boolean variable for starting the special_termination
	global end_optimization
	
	begin = time.time()

	# line for printing the checked new_cand and information of the generation process
	candidate_list = []

	# if termination of optimization process is started
	# then fill the remaining individuals of the population with default individual
	if end_optimization == True:
		new_cand = start_individual
	# check if candidate is part of the tabu memory
	# if True -> create a new candidate
	# if False -> check feasibility -> if True then use the candidate else create a new candidate
	elif check_impossible_candidates(candidate) == True or (check_impossible_candidates(candidate) == False and individual_filter(candidate) == False):
		# new candidate list 
		new_cand = [] 
		# create individual according to land use transition and/or total area constraints
		if cfg.mapConfig.file_transformation != 'None' or  cfg.mapConfig.file_difference != 'None':				
			WriteLogMsg("create individual based on tabu memory and in accordance with transition and total area constraints") 
		else:				
			WriteLogMsg("create individual only based on tabu memory (no constraints defined!)")				  

		# create new candidate	  
		count = 0	
		# priorisation of one special land use based on the plausibility checks
		next_position = -1
		# if last items of new_cand were deleted based on the plausibilty checks 
		# then reset the information of last critical index and excluded land uses
		ind_critical = -1
		excluded_landuse = []

		while len(new_cand) < len(start_individual):  
			count += 1
			# if signal_del_more = True then delete more than one position of the new_cand
			signal_del_more = False
			# signal that new_cand was changed during the plausibility checks
			change_new_cand = False
			# if last items of new_cand were deleted based on the plausibilty checks 
			# then reset the information of last critical index and excluded land uses
			if ind_critical >= len(new_cand):
				ind_critical = -1
				excluded_landuse = []	  

			# determine possible land use options for the next new_cand element
			if cfg.mapConfig.file_transformation != 'None':
				pos_el_new_cand = possible_elements[start_individual[len(new_cand)]]			
				option_elements = determine_possible_elements(pos_el_new_cand, new_cand)
			else:
				option_elements = determine_possible_elements(np.arange(1,cfg.modelConfig.max_range+1), new_cand)
			# append the next element of new_cand			  
			if len(option_elements) > 0:

				if count % 100000 == 0 or (count % 1000 == 0 and time.time()-begin > 300 and time.time()-begin < 600 ):
					WriteLogMsg("count: %s" %count)
					WriteLogMsg("new_cand: %s" %new_cand)
					# line for printing the checked new_cand and information of the generation process
					candidate_list.append("%s" %['preliminary result'])
					if cfg.ea.write_tabu_memory == True:
						WriteCandidateList(candidate_list) 

				# use the land use which was determined in a plausibility check before, if it is possible
				if next_position > 0 and next_position in option_elements:
					new_cand.append(next_position)	
				# try to use as much as possible of the candidate from the optimization algorithm
				elif cfg.ea.priority == 'True' and first_generation == 'False' and (candidate[len(new_cand)] in option_elements):
					new_cand.append(candidate[len(new_cand)])
				# choose randomly one of the possible land use classes for the next new_cand element
				else: 
					clbrValue = random.choice(option_elements.tolist())
					new_cand.append(clbrValue) 
				# delete old priority of a land use for the next plausibility checks 
				if next_position > -1:
					next_position = -1

				# line for printing the checked new_cand and information of the generation process
				wert = "%s" %new_cand
				candidate_list.append(wert)
 
				# check if total area constraints for the given land use are satisfied
				if cfg.mapConfig.file_difference != 'None':
					plausible_value = True
					diff_min = []
					diff_min_landuse = [] 
					dict_sum_landuse = {}
					for i in np.unique(new_cand):
						sum_new = 0
						for m in np.nonzero(new_cand==i)[0]:
							# index m+1 because of the exclusion of patch ID 0
							sum_new += map_proportion[m+1] 
						# save area of the given land use in dictionary dict_sum_landuse
						dict_sum_landuse.update({i:sum_new})
						# if total area exceeds the defined maximum 
						if sum_new > min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]:
							plausible_value = False

							# if only one land use class is possible at the last position of new_cand then look for a previous position 
							# which can be changed into another class to make the last position possible
							if cfg.mapConfig.file_transformation != 'None' and new_cand[-1] == start_individual[len(new_cand)-1] and len(possible_elements[i]) == 1:
								sum_mandatory = sum_new - min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]
								# determine the position that should be changed
								for k in range(2,len(new_cand)+1):
									ind_index = len(new_cand)-k
									if new_cand[ind_index] == i and start_individual[ind_index] != i:
										sum_mandatory -= map_proportion[ind_index+1]
										if sum_mandatory <= 0:
											# change new_cand
											test_cand = new_cand[:ind_index] 
											if check_impossible_candidates(test_cand) == False:
												# change new_cand and save the old one in tabu memory
												check_redundancy(new_cand[:ind_index+1])
												new_cand = test_cand[:]
												change_new_cand = True
												# line for printing the checked new_cand and information of the generation process
												candidate_list.append("%s" %['1a', i, 'landuse i > max and the last position is mandatory', ind_index, 'ind_index', 'last new_cand elements are deleted until index ind_index'])
												break
											else:
												# line for printing the checked new_cand and information of the generation process
												candidate_list.append("%s" %['information: case 1a with no better break rule. Please check the logical_variator at this point.'])
								if change_new_cand == False:
									# line for printing the checked new_cand and information of the generation process
									candidate_list.append("%s" %['information: case 1a and new_cand did not change. Please check the logical_variator at this point.'])

							# if last position of new_cand is marked as critical and number of excluded land use classes
							# is equal to the possible land use classes then look for the next position 
							# which is one of the excluded land uses but not mandatory
							elif cfg.mapConfig.file_transformation != 'None' and len(excluded_landuse) == len(possible_elements[start_individual[len(new_cand)-1]]):
								# determine the position which should be changed
								for k in range(2,len(new_cand)+1):
									ind_index = len(new_cand)-k
									if new_cand[ind_index] in excluded_landuse and not (new_cand[ind_index] == start_individual[ind_index] and len(possible_elements[start_individual[ind_index]]) == 1):
										# change new_cand
										test_cand = new_cand[:ind_index] 
										if check_impossible_candidates(test_cand) == False:
											# change new_cand and save the old one in the tabu memory
											check_redundancy(new_cand[:ind_index+1])
											new_cand = test_cand[:]
											change_new_cand = True
											# line for printing the checked new_cand and information of the generation process
											candidate_list.append("%s" %['1b', i, 'landuse i > max and all possible land uses for the last position are excluded', excluded_landuse, 'excluded_landuse', ind_index, 'ind_index', 'last new_cand elements are deleted until index ind_index'])
											break
										else:
											# line for printing the checked new_cand and information of the generation process
											candidate_list.append("%s" %['information: case 1b with no better break rule. Please check the logical_variator at this point.'])
								if change_new_cand == False:
									# line for printing the checked new_cand and information of the generation process
									candidate_list.append("%s" %['information: case 1b and new_cand did not change. Please check the logical_variator at this point.'])

							else:
								# save last index as critical and excluded land use
								ind_critical = len(new_cand)-1
								excluded_landuse.append(new_cand[-1]) 
								# line for printing the checked new_cand and information of the generation process
								candidate_list.append("%s" %[1, i, 'landuse i > max'])

							break

						# save difference between the total area in new_cand and the defined minimum
						if sum_new < min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]]:
							diff_min.append(min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]]-sum_new)
							diff_min_landuse.append(i)

						# check that sum_new is conform with the defined minimum area for this land use class
						if plausible_value == True and len(new_cand) == len(start_individual) and sum_new < min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]]:
							plausible_value = False
							# set priority for land use i
							next_position = i

							# line for printing the checked new_cand and information of the generation process
							candidate_list.append("%s" %[2, i, 'landuse i < min and individual complete'])

							break

					# line for printing the checked new_cand and information of the generation process
					if len(new_cand) == (len(start_individual)-1):
						candidate_list.append("%s" %[diff_min,'diff_min',diff_min_landuse,'diff_min_landuse'])
	
					# if new_cand is not excluded and the length of new_cand is equal to the length of the start individual
					# check min/max rules for non-static land use types which are not included in new_cand 
					# -> min should be zero
					if plausible_value == True and len(new_cand) == len(start_individual):
						for i in nonstatic_elements:
							if (i not in np.unique(new_cand)) and min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] > 0:
								plausible_value = False

								# line for printing the checked new_cand and information of the generation process
								candidate_list.append("%s" %[3, i, 'landuse i min > 0, but i not in new_cand'])

								break
					# if new_cand is not excluded and length of new_cand is smaller then the length of the start individual
					# check min area rules for non-static land use types
					elif plausible_value == True and len(new_cand) < len(start_individual):
						# add the missing area for these land use classes to the diff_min list
						for i in nonstatic_elements:
							# element is not in new_cand but min > 0
							if (i not in np.unique(new_cand)) and min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] > 0:
								diff_min.append(min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]])
								diff_min_landuse.append(i)
								# save total land use area in dictionary dict_sum_landuse
								dict_sum_landuse.update({i:0})
						# check if enough free patches exist and if the respective area could be enough 
						if len(diff_min) > (len(start_individual)-len(new_cand)) or sum(diff_min) > sum(map_proportion[len(new_cand)+1:]):
							plausible_value = False
							# set priority for first land use in diff_min_landuse
							for i in diff_min_landuse:
								if i in option_elements:
									next_position = i
							candidate_list.append("%s" %[next_position,'next_position',option_elements,'option_elements'])

							# line for printing the checked new_cand and information of the generation process
							if len(new_cand) != (len(start_individual)-1):
								candidate_list.append("%s" %[diff_min,'diff_min',diff_min_landuse,'diff_min_landuse'])
							candidate_list.append("%s" %[4, 'not enough patches left for missing land uses or not enough area left to reach required min area'])

						# check if missing land use areas are possible in the last part of new_cand
						# with respect to start_individual and transition rules
						elif cfg.mapConfig.file_transformation != 'None':
							index_landuse = 0
							for i in diff_min_landuse:
								sum_free = 0
								# determine possible area left for land use i
								for j in range(len(new_cand),len(start_individual)):
									if i in possible_elements[start_individual[j]]:
										sum_free += map_proportion[j+1]
								# if not enough area is left for land use i and i is not element of option_elements 
								# then try to go back to a previous position where land use i is not set but possible
								if sum_free < diff_min[index_landuse] and i not in option_elements:
									# determine the previous position where i is not set but possible
									# (begin search at the end of new_cand)
									sum_mandatory = sum_free
									for k in range(1,len(new_cand)+1):
										ind_index = len(new_cand)-k
										if i in possible_elements[start_individual[ind_index]]:
											sum_mandatory += map_proportion[ind_index+1]
											if sum_mandatory >= min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] and new_cand[ind_index] != i:
												# check if change at this position is not excluded by tabu memory
												test_cand = new_cand[:ind_index]
												test_cand.append(i)
												if check_impossible_candidates(test_cand) == False:
													# change new_cand and save the old one in tabu memory
													plausible_value = False
													del new_cand[-1]
													check_redundancy(new_cand)
													new_cand = test_cand[:]
													change_new_cand = True
													candidate_list.append("%s" %[diff_min,'diff_min',diff_min_landuse,'diff_min_landuse'])
													candidate_list.append("%s" %[5,i,'i',sum_free,'sum_free',ind_index,'ind_index','new_cand is changed at position ind_index to i and the rest of new_cand is deleted'])
													break
									# maybe in this case a new logical rule is also possible
									if change_new_cand == False:
										candidate_list.append("%s" %['information: check if a new additional logical rule is possible (no previous position for land use i possible; check impossible_cand)'])

								# missing land uses could be added if one or two possitions of new_cand will be deleted
								if sum_free < diff_min[index_landuse] and change_new_cand == False:
									plausible_value = False
									# set priority for first land use in diff_min_landuse which is in option_elements
									for k in diff_min_landuse:
										if k in option_elements:
											next_position = k

									# diff_min_landuse are not in option_elements
									if next_position == -1:
										# if signal_del_more = True then delete more than one position of the new_cand
										signal_del_more = True

									# line for printing the checked new_cand and information of the generation process
									candidate_list.append("%s" %[diff_min,'diff_min',diff_min_landuse,'diff_min_landuse'])
									candidate_list.append("%s" %[6,i,'i',sum_free,'sum_free','not enough area left for land use i within the rest of new_cand given the transition rules'])
									break
								index_landuse += 1			

					# if current new_cand is implausible and 
					# new_cand was not changed during the plausibility checks
					# then save new_cand as impossible candidate in tabu memory
					# and delete the last new_cand element
					if plausible_value == False and change_new_cand == False:
						# check if missing landuses are possible at the last position
						# if not delete the last two elements of new_cand
						if (len(diff_min_landuse) > len(start_individual)-len(new_cand)) or signal_del_more == True:
							# check if missing land use areas are possible
							del_more = False
							if signal_del_more != True:
								del_more = True
								for i in diff_min_landuse:
									if i in option_elements:
										del_more = False

							if del_more == True or signal_del_more == True:
								del new_cand[-1]

								# line for printing the checked new_cand and information of the generation process
								candidate_list.append("%s" %[7, option_elements,'option_elements', 'missing land uses not in option_elements for last position'])

						check_redundancy(new_cand)
						del new_cand[-1]
			
			# all possible individuals have been used -> terminate the optimization process
			elif len(option_elements) == 0 and len(new_cand) == 0:
				end_optimization = True
				WriteLogMsg("No more plausible individuals exist. The termination is started.") 
				# if termination of optimization process is started
				# then fill the rest of the population with start individual
				new_cand = start_individual
				break
			# no plausible land use class for the next element exists
			else:
				# add part of candidate to tabu memory and try another last element for new_cand
				check_redundancy(new_cand)
				del new_cand[-1]

				# line for printing the checked new_cand and information of the generation process
				candidate_list.append("%s" %[8, option_elements, 'option_elements is empty'])
			if cfg.ea.max_repair_trials > 0:
				if count % cfg.ea.max_repair_trials == 0:
					WriteLogMsg("Max trials for candidate generation have been reached. Continue with start individual")
					new_cand = start_individual

					break
			 
	# use the candidate given by the optimization algorithm
	else:
		new_cand = candidate[:]
	
	# add accepted new candidate to tabu memory 
	# to avoid that this candidate is used again as an individual
	if end_optimization != True:
		check_redundancy(new_cand)
		
	end = time.time()
	WriteLogMsg("The generation of a new candidate needed %d seconds." %(end-begin))
	if len(start_individual) < 100:
		WriteLogMsg("%r" % new_cand)


	# line for printing the checked new_cand and information of the generation process
	if end-begin >= 120:
		WriteCandidateList(candidate_list)

	return new_cand 

#------------------------------------------------------------------------------	 
#	Generate first parameter for optimization algorithm
#------------------------------------------------------------------------------	 
def create_extreme_seed(land_use, maximize):  
	"""Generate new candidate with extreme land cover for given land use according to transition 
		and total area constraints.
		Return the first individual satisfying the constraints.
	  
		input: 
			land_use is the land use class for extreme land cover
			maximize is the direction of the extreme land cover
	""" 

	# original start individual derived from input data
	global start_individual
	# array for transition constraints
	global trans_matrix
	# array for total area constraints
	global min_max_diff
	# dictionary for possible land use options per element
	global possible_elements
	# array for non-static land use classes
	global nonstatic_elements
	# array for the area percentage of the patches
	global map_proportion
	# dictionary with patch indices and total area per land use class
	global start_land_cover

	begin = time.time()

	# line for printing the checked new_cand and information of the generation process
	candidate_list = []

	candidate_list.append("land_use: %s, maximize: %s/n" %(land_use, maximize))

	# new candidate list 
	new_cand = start_individual[:] 
	# create individual according to transition and total area constraints
	# if no constraints or only transition constraints are defined
	# then replace all land use patches with other land use options for min
	# or replace all land use patches with the land use at hand for max
	# according to the transition constraints
	if cfg.mapConfig.file_difference == 'None':				
		WriteLogMsg("create individual with create_extreme_seed without total area constraints but transition matrix if given")
		if maximize == True:
			for i in range(len(start_individual)):
				if i not in start_land_cover[land_use][0]:
					if cfg.mapConfig.file_transformation == 'None' or (cfg.mapConfig.file_transformation != 'None' and land_use in possible_elements[start_individual[i]]):
						new_cand[i] = land_use
		else:
			if cfg.mapConfig.file_transformation == 'None':
				option_elements = np.asarray(range(1,cfg.modelConfig.max_range+1))
				index = np.nonzero(option_elements == land_use)[0]
				option_elements = np.delete(option_elements, index)
			else:
				index = np.nonzero(possible_elements[land_use] == land_use)[0]
				option_elements = np.delete(possible_elements[land_use], index)
			for i in start_land_cover[land_use][0]:
				if len(option_elements) != 0:
					clbrValue = random.choice(option_elements.tolist())
					new_cand[i] = clbrValue
				else:
					# no other land use option possible -> no extreme seed possible for this land use
					break

	else:				
		WriteLogMsg("create individual with create_extreme_seed based on total area constraints and transition matrix if given")  

		cand_accept = False
		count_loops = 0
		
		# determine all possible land_use options
		if cfg.mapConfig.file_transformation == 'None':
			option_positions = range(len(start_individual))
		else:
			option_positions = []
			for key in possible_elements:
				if land_use in possible_elements[key]:
					index = np.nonzero(np.asarray(start_individual) == key)[0]
					for i in index:
						option_positions.append(i)
						
		candidate_list.append("option_positions: %s" %option_positions)
						
		# check if modifications are possible
		if (maximize == True and len(option_positions) > len(start_land_cover[land_use][0])) or \
			(maximize == False and cfg.mapConfig.file_transformation == 'None') or \
			(maximize == False and cfg.mapConfig.file_transformation != 'None' and len(possible_elements[land_use]) > 1):
			
			candidate_list.append("modification possible")
		
			# memory for checked individuals with land use cover 
			memory = {}
			# indices of option_positions which have been checked before and were rejected
			rejected_indices = []
			# count starts if start_individual
			count_new_start = 1

			# check if no individual was found and no other one exists (start individual is excluded) 
			while cand_accept == False and len(memory) < pow(2,len(option_positions))-1:
				
				candidate_list.append("while loop started; new_cand %s" %new_cand)				
				
				count_loops += 1
				
				sum_old = 0
				sum_new = sum_old
				# calculate land use cover of new_cand
				for m in np.nonzero(np.asarray(new_cand) == land_use)[0]:
					# index m+1 because of the exclusion of patch ID 0
					sum_old += map_proportion[m+1] 
					
				# check if land_use is in all elements/no element of option_positions 
				# then no more/less land use cover is possible
				all_land_use = True
				if maximize == True:
					for i in option_positions:
						if new_cand[i] != land_use:
							all_land_use = False
							break
				else:
					for i in option_positions:
						if new_cand[i] == land_use:
							all_land_use = False
							break
				if all_land_use == True:
					# found candidate
					cand_accept = True
					candidate_list.append("all option_position with/without land_use")
					candidate_list.append("found feasible solution with min/max land use: %s, %s" %(new_cand, sum_new))
					break
				
				# indices of option_positions that have been changed before
				changed_indices = []				
				
				# determine changed_indices 
				for i in range(len(new_cand)):
					if new_cand[i] != start_individual[i]:
						changed_indices.append(i)
						
				# mark the checked indices with ones and all other with zeros
				modifications = np.zeros(len(start_individual), dtype = bool)
				for i in changed_indices:
					modifications[i] = 1
					
				# add new_cand to memory
				if tuple(modifications) not in memory:
					memory.update({tuple(modifications): [new_cand, sum_old]})

				if count_loops % 10000 == 0:
					WriteLogMsg("count_loops: %s" %count_loops)
					WriteLogMsg("count_new_start: %s" %count_new_start)
					WriteLogMsg("new_cand: %s" %new_cand)
					candidate_list.append("%s" %['preliminary result'])
					WriteCandidateList(candidate_list)
					
				candidate_list.append("maximize %s, sum_old %s, min %s, max %s" %(maximize,sum_old,min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]], min_max_diff[2][np.nonzero(min_max_diff[0] == land_use)[0][0]])) 
				
				# select one index which could be important for minimization/maximization of land use cover
				if (maximize == True and sum_old < min_max_diff[2][np.nonzero(min_max_diff[0] == land_use)[0][0]]) or \
					(maximize == False and sum_old < min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]):
					# exclude indices which have been changed before 
					# (only for increasing -> return to the last checked new_cand must be possible)
					excluded_by_memory = [] 
					for element in memory:
						element_ones = np.nonzero(np.asarray(element) == True)[0]
						# check only memory elements with changed_indices length +1
						if len(element_ones) == len(changed_indices)+1:
							new_index = [i for i in element_ones if i not in changed_indices]
							# check only memory elements with difference in one element
							if len(new_index) == 1:
								excluded_by_memory.append(new_index[0])
						
					possible_indices = [i for i in option_positions if i not in rejected_indices and i not in np.nonzero(np.asarray(new_cand) == land_use)[0] and i not in changed_indices and i not in excluded_by_memory]
					candidate_list.append("option_positions: %s, rejected_indices %s, indices %s of new_cand with land_use %s" %(option_positions, rejected_indices, np.nonzero(np.asarray(new_cand) == land_use)[0], land_use))
					candidate_list.append("changed_indices %s, excluded_by_memory %s" %(changed_indices,excluded_by_memory))					
					candidate_list.append("possible_indices for increasing land cover: %s" %possible_indices)

						
				else:
					possible_indices = [i for i in option_positions if i not in rejected_indices and i in np.nonzero(np.asarray(new_cand) == land_use)[0]]
					candidate_list.append("option_positions: %s, rejected_indices %s, indices %s of new_cand with land_use %s" %(option_positions, rejected_indices, np.nonzero(np.asarray(new_cand) == land_use)[0], land_use))
					candidate_list.append("possible_indices for reduction of land cover: %s" %possible_indices)
					
				candidate_list.append("possible_indices %s" %possible_indices)					 
					
				if len(possible_indices) > 0:
					check_index = random.choice(possible_indices)
					candidate_list.append("check_index: %s" %check_index)
				else:
					# try to optimize the current new_cand for minimization
					if maximize == False and sum_old > min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]:
						diff_cover = sum_old-min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]
						fixed_indices = []
						len_fix_indices = len(fixed_indices)
						# determine all indices which can have the land use option at hand
						indices_land_use = np.nonzero(np.asarray(new_cand) == land_use)[0]
						indices_other_land_use = []
						for i in range(len(new_cand)):
							if i not in indices_land_use:
								indices_other_land_use.append(i)
						indices_other_land_use = np.asarray(indices_other_land_use)

						# look for a patch with map_proportion = min
						for i in range(len(new_cand)):
							if map_proportion[i+1] == min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]:
								# determine other possible land use options
								for k in indices_land_use:
									option_land_use = []
									for check_land_use in nonstatic_elements:
										if check_land_use != land_use:
											indices_check_land_use = np.nonzero(np.asarray(new_cand) == check_land_use)[0]
											sum_check_land_use = 0
											for j in indices_check_land_use:
												sum_check_land_use += map_proportion[j+1] 
											if min_max_diff[2][np.nonzero(min_max_diff[0] == check_land_use)[0][0]] - sum_check_land_use >= map_proportion[k+1]:
												option_land_use.append(check_land_use)
									# and assign those other land use options
									new_cand[k] = random.choice(option_land_use)
								new_cand[i] = land_use
								fixed_indices.append(i)
								break
						# if no patch was found then look for one patch with another land use 
						# so that min is reached by changing this patch
						if len_fix_indices == len(fixed_indices):
							for i in indices_land_use:
								if len_fix_indices != len(fixed_indices):
									break
								if map_proportion[i+1] == diff_cover:
									# determine other possible land use options
									option_land_use = []
									for check_land_use in nonstatic_elements:
										if check_land_use != land_use:
											indices_check_land_use = np.nonzero(np.asarray(new_cand) == check_land_use)[0]
											sum_check_land_use = 0
											for j in indices_check_land_use:
												sum_check_land_use += map_proportion[j+1] 
											if min_max_diff[2][np.nonzero(min_max_diff[0] == check_land_use)[0][0]] - sum_check_land_use >= map_proportion[i+1]:
												option_land_use.append(check_land_use)
									# and assign those other land use options
									new_cand[i] = random.choice(option_land_use)
									break
								elif map_proportion[i+1] > diff_cover:
									for j in indices_other_land_use:
										if map_proportion[j+1] == map_proportion[i+1]-diff_cover:
											# determine other possible land use options
											option_land_use = []
											for check_land_use in nonstatic_elements:
												if check_land_use != land_use:
													indices_check_land_use = np.nonzero(np.asarray(new_cand) == check_land_use)[0]
													sum_check_land_use = 0
													for k in indices_check_land_use:
														sum_check_land_use += map_proportion[k+1] 
													if min_max_diff[2][np.nonzero(min_max_diff[0] == check_land_use)[0][0]] - sum_check_land_use >= map_proportion[i+1]:
														option_land_use.append(check_land_use)
											# and assign those other land use options
											new_cand[i] = random.choice(option_land_use)
											new_cand[j] = land_use
											fixed_indices.append(j)
											break
						# if no single patch could be found then look for two patches 
						# to reach the minimum area
						if len_fix_indices == len(fixed_indices):
							# try to replace two smaller land use patches to reach min
							smaller_patches = []
							indices_smaller_patches = []
							for i in range(len(new_cand)):
								if len_fix_indices == len(fixed_indices):
									if map_proportion[i+1] < min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]:
										if len(smaller_patches)==0:
											smaller_patches.append(map_proportion[i+1])
											indices_smaller_patches.append(i)
										else:
											for index in range(len(smaller_patches)):
												if smaller_patches[index] + map_proportion[i+1] == min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]:
													# determine other possible land use options
													for k in indices_land_use:
														option_land_use = []
														for check_land_use in nonstatic_elements:
															if check_land_use != land_use:
																indices_check_land_use = np.nonzero(np.asarray(new_cand) == check_land_use)[0]
																sum_check_land_use = 0
																for j in indices_check_land_use:
																	sum_check_land_use += map_proportion[j+1] 
																if min_max_diff[2][np.nonzero(min_max_diff[0] == check_land_use)[0][0]] - sum_check_land_use >= map_proportion[k+1]:
																	option_land_use.append(check_land_use)
														# and assign those other land use options
														new_cand[k] = random.choice(option_land_use)
													new_cand[i] = land_use
													new_cand[index] = land_use
													fixed_indices.append(i)
													fixed_indices.append(index)
													break
											if len_fix_indices == len(fixed_indices):
												smaller_patches.append(map_proportion[i+1])
												indices_smaller_patches.append(i)
												
						# check if min has been reached
						sum_new = 0
						for m in np.nonzero(np.asarray(new_cand) == land_use)[0]:
							# index m+1 because of the exclusion of patch ID 0
							sum_new += map_proportion[m+1] 
						if sum_new == min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]] and individual_filter(new_cand) == True:
							# found candidate
							cand_accept = True
							candidate_list.append("try to optimize last new_cand by replacing: found feasible solution with min land use: %s, %s" %(new_cand, sum_new))
							continue
						else:
							candidate_list.append("try to optimize last new_cand: maybe an other optimize algorithm could be implemented for this case 0! new_cand %s, sum_new %s" %(new_cand, sum_new))
					
					# try to optimize the current new_cand for maximization
					if maximize == True and sum_old < min_max_diff[2][np.nonzero(min_max_diff[0] == land_use)[0][0]]:
						# determine min possible land cover
						sum_max = 100 - map_proportion[0]
						not_min_land_use = []
						fixed_indices = []
						for check_land_use in nonstatic_elements:
							if check_land_use != land_use:
								if min_max_diff[1][np.nonzero(min_max_diff[0] == check_land_use)[0][0]] > 0:
									sum_max -= min_max_diff[1][np.nonzero(min_max_diff[0] == check_land_use)[0][0]]
									# determine land use options that could be optimized 
									# from check_land_use in land_use and vice versa
									sum_current_cover = 0
									for m in np.nonzero(np.asarray(new_cand) == check_land_use)[0]:
										# index m+1 because of the exclusion of patch ID 0
										sum_current_cover += map_proportion[m+1] 
									if sum_current_cover > min_max_diff[1][np.nonzero(min_max_diff[0] == check_land_use)[0][0]]:
										not_min_land_use.append([check_land_use,sum_current_cover-min_max_diff[1][np.nonzero(min_max_diff[0] == check_land_use)[0][0]]])
						# if sum_old is min possible land cover then candidate is optimal
						if sum_old == sum_max and individual_filter(new_cand) == True:
							# found candidate
							cand_accept = True
							candidate_list.append("try to optimize last new_cand: found feasible solution with max/min land use: %s, %s" %(new_cand, sum_old))
							continue
						elif sum_old < sum_max and len(not_min_land_use) > 0:
							for element in not_min_land_use:
								check_land_use = element[0]
								diff_cover = element[1]
								check_option_positions = []
								# determine all indices which can have the land use option at hand
								indices_land_use = np.nonzero(np.asarray(new_cand) == land_use)[0]
								indices_check_land_use = np.nonzero(np.asarray(new_cand) == check_land_use)[0]
								for i in indices_land_use:
									check_option_positions.append(i)
								for i in indices_check_land_use:
									check_option_positions.append(i)
								len_fix_indices = len(fixed_indices)
								# among those indices look for a patch with map_proportion = min
								for i in check_option_positions:
									if map_proportion[i+1] == min_max_diff[1][np.nonzero(min_max_diff[0] == check_land_use)[0][0]]:
										# use this one
										for j in indices_check_land_use:
											new_cand[j] = land_use
										new_cand[i] = check_land_use
										fixed_indices.append(i)
										break
								# if no patch was found then look for a patch 
								# which could be changed that land cover is min
								if len_fix_indices == len(fixed_indices):
									# look for patches equal to or larger than diff_cover and try to replace it
									for i in indices_check_land_use:
										if len_fix_indices != len(fixed_indices):
											break
										if map_proportion[i+1] == diff_cover:
											new_cand[i] = land_use 
											break
										elif map_proportion[i+1] > diff_cover:
											for j in indices_land_use:
												if map_proportion[j+1] == map_proportion[i+1]-diff_cover:
													new_cand[i] = land_use
													new_cand[j] = check_land_use
													fixed_indices.append(j)
													break
											# try to replace two smaller land use patches to reach min
											smaller_patches = []
											indices_smaller_patches = []
											for j in indices_land_use:
												if len_fix_indices == len(fixed_indices):
													if map_proportion[j+1] < map_proportion[i+1]-diff_cover:
														if len(smaller_patches)==0:
															smaller_patches.append(map_proportion[j+1])
															indices_smaller_patches.append(j)
														else:
															for index in range(len(smaller_patches)):
																if smaller_patches[index] + map_proportion[j+1] == map_proportion[i+1]-diff_cover:
																	new_cand[i] = land_use
																	new_cand[j] = check_land_use
																	new_cand[index] = check_land_use
																	fixed_indices.append(j)
																	fixed_indices.append(index)
																	break
															if len_fix_indices == len(fixed_indices):
																smaller_patches.append(map_proportion[j+1])
																indices_smaller_patches.append(j)
								candidate_list.append("try to optimize: new_cand %s, changed land_use %s" %(new_cand, element))
							# check if min is reached
							sum_new = 0
							for m in np.nonzero(np.asarray(new_cand) == land_use)[0]:
								# index m+1 because of the exclusion of patch ID 0
								sum_new += map_proportion[m+1] 
							if sum_new == sum_max and individual_filter(new_cand) == True:
								# found candidate
								cand_accept = True
								candidate_list.append("try to optimize last new_cand by replacing: found feasible solution with max/min land use: %s, %s" %(new_cand, sum_new))
								continue
							else:
								candidate_list.append("try to optimize last new_cand: maybe an other optimize algorithm could be implemented for this case 1! new_cand %s, sum_new %s" %(new_cand, sum_new))
											
						else:
							candidate_list.append("try to optimize last new_cand: maybe an other optimize algorithm could be implemented for this case 2!")
					
					# if no optimization of new_cand was successfully
					if cand_accept != True:
						# for this new_cand all changes for increasing are checked before
						# continue with start_individual
						new_cand = start_individual[:] 
						# delete rejected indices
						rejected_indices = []
						# count for test the start with start_individual
						count_new_start += 1
						candidate_list.append("no possible indices for new_cand! continue with start_individual, count new start: %s" %count_new_start)
						# start next loop
						continue
					
				# determine possible land use classes for selected index
				if cfg.mapConfig.file_transformation != 'None':
					possible_land_use = possible_elements[start_individual[check_index]]			
					option_elements = determine_possible_landuse(possible_land_use, check_index, new_cand)
				else:
					option_elements = determine_possible_landuse(np.arange(1,cfg.modelConfig.max_range+1), check_index, new_cand)
				# delete checked land use from possible land use classes if land cover should be reduced
				if (maximize == False and sum_old > min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]) or \
					(maximize == True and sum_old > min_max_diff[2][np.nonzero(min_max_diff[0] == land_use)[0][0]]):
					if land_use in option_elements: 
						index = np.nonzero(option_elements == land_use)[0]
						option_elements = np.delete(option_elements, index)
				
				# if option_elements is empty then continue with other index
				if len(option_elements) == 0:
					modifications[check_index] = 1
					# add checked candidate to memory with worst land use cover
					if maximize == True:
						memory.update({tuple(modifications): [new_cand, -1]})
					else:
						memory.update({tuple(modifications): [new_cand, 110]})
					# add index to rejected indices
					rejected_indices.append(check_index)
					#candidate_list.append("no option_elements for checked_index! modifications %s, memory:" %modifications)
					candidate_list.append("no option_elements for checked_index! modifications %s:" %modifications)
					#for element in memory:
						#candidate_list.append("%s" %(element,))
					# start next loop
					continue
					
				changed_new_cand = False
				old_land_use = new_cand[check_index]
				#if maximize == False:
				# for reduction of land cover change checked land use in an other one
				if (maximize == False and sum_old > min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]) or \
					(maximize == True and sum_old > min_max_diff[2][np.nonzero(min_max_diff[0] == land_use)[0][0]]):
					# check if one of the new candidates are feasible and improve the land use cover
					for check_landuse in option_elements:
						new_cand[check_index] = check_landuse
						if individual_filter(new_cand) == True:
							changed_new_cand = True
							sum_new = sum_old - map_proportion[check_index+1]
							if (maximize == False and sum_new == min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]) or \
								(maximize == True and sum_new == min_max_diff[2][np.nonzero(min_max_diff[0] == land_use)[0][0]]):
								# found candidate
								cand_accept = True
								candidate_list.append("found feasible solution with min/max land use: %s, %s" %(new_cand, sum_new))
							break
						else:
							#add infeasible candidate to impossible candidates
							check_redundancy(new_cand)
				# for increasing of land cover change land use of check_index in the needed land_use
				else:
					new_cand[check_index] = land_use  
					if individual_filter(new_cand) == True:
						changed_new_cand = True
						sum_new = sum_old + map_proportion[check_index+1]
						if (maximize == True and sum_new == min_max_diff[2][np.nonzero(min_max_diff[0] == land_use)[0][0]]) or \
							(maximize == False and sum_new == min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]):
							# found candidate
							cand_accept = True
							candidate_list.append("found feasible solution with max/min land use: %s, %s" %(new_cand, sum_new))
					else:
						#add infeasible candidate to impossible candidates
						check_redundancy(new_cand)
						
				# no land use generates feasible individual
				if changed_new_cand == False:
					# reverse last change of new_cand
					new_cand[check_index] = old_land_use
					# add checked candidate to memory with worst land use cover
					# if it was not checked before
					if modifications[check_index] == 1:
						modifications[check_index] = 0
					else:
						modifications[check_index] = 1
					if tuple(modifications) not in memory:
						if maximize == True:
							memory.update({tuple(modifications): [new_cand, -1]})
						else:
							memory.update({tuple(modifications): [new_cand, 110]})
					# add index to rejected indices
					rejected_indices.append(check_index)
					candidate_list.append("no feasible new cand for checked_index! modifications %s, rejected_indices %s" %(modifications, rejected_indices))
					#candidate_list.append("memory:")
					#for element in memory:
					#	candidate_list.append("%s"%(element,))					
					# start next loop
					continue
					
				else: 
					# add checked candidate to memory 
					# if it was not checked before
					if modifications[check_index] == 1:
						modifications[check_index] = 0
					else:
						modifications[check_index] = 1
					if tuple(modifications) not in memory:
						# add checked candidate to memory with land use cover
						memory.update({tuple(modifications): [new_cand, sum_new]})
					# start next loop
					continue

			# end of loop
			# current new_cand have min/max land cover range
			if cand_accept == True:
				candidate_list.append("%r individuals are checked." % count_loops)
			# no checked new_cand match the min/max land cover range
			else:
				new_cand = start_individual[:]
				if maximize == True:
					best_land_cover = 0
					for element in memory:
						if memory[element][1] > best_land_cover:
							new_cand = memory[element][0]
							best_land_cover = memory[element][1]
				else:
					best_land_cover = 100
					for element in memory:
						if memory[element][1] < best_land_cover:
							new_cand = memory[element][0]
							best_land_cover = memory[element][1]
				candidate_list.append("best new_cand %s with land cover %s from memory:" % (new_cand, best_land_cover))
				for element in memory:
					WriteLogMsg(element)	  

	# end
	end = time.time()

	# if generated individual is infeasible or excluded by the tabu memory
	# then return an empty list
	if individual_filter(new_cand) == False or check_impossible_candidates(new_cand) == True:
		WriteLogMsg("create_extreme_seeds: No new extreme cover for land_use %s and maximize %s was found." % (land_use, maximize))
		candidate_list.append("%s rejected: infeasible %s or excluded by impossible canddiate list %s" %(new_cand,individual_filter(new_cand),check_impossible_candidates(new_cand)))
		new_cand = []	
	else:
		# add accepted new candidate to impossible candidates 
		# to avoid the repeatedly usage of the accepted candidate as individual
		check_redundancy(new_cand)
		WriteLogMsg("The generation of a new candidate needed %d seconds." %(end-begin))
		#WriteLogMsg("%r" % new_cand)

	candidate_list.append("%s end create_extreme_seeds/n" %new_cand)
		
	# line for printing the checked new_cand and information of the generation process
	if (end-begin >= 120 and cfg.ea.write_tabu_memory == True):
		WriteCandidateList(candidate_list)

	return new_cand
		  
#------------------------------------------------------------------------------	 
#	Generate first parameter for optimization algorithm
#------------------------------------------------------------------------------
def generate_parameter( random, args, custom_individual):	
	"""Generates first set of candidates for algorithm. 
		Return one population of the first set per call of this function.
	"""

	# first_ind is True for the first call of the generate_parameter function
	global first_ind
	# original start individual of the input data
	global start_individual
	# Array for min/max proportional difference of land use type
	global min_max_diff
	# Array for static land use types
	global static_elements
	# Array for non static land use types
	global nonstatic_elements
	# Array for the area percentage of the patches
	global map_proportion
	# Dictionary with patch indices and total area per land use category
	global start_land_cover
	# Dictionary with information per land use category if extreme seeds were created before
	global extreme_seeds_dict
	
	global custom_index


	new_cand = []

	# Use custom individuals first
	if cfg.ea.start_from_previous_gen == True and custom_index < len(custom_individual):
		print(len(custom_individual))
		WriteLogMsg(f"Using custom individual #{custom_index + 1}")
		new_cand.extend(custom_individual[custom_index])
		custom_index += 1
		first_ind = False
		print(f"Starting  with the individual: {new_cand}")
		return new_cand


	# analyse land cover of start individual for extreme seeds 
	if cfg.ea.extreme_seeds == True and first_ind != True and len(start_land_cover) == 0:
		for i in range(1,cfg.modelConfig.max_range+1):
			if i in nonstatic_elements:
				sum_new = 0
				for m in np.nonzero(np.array(start_individual)==i)[0]:
					# index m+1 because of the exclusion of patch ID 0 of the optimization
					sum_new += map_proportion[m+1]
				start_land_cover.update({i: [np.nonzero(np.array(start_individual)==i)[0],sum_new]})

		WriteLogMsg("start_land_cover: %s" % start_land_cover)
		
	
	# read information of min/max proportional deviation of land use types
	if cfg.mapConfig.file_difference != 'None' and len(min_max_diff)==0:
		WriteLogMsg("Read min/max limits ...")
		min_max_diff = np.genfromtxt(os.path.join(wrkDir, 'input', cfg.mapConfig.file_difference), dtype=int, filling_values = '0')[:,1:]
		if np.amin(min_max_diff) < 0:
			msg = "Error: One or more min/max proportional elements are negativ (%s). Please check the input file." %np.amin(min_max_diff)
			WriteLogMsg(msg) 
			raise SystemError("Error: One or more min/max proportional elements are negativ (%s). Please check the input file." %np.amin(min_max_diff))
			req.close_window
		WriteLogMsg("min_max_diff: %s" %min_max_diff)

	# First population in the first generation is the start individual of the input data
	if first_ind == True:	
		WriteLogMsg("Generator candidates are: ")
		for i in range(0, len(start_individual)):
			# Add the original population elements					
			new_cand.append(start_individual[i]) 
		# for logical and filter_mutation to avoid that the start individual is repeatedly used as individual
		check_redundancy(new_cand)
		# check if start individual match min/max proportional
		if cfg.mapConfig.file_difference != 'None':
		# check min/max rules for land use types are included in start_individual
			for i in np.unique(new_cand):
				sum_new = 0
				#sum_old = 0
				for m in np.nonzero(new_cand==i)[0]:
					# index m+1 because of the exclusion of patch ID 0 of the optimization
					sum_new += map_proportion[m+1]
				if not (min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]]-.00001 <= sum_new <= min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]+.00001):
					msg = "Error: Land use ID %s covers %s percent of the map and doesn`t match the min/max rules. Please check the input data." %(i, sum_new)
					WriteLogMsg(msg) 
					raise SystemError("Error: Land use ID %s covers %s percent of the map and doesn`t match the min/max rules. Please check the input data." %(i, sum_new))
					req.close_window
			# check min/max rules for land use types are not included in start_individual -> min should be zero
			for i in range(1,cfg.modelConfig.max_range+1):
				if (i not in static_elements) and (i not in np.unique(new_cand)):
					if min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] > 0:
						msg = "Error: Land use ID %s covers 0 percent of the map and doesn`t match the min/max rules. Please check the input data." %i
						WriteLogMsg(msg) 
						raise SystemError("Error: Land use ID %s covers 0 percent of the map and doesn`t match the min/max rules. Please check the input data." %i)
						req.close_window

	# Create extreme seeds of land cover per land use category
	elif cfg.ea.extreme_seeds == True and first_ind != True and (len(extreme_seeds_dict) < len(start_land_cover) or (len(extreme_seeds_dict) == len(start_land_cover) and extreme_seeds_dict[max(extreme_seeds_dict.keys())][1] == False)) and cfg.modelConfig.max_range > 1:

		# if extreme seed is generated or all extreme land covers are generated/checked
		# then break the loop
		while len(new_cand) == 0 and (len(extreme_seeds_dict) < len(start_land_cover) or (len(extreme_seeds_dict) == len(start_land_cover) and extreme_seeds_dict[max(extreme_seeds_dict.keys())][1] == False)):
			land_use = 0

			#WriteLogMsg("extreme_seeds_dict: %s" % extreme_seeds_dict)

			# determine next land use category and direction of the extreme value
			if	len(extreme_seeds_dict) != 0 and extreme_seeds_dict[max(extreme_seeds_dict.keys())][1] == False:
				land_use = max(extreme_seeds_dict.keys())
				maximize = True
			elif len(extreme_seeds_dict) != len(start_land_cover) and len(start_land_cover) != 0:
				for key in start_land_cover.keys():
					if key not in extreme_seeds_dict:
						land_use = key
						break
				maximize = False

			# continue if one more extreme seed could exist
			if land_use > 0:
				# check if current land cover is one of the extreme individuals
				if cfg.mapConfig.file_difference == 'None':
					if maximize == False:
						# start individual is one map with minimal land cover for checked land use category
						# then update extreme_seeds_dict else look for a individual with extreme land cover
						if start_land_cover[land_use][1] != 0:
							# if no new individual with extreme land cover exists
							# then update extreme_seeds_dict and continue with next land_use/maximize combination
							new_cand = create_extreme_seed(land_use, maximize)
						extreme_seeds_dict.update({land_use:[True, False]})

					else:
						# start individual is one map with maximal land cover for checked land use category 
						# then update extreme_seeds_dict else look for a individual with extreme land cover
						if start_land_cover[land_use][1] != 100-map_proportion[0]:
							new_cand = create_extreme_seed(land_use, maximize)
						extreme_seeds_dict.update({land_use:[True, True]})
				else:
					if maximize == False:
						# start individual is one map with minimal land cover for checked land use category
						# then update extreme_seeds_dict else look for a individual with extreme land cover
						if start_land_cover[land_use][1] != min_max_diff[1][np.nonzero(min_max_diff[0] == land_use)[0][0]]:
							# if no new individual with extreme land cover exists
							# then update extreme_seeds_dict and continue with next land_use/maximize combination
							new_cand = create_extreme_seed(land_use, maximize)
						extreme_seeds_dict.update({land_use:[True, False]})

					else:
						# start individual is one map with maximal land cover for checked land use category 
						# then update extreme_seeds_dict else look for a individual with extreme land cover
						if start_land_cover[land_use][1] != min_max_diff[2][np.nonzero(min_max_diff[0] == land_use)[0][0]]:
							new_cand = create_extreme_seed(land_use, maximize)
						extreme_seeds_dict.update({land_use:[True, True]})
		
	# if no new_cand was generated before
	if len(new_cand) == 0:
		# Create individuals which takes account of the transformation information and the min/max rules
		if first_ind != True and (cfg.mapConfig.file_transformation != 'None' or  cfg.mapConfig.file_difference != 'None'):
			# call logical_variator if repair_mutation is selected as variator 
			# or a feasible first population is desired 
			if 'repair_mutation' in cfg.ea.variator or cfg.ea.feasible_first_pop == 'True':
				new_cand = logical_variator(start_individual,'True')

			# if no repair_mutation but filter_mutation than the first set is generated by the filter function
			elif 'filter_mutation' in cfg.ea.variator:
				new_cand = filter_variator(start_individual)
		
			else:
				for i in range(0, len(start_individual)):
					# generates random integers 1 <= N <= max_range
					clbrValue = random.randint(1, cfg.modelConfig.max_range)
					new_cand.append(clbrValue)
				WriteLogMsg("%r" % new_cand)
				# for logical and filter_mutation to avoid that the new candidate is repeatedly used as individual
				check_redundancy(new_cand)
		# Create individuals randomly
		else:
			for i in range(0, len(start_individual)):
				# generates random integers 1 <= N <= max_range
				clbrValue = random.randint(1, cfg.modelConfig.max_range)
				new_cand.append(clbrValue)
			WriteLogMsg("%r" % new_cand)
			# for logical and filter_mutation to avoid that the new candidate is repeatedly used as individual
			check_redundancy(new_cand)
	
	if first_ind == True:
		first_ind = False	

	return new_cand

#------------------------------------------------------------------------------	 
#	Starts the variator function with filter method
#------------------------------------------------------------------------------
def generate_cand_filter(random, candidate, args):	
	"""Starts the filter variator for generation of a new candidate. 
		Return one individual per call of this function.

		New mutator is implemented in the following files of the inspyred library:
			inspyred.ec.variators __init__.py
			inspyred.ec.mutators  mutators.py
		
		input data: 
			candidate is the given individual which is usually used for the 
						 creation of a new candidate, but in this special variator 
						 only feasible solutions are returned (rejection of infeasible solutions)
	
	"""
	#WriteLogMsg("Input parameter candidate in generate_cand_filter: %r " % candidate)
	new_cand = filter_variator(candidate)
	return new_cand

#------------------------------------------------------------------------------	 
#	Generate new candidate for optimization algorithm with logical function
#------------------------------------------------------------------------------
def generate_cand_logical(random, candidate, args): 
	"""Starts the logical variator for generation of a new candidate. 
		Return one individual per call of this function.

		New mutator is implemented in the following files of the inspyred library:
			inspyred.ec.variators __init__.py
			inspyred.ec.mutators  mutators.py
		
		input data: 
			candidate is the given individual which is usually used for the 
						 creation of a new candidate, but this special variator 
						 repairs the candidate at the elements where constraints are violated
	"""
	
	#WriteLogMsg("Input parameter candidate in generate_cand_filter: %r " % candidate)
	new_cand = logical_variator(candidate)
	return new_cand

#----------------------------------------------------------------------------------	 
#	Determine the informations for penalty function of constraint_tourn_selection
#----------------------------------------------------------------------------------
def analyse_violation(candidate):
	"""Analyse where the individual violate the constraints and collect the 
		information for the penalty function of the constraint tournament selection.
		Return an array with
			[0] Number of land use categories which violate the cover ranges
			[1] Sum of differences between real area and cover ranges
				for all land use categories which violate the cover ranges
			[2] Number of non static land use categories
			[3] Number of patches which violate the constraints of the 
				transformation information 
			[4] Sum of the area of all patches which violate the 
				constraints of the transformation information
		
		input data: 
			candidate is an individual for checking constraint violation
	"""
	
	# Array for transformation information
	global trans_matrix
	# Array for the area percentage of the patches
	global map_proportion
	# Array for min/max proportional difference of land use type
	global min_max_diff
	# Array for the start individual
	global start_individual
	# Array for static land use types
	global static_elements
	# Dictionary for area percentage of the static land use categories
	global static_area 
	# Array for original ASCII map
	global map

	# if static_area is empty than fill it with the area of the land use in the original map
	if len(static_area) == 0:
		all_cells = map.shape[0]*map.shape[1]
		for row in range(0,map.shape[0]):
			for col in range(0,map.shape[1]):
				if map[row][col] in static_elements:
					if map[row][col] in static_area:
						static_area[map[row][col]] += 1/all_cells*100
					else:
						static_area.update({map[row][col]: 1/all_cells*100})
	
	compare_individual = start_individual

	#WriteLogMsg("candidate in analyse_violation: %s " % candidate)
	
	number_landuse = 0
	sum_landuse = 0
	number_patches = 0
	sum_patches = 0

	# check if transformation of the individual elements are allowed
	if cfg.mapConfig.file_transformation != 'None':				
		for i in range(0,len(candidate)):
			if trans_matrix[np.nonzero(np.unique(trans_matrix[:,:1]) == compare_individual[i])[0][0]][np.nonzero(trans_matrix[0] == candidate[i])[0][0]] != 1:
				number_patches += 1
				# index i+1 because of the exclusion of patch ID 0 of the optimization
				sum_patches += map_proportion[i+1]
	
	# check if minimal and maximal difference of the proportion of land use classes are plausibel
	if cfg.mapConfig.file_difference != 'None':
		# land use exists in candidate
		for i in np.unique(candidate):
			sum_new = 0		
			for m in np.nonzero(candidate==i)[0]:
				# index m+1 because of the exclusion of patch ID 0 of the optimization
				sum_new += map_proportion[m+1]
			if i in static_area:
				sum_new += static_area[i]
			# for plausibility: min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] <= sum_new <= min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]
			if sum_new < min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]]:
				number_landuse += 1
				sum_landuse += min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] - sum_new 
			elif sum_new > min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]:
				number_landuse += 1
				sum_landuse += sum_new - min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]

		# check min/max rules for land use types are not included in candidate 
		# for static land use categories -> check min and max
		for i in range(1,cfg.modelConfig.max_range+1):
			# for dynamic land use categories -> min should be zero
			if (i not in static_elements) and (i not in np.unique(candidate)):
				if min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] > 0:
					number_landuse += 1
					sum_landuse += min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]]
			# for static elements with area > 0
			elif (i in static_area) and (i not in np.unique(candidate)):
				if static_area[i] < min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]]:
					number_landuse += 1
					sum_landuse += min_max_diff[1][np.nonzero(min_max_diff[0] == i)[0][0]] - static_area[i] 
				elif static_area[i] > min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]:
					number_landuse += 1
					sum_landuse += static_area[i] - min_max_diff[2][np.nonzero(min_max_diff[0] == i)[0][0]]

	violation_info = [number_landuse,sum_landuse,number_patches,sum_patches]
	#WriteLogMsg("violation_info[number_landuse,sum_landuse,number_patches,sum_patches]: %s " % violation_info)
	
	return violation_info

#-------------------------------------------------------------------------------------	
#	Generate a tournament sampling of individuals in consideration of the constraints
#-------------------------------------------------------------------------------------
def constraint_tourn_selection(random, population, args):
	"""Return a tournament sampling of individuals in consideration of the constraint 
		violation from the population.
	
		This function selects ``num_selected`` individuals from the population. 
		It selects each one by using random sampling without replacement
		to pull ``tournament_size`` individuals and adds the best of the
		tournament as its selection. If ``tournament_size`` is greater than
		the population size, the population size is used instead as the size
		of the tournament.

		New selection is an extension of the tournament_selection and is implemented 
		in the following files of the inspyred library:
			inspyred.ec selectors.py
	
		input data:
			random -- the random number generator object
			population -- the population of individuals
			args -- a dictionary of keyword arguments

		Optional keyword arguments in args:
	
			- *num_selected* -- the number of individuals to be selected (default 1)
			- *tournament_size* -- the tournament size (default 2)
	
	"""

	# Dictionary as memory for the analysis results of the candidate violation
	global violation_memory

	num_selected = cfg.ea.num_selected
	tournament_size = cfg.ea.tournament_size

	# penalty_function != 0 for tournament_selection with constraint handling
	penalty_function = cfg.ea.penalty_function

	if tournament_size > len(population):
		tournament_size = len(population)
	selected = []
	for _ in range(num_selected):
		tourn = random.sample(population, tournament_size)

		# tournament_selection without penalty function
		if penalty_function == 0:
			selected.append(max(tourn))
		# tournament_selection with penalty function
		elif penalty_function == 1 or penalty_function == 2:
			# list with violation_values
			violation_value = []
			for element in tourn:

				# candidate not listed in violation_memory -> determine violation value
				if tuple(element.candidate) not in violation_memory:
					# violation_info[number_landuse,sum_landuse,number_patches,sum_patches]
					violation_info = analyse_violation(element.candidate)

					if penalty_function == 1:
						# violation value = add sum of area differences of land uses from cover ranges
						# to the sum of the areas of the patches with transformation violation
						value = violation_info[1] + violation_info[3]
					# penalty_function == 2
					else:
						# violation value = add quotient of the sum of area differences of land uses from cover ranges
						# and number of violating land use categories muliplied with 100 to the quotient of the number 
						# of the patches with transformation violation and individuen length (number of all patches)
						if violation_info[0] > 0:
							value1 = violation_info[1]/violation_info[0] 
						else:
							value1 = 0
						value = value1+violation_info[3]
					violation_value.append(value)
					# dictionary key must be immutable (tuple but no list is immutable)
					violation_memory.update({tuple(element.candidate): value})

				# candidate listed in violation_memory -> use violation value of memory
				else:
					violation_value.append(violation_memory[tuple(element.candidate)])					
					
			# if both feasible -> tournament selection
			if max(violation_value) == 0:
				selected.append(max(tourn))
			# if one or more are infeasible -> select candidate with minimum of violation
			else:
				selected.append(tourn[violation_value.index(min(violation_value))])
				#WriteLogMsg("min. one infeasible: violation_value %s, selected one %s " % (violation_value,tourn[violation_value.index(min(violation_value))]))
			
		# wrong penalty function setting
		else:
			selected.append(max(tourn))
			WriteLogMsg("For integer %s exist no penalty function. Please check your config.ini settings." % penalty_function)
			raise SystemError("Error: For integer %s exist no penalty function. Please check your config.ini settings." % penalty_function)
			req.close_window

	return selected

#-------------------------------------------------------------------------------------	
#	Transfer the variables for map creation to the optiAlgorithm.py
#-------------------------------------------------------------------------------------
def get_from_maphandler():
	"""Transfer the variables for map creation to the optiAlgorithm.py.
	
	"""
	
	# Array for original ASCII map
	global map
	# Array for the patch ID map
	global patchID_map
	# Array for the header information from ASCII map
	global header_all
	
	return map, patchID_map, header_all

#-------------------------------------------------------------------------------------	
#	Generate ascii-map from individual
#-------------------------------------------------------------------------------------
def transform_individual_ascii_map(individual, modelfolder=False, ind_number=0, map_info=None, patchID_map_info=None, header_all_info=None, feasible=True):
	"""Transform the genome back to an ascii file.
	
		input data: 
			individual as basic for the new ascii map
			modelfolder - if True than write the ascii-map in the modelfolder
						  else write it in the output-folder
			ind_number - describe in which model folder the map should be saved
						 or if best values are saved in the output folder than 
						 the number attributed the individual number to the map	 
			feasible - information for constraint_tournament_selection if individual is feasible
	"""
	
	#WriteLogMsg("Update ascii-map for %r" %individual)
	if map_info is None or patchID_map_info is None or header_all_info is None:
		# Array for original ASCII map
		global map
		# Array for the patch ID map
		global patchID_map
		# Array for the header information from ASCII map
		global header_all
	else:
		map = np.copy(map_info)
		patchID_map = np.copy(patchID_map_info)
		header_all = header_all_info
	
	# create an independent copy of the patchID_map for the new map
	ascii_map = np.copy(patchID_map)
	for row in range(ascii_map.shape[0]):
		for col in range(ascii_map.shape[1]):
			# replace static elements with original elements of the original map
			if ascii_map[row][col] == 0:
				ascii_map[row][col] = map[row][col]
			# replace non static elements with the land use typ of the individual
			else:
				index = ascii_map[row][col]-1
				ascii_map[row][col] = individual[index]
			  
	# write header information and ascii_map in new ascii-file
	WriteMap(header_all, ascii_map, modelfolder, ind_number, feasible)
  
#------------------------------------------------------------------------------	 
#	Test functions
#------------------------------------------------------------------------------ 
#if __name__ == "__main__":
	#cfg.mapConfig.file_HRU == "C:/SWAT_POF/Lossa/Optimization/hru_patchID_map"
	#max_range = 8
	#file_HRU = 'None'
	#file_HRU = 'Y:/Home/wegwitz/EclipseWorkspace/ToolCombiModels/ToolCombiModels/input/Size_HRUs_ff.csv'
	#map_file = 'None'
	#map_file = "Y:/Home/wegwitz/EclipseWorkspace/ToolCombiModels/ToolCombiModels/input/lu_saale_ganz_klein.asc"
	#trans_file = 'None'
	#trans_file = 'Y:/Home/wegwitz/EclipseWorkspace/ToolCombiModels/ToolCombiModels/input/transformation_info3.txt' 
	#file_difference = 'None'
	#file_difference = "Y:/Home/wegwitz/EclipseWorkspace/ToolCombiModels/ToolCombiModels/input/min_max2.txt"
	#patchIDmap_file = 'None'
	#patchIDmap_file = 'Y:/Home/wegwitz/EclipseWorkspace/ToolCombiModels/ToolCombiModels/input/test_lu_saale_klein_patch_ID_map_transinfo2_8neighbors.asc'
	#four_neighbours = 'False'
	#genom = generateASCIIgenom(max_range, map_file, trans_file, patchIDmap_file, four_neighbours)
	#print("patchID_map:\n%s" %patchID_map)
	#print("Genomlaenge: %s; genom: %s" %(len(genom),genom))
	
	#msg = "Read min/max proportional difference of the land use types ..."
	#WriteLogMsg(msg) 
	# write information in a matrix
	#min_max_diff = np.genfromtxt(file_difference, dtype=int, filling_values='0')
	#header_all = open("Y:/Home/wegwitz/EclipseWorkspace/ToolCombiModels/ToolCombiModels/input/lu_saale_ganz_klein.asc",'rb').readlines()[0:6]
	#generate_genom(max_range, file_HRU, map_file, trans_file, patchIDmap_file, four_neighbours)
	#possible_elements = np.array([1,3,4,6,8])
	#new_cand = [1,1,3]
	#impossible_cand = [[3,3,3],[1,1,3,8],[1,2,3],[1,1,3,4]]
	#print("determine_possible_elements: %s" %determine_possible_elements(possible_elements, new_cand, impossible_cand))
	#print("type, header_all: %s, %s" %(type(header_all), header_all))
	#global map_proportion
	#print ("map_proportion: %r" %map_proportion)
	#print ("sys.path: %r" %sys.path)
	#print(individual_filter([8, 4, 2, 6, 6, 1, 5]))
	
	#begin = time.time()
	#help_map = np.zeros([3000, 3000], int) # default is 0
	#for row in range(help_map.shape[0]):
	#	for col in range(help_map.shape[1]):
	#		help_map[row][col] = random.randint(1, cfg.modelConfig.max_range)
	#		
	#end = time.time()
	#print("The total runtime of the tool: %d seconds." %(end-begin))
#------------------------------------------------------------------------------
#
#	EOF
#
#------------------------------------------------------------------------------
