# __CoMOLA__

CoMOLA is a free Python tool to optimize the allocation of land use for multiple objectives. It builds upon the open source “inspyred” Python library and includes functions for reading, encoding and writing land use maps as well as genome generation and repair mutation algorithms for considering constraints during the optimization procedure. It runs on Windows and Linux and allows for the integration of any model whose prediction (e.g. a value for an ecosystem service) is based on a land use raster map. In its basic form, CoMOLA can be used immediately by inputting a raster map representing the status-quo land use, ready-to-run models written in R including their input data, and (optional) information on constraints. As constraints, the tool can consider (1) transition rules defining which type of land use can be converted into which other type and (2) minimum and maximum area proportions of each land use type within the study area. All relevant settings, such as paths to input data and models as well as optimization-specific parameters (e.g. population size, crossover and mutation rates) and settings related to constraint-handling and raster map-analysis are managed in one single control file ("config.ini").

## __Installation requirements__

CoMOLA was developed and tested for Python 2.7.

* Required Python packages
  * matplotlib
  * numpy
  * pylab

Furthermore you need to install R to run external models.

## __Input__
*(see example files in input folder)*

### __Maps__

__Land use map (required)__

Provide a land use raster map in ascii format (with consecutive integer values representing the different land use classes, starting with value 1). If no user-defined patch ID map is given, CoMOLA generates a patch map where neighboring raster cells of the same type are aggregated.

Example *(land\_use.asc)*: <pre>
ncols         10
nrows         10
xllcorner     4376461.4080843
yllcorner     5553063.2189224
cellsize      75
NODATA_value  -2
2 6 6 1 1 2 2 3 4 4
4 6 6 1 5 5 3 6 7 8
4 2 1 1 4 4 3 3 5 5
5 2 6 6 8 8 1 1 1 1
5 7 7 2 5 8 8 4 3 3
1 7 7 4 6 8 8 4 3 1
5 2 2 2 6 8 7 5 3 1
5 5 6 6 3 4 7 5 1 1
2 2 6 6 3 4 7 3 3 3
2 6 7 7 3 2 2 1 1 1</pre>

__Patch ID map (optional)__

If appropriate provide a Patch ID map as ascii file to delineate the spatial optimization units as needed in your specific case (with consecutive integer values representing the different patches, starting with value 1). For a cell-level optimization, an individual ID must be assigned to each cell. The spatial resolution must be the same as for the land use map. 

Example *(patch\_IDmap\_eachcell\_constraint.asc)*: <pre>
ncols         10
nrows         10
xllcorner     4376461.4080843
yllcorner     5553063.2189224
cellsize      75
NODATA_value  -2
1 2 3 4 5 6 7 8 9 10
11 12 13 14 15 16 17 18 19 0
20 21 22 23 24 25 26 27 28 29
30 31 32 33 0 0 34 35 36 37
38 39 40 41 42 0 0 43 44 45
46 47 48 49 50 0 0 51 52 53
54 55 56 57 58 0 59 60 61 62
63 64 65 66 67 68 69 70 71 72
73 74 75 76 77 78 79 80 81 82
83 84 85 86 87 88 89 90 91 92</pre>

### __Constraints__

CoMOLA can handle two types of land use change constraints (simultaneously or stand-alone):

__(1) Transition matrix__

Constraint defining which type of land use (given in rows) can be converted into which other type (given in columns). 1 = transition is possible, 0 = transition is not possible. 

Example *(transition\_matrix.txt)*: <pre>
-2 1 2 3 4 5 6 7 8
1 1 1 1 1 1 1 1 0
2 1 1 1 1 1 1 1 0
3 1 1 1 1 1 1 1 0
4 1 1 1 1 1 1 1 0
5 1 1 1 1 1 1 1 0
6 0 0 0 0 0 1 1 0
7 0 0 0 0 0 0 1 0
8 0 0 0 0 0 0 0 1</pre>


__(2) Total area__

Constraint defining minimum and maximum area proportions of each land use type within the study area.

Example *(min\_max.txt)*: <pre>
land_use 1 2 3 4 5 6 7 8
min 0 0 0 0 0 10 10 0
max 100 100 100 100 100 25 30 100</pre>

## __External models__

CoMOLA handles up to four external models which must be provided as R scripts and stored in a separate directory within the models folder. Each model evaluates the land use map (see above) for a specific objective that should be maximized (e.g. a certain ecosystem service). The output of each model is a single value representative for the whole study area (e.g. total agricultural yield) and needs to be written in a .csv file. If the objective value should be minimzed during optimization, multiply with -1.

Example *(SYM.R)*: <pre>
setwd("C:/+STRAUCH+/+PAPER\_WORK+/Opti-Tool/CoMOLA\_basic/models/SYM")
sink("C:/+STRAUCH+/+PAPER\_WORK+/Opti-Tool/CoMOLA\_basic/models/SYM/console.txt", append=FALSE)
##########################################################################################
#
#     ~ ~ ~ Simple Yield Model (SYM) ~ ~ ~
#     ~ ~ ~ this is just a toy model ~ ~ ~
#
#     ~ ~ ~ Input data ~ ~ ~
#    land_use.asc        |land use map containing the following classes
#                        |1,2,3,4,5 = arable land with increasing intensity from 1 to 5
#                        |6 = forest
#                        |7 = pasture
#                        |8 = urban area
#                        |-2 = no data
#
#    soil_fertility.asc  |map on soil fertility which can range from 0.1 to 1
#
#    Objective: Maximize crop yield
#
##########################################################################################

# set working directory

# read in ascii files
lu.map <- read.table("map.asc", h=F, skip=6, sep=" ")
fert.map <- read.table("soil\_fertility.asc", h=F, skip=6, sep=" ")

# array index for arable land
arable.idx <- which(lu.map <= 5 & lu.map > 0, arr.ind=T)

# calculate crop yield as logarithmic function of intensity and soil fertility
yield <- log(lu.map[arable.idx] * (1 + fert.map[arable.idx]))
yield[is.na(yield)] <- 0
yield.sum <- sum(yield)

# write model output
write.table(yield.sum , "SYM\_output.csv",append=FALSE ,sep =";",col.names=FALSE ,row.names=FALSE)

sink()</pre>

## __Configuration and optimization settings (config.ini)__

All relevant settings, such as paths to input data and models as well as optimization-specific parameters and settings related to constraint-handling and raster map-analysis are managed in one single control file called "config.ini".

## __Running CoMOLA__

Call from the console (within your CoMOLA folder): <pre>
python \_\_init\_\_.py </pre>

You can limit the maximum number of threads to be used for parallel computation by adding "-t x" to the command, where x is the maximum number of threads, e.g. <pre>
python \_\_init\_\_.py -t 2 </pre>