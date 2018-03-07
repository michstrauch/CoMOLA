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

__Land use map (required)__ *(see "land_use.asc")*

Provide a land use raster map in ascii format (with consecutive integer values representing the different land use classes, starting with value 1). If no user-defined patch ID map is given, CoMOLA generates a patch map where neighboring raster cells of the same type are aggregated.

__Patch ID map (optional)__ *(see "patch_IDmap_eachcell_constraint.asc")*

If appropriate provide a Patch ID map as ascii file to delineate the spatial optimization units as needed in your specific case (with consecutive integer values representing the different patches, starting with value 1). For a cell-level optimization, an individual ID must be assigned to each cell. The spatial resolution must be the same as for the land use map. 

### __Constraints__

CoMOLA can handle two types of land use change constraints (simultaneously or stand-alone):

__Transition matrix__ *(see "transition_matrix.txt")*

Constraint defining which type of land use (given in rows) can be converted into which other type (given in columns). 1 = transition is possible, 0 = transition is not possible. 

_Total area_ *(see "min_max.txt")*

Constraint defining minimum and maximum area proportions of each land use type within the study area.

## __External models__

CoMOLA handles up to four external models which must be provided as R scripts. Each model evaluates the land use map (see above) for a specific objective that should be maximized (e.g. a certain ecosystem service). The output of each model is a single value representative for the whole study area (e.g. total agricultural yield) and needs to be written in a .csv file. If the objective value should be minimzed during optimization, multiply with -1.

## __Configuration and optimization settings (config.ini)__

All relevant settings, such as paths to input data and models as well as optimization-specific parameters and settings related to constraint-handling and raster map-analysis are managed in one single control file called "config.ini".

## __Running CoMOLA__

Call <pre>
"python __init__.py" 

from the console within your CoMOLA folder. You can limit the maximum number of threads to be used for parallel computation by adding "-t x" to the command, where x is the maximum number of threads, e.g. "python __init__.py -t 2" to use two cores in maximum.