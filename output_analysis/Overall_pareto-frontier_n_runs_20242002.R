##########################################################################################################
#
#                                  Overall Pareto-Frontier Analysis
#
# Description: Identification of the overall Pareto-frontier for multiple optimization runs with CoMOLA and
#              identification of single maxima and compromise solution, plot of results. The script handles
#              m runs and n objectives.
#
#
# Authors: Andrea Kaim & Victor Steffens
# Date: 20-02-2024
#
##########################################################################################################

# Install required packages
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("viridis")
# install.packages("emoa")
# install.packages("purrr")
# install.packages("here")
# install.packages("raster")
# install.packages("Rcpp")
# install.packages("plotly")
# install.packages("reshape")
# install.packages("GGally")
# load packages

here::i_am('output_analysis/Overall_pareto-frontier_n_runs_20242002.R')

library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(emoa)
library(purrr)
library(here)
library(raster)
library(plotly)
library(reshape)
library(GGally)

#####################################
#         Set prerequisites         #
#####################################


#The order of modelnames corresponds to the order in the config.ini-file from the optimisation-process
modelnames <- c("HabStruct", "SYM", "WYLD", "SAR")
pal <- c('#FFF9C4', '#FFF176', '#FBC02D','goldenrod', 'darkgoldenrod', '#8BC34A', 'darkgreen', 'red' )
rasternames <- c('Cropland 1','Cropland 2','Cropland 3','Cropland 4','Cropland 5','Pasture','Forest','Urban')
location <- here('output') #paste the folder, where you put the results of the optimization process
#location <- "C:/Users/kaim/Nextcloud/Cloud/CoMOLA_Blockseminar/WS2324/Output_virtual_study/output" # use this if here does not work for you, paste directory of your output folder

#old function for extracting fitness values into new dataframe
#function for extracting fitness-scores for 4 objectives in given data frame, can be adjusted to 3 objectives respectively
#extract_fitness <- function(df){
#  df <- cbind.data.frame(as.numeric(gsub('\\[', '', df[,ncol(df)-3])), #for 3 objectives, change to '-2'
#                         as.numeric(df[,ncol(df)-2]), #for 3 objectives, delete line
#                         as.numeric(df[,ncol(df)-1]),
#                         as.numeric(gsub('\\]', '', df[,ncol(df)])))
#}

#####################################
#         Set functions             #
#####################################

#new function to extract fitness values into new dataframe -> no need to adjust number of variables anymore
extract_fitness <- function(bs){
  df <- data.frame()
  for (x in 1:nrow(bs)){
    df <- rbind(df, as.numeric(unlist(strsplit(as.vector(regmatches(toString(bs[x,]), gregexpr("(?<=\\[).*?(?=\\])", toString(bs[x,]), perl=T))[[1]][2]),','))))
  }
  bs <- df
}
#function for parsing a vector, consisting of 2 numbers delimited by a point, to the corresponding mapnames
vectomap <- function(v,f,l){
  runinddf <- sapply(strsplit(v,'[.]'),'[')
  runind <- sapply(runinddf[1,], FUN=function(x){tail(strsplit(f[as.numeric(x)],'/')[[1]],n=1)})
  runind <- as.data.frame(sapply(strsplit(runind,'_'),'[',1:2))
  runind <- unname(apply(as.data.frame(runind[c(1:nrow(runind)),]), MARGIN=2, paste, collapse = "_"))
  mapnames <- paste(runind,'_best_ascii_map',runinddf[2,],'.asc',sep='')
  mapnames <- paste(l, mapnames, sep='/')
  return(mapnames)
}
#function for retrieving mean values of a dataframe
which.mean <- function(df){
  compromise <- rep(c(0), times = nrow(df))
  for (x in c(1:ncol(df))){
    compromise <- compromise + abs(df[,x]-mean(df[,x]))
  }
  return(which.min(compromise))
}

#####################################
#       Read best solutions         #
#####################################

# This code creates data frames with fitness values and their corresponding solution IDs for each run and stores them in a list

# Parse all best_solution.csv-files from the output-folder into a list called files
filenames <- list.files(location, pattern = 'best_solutions.csv')

#check if only the best solutions you want to analyze are stored inside filenames, if not, then subset filenames: filenames <- filenames[c(x,x,x)]
print(filenames)
filenames <- paste(location,filenames, sep ='/')
files <- lapply(filenames, read.csv, h=F, skip=1, as.is=T)

#process list of files, the extract_fitness-function could take a while
files <- lapply(files, extract_fitness)
files <- lapply(files, setNames, modelnames)

#apply ID numbers
files <- lapply(1:length(files), \(x) transform(files[[x]], ID = paste0(x, ".", 1:nrow(files[[x]]))))


#####################################
#      Overall Pareto-Frontier      #
#####################################
#If several runs are selected, this part of the script determines the combined paretofront of all runs

# Create list with Pareto-optimal solutions of all optimization runs
S <- subset(do.call(rbind, files), select = -ID)

# transpose S
S_trans <- t(S)

# Identify non-dominated set of solutions (multiplication by -1 because nondominated_points uses minimization as default but here, we maximize)
# Check documentation of emoa package for details
P <- nondominated_points(as.matrix(S_trans * -1)) * -1
P <- as.data.frame(t(P))


#####################################
#             Get IDs               #
#####################################


# Identify IDs from all runs
joinlist <- append(list(P), files, 1)
overall_best_solutions <- reduce(joinlist, left_join, by = modelnames)

#prepare colnames and assign them to overall_best_solutions
runnumber <- c(1:length(filenames))
IDnames <- paste(rep("ID"),runnumber, sep='_')
colnames(overall_best_solutions) <- c(modelnames,IDnames)


#####################################
#         Analyze solutions         #
#####################################

##Initial Individual
#if Values for the initial individual is known, you can paste these inside the current-values-variable
#current_values <- c(1,2,3,4)

#####################################
## Single maxima
#Store the single maxima in a dataframe
Solutions <- sapply(overall_best_solutions[, c(1:length(modelnames))], FUN=function(x){which(x == max(x) )})
Solutions <- lapply(Solutions, FUN = function(x){overall_best_solutions[x,]})
Solutions <- append(Solutions, list(overall_best_solutions[which.mean(overall_best_solutions[,c(1:length(modelnames))]),]))

#Solutions <- append(list(current_values)

names(Solutions) <- c(paste('Max_', modelnames, sep=''),'Compromise'#,'Current')
)

#Extract and store ID's in list and prepare for plotting
IDs <- lapply(Solutions, FUN = function(x){x[, (length(modelnames)+1):ncol(x)]})
IDs <- lapply(IDs, FUN=function(x){x[!is.na(x)]})
IDs <- lapply(IDs, FUN=unique)
names(IDs) <- c(paste('Max_', modelnames, sep=''),'Compromise')
mapnames <- lapply(IDs, FUN=function(x){vectomap(x,filenames,location)})

maps <- lapply(mapnames, FUN=function(x){lapply(x, raster)})

#####################################
#             Plot                  #
#####################################

sol <- overall_best_solutions[,(1:length(modelnames))]

# You can change the axes on which objectives should be displayed - 
# don't forget to adapt axis labels, geom_point-data and legend-titles respectively

# Plot with x-axis HabStruct, y-axis SAR, dots WYLD, color SYM
ggplot(sol, aes(x=HabStruct, y=SAR)) +
  geom_point(aes(size = WYLD, fill=SYM), shape=21, col="black")+
  
  #highlight the maximised solutions
  geom_point(data = Solutions[['Max_HabStruct']], aes(x=HabStruct, y = SAR, shape = 'Max_HabStruct'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Max_SAR']], aes(x=HabStruct, y = SAR, shape = 'Max_SAR'), fill='red', col='black', size = 5)+
  geom_point(data = Solutions[['Max_SYM']], aes(x=HabStruct, y = SAR, shape = 'Max_SYM'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Max_WYLD']], aes(x=HabStruct, y = SAR, shape = 'Max_WYLD'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Compromise']], aes(x=HabStruct, y = SAR, shape = 'Compromise'), fill='red',col='black', size = 5)+
  #geom_point(data = Solutions[['Current']], aes(x=HabStruct, y = SAR, shape = 'Current'), fill='white', col='black', size = 5)+
  
  xlab("Habitat heterogeneity")+
  ylab("Species richness")+
  
  scale_fill_gradientn(colours = viridis(100))+
  labs(size = "Water Yield", fill = "Agricultural Yield")+
  
  scale_shape_manual(name = "Shapes",
                     values = c("Compromise" = 24,"Current" = 4, "Max_HabStruct" = 25, "Max_SAR" = 21, "Max_SYM" = 22, "Max_WYLD" = 23),
                     labels = c("Compromise", "Max_HabStruct", "Max_SAR", "Max_SYM", "Max_WYLD"))+
  
  guides(shape = guide_legend(title = "Shapes", order = 1),
         size = guide_legend(title = "Water Yield", order = 2),
         color = guide_legend(title = "Agricultural Yield", order =3))
#####################################
# Plot with x-axis HabStruct, y-axis WYLD, dots SAR, color SYM
ggplot(sol, aes(x=HabStruct, y=WYLD)) +
  geom_point(aes(size = SAR, fill=SYM), shape=21, col="black")+
  
  #highlight the maximised solutions
  geom_point(data = Solutions[['Max_HabStruct']], aes(x=HabStruct, y = WYLD, shape = 'Max_HabStruct'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Max_SAR']], aes(x=HabStruct, y = WYLD, shape = 'Max_SAR'), fill='red', col='black', size = 5)+
  geom_point(data = Solutions[['Max_SYM']], aes(x=HabStruct, y = WYLD, shape = 'Max_SYM'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Max_WYLD']], aes(x=HabStruct, y = WYLD, shape = 'Max_WYLD'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Compromise']], aes(x=HabStruct, y = WYLD, shape = 'Compromise'), fill='red',col='black', size = 5)+
  #geom_point(data = Solutions[['Current']], aes(x=HabStruct, y = SAR, shape = 'Current'), fill='white', col='black', size = 5)+
  
  xlab("Habitat heterogeneity")+
  ylab("Water Yield")+
  
  scale_fill_gradientn(colours = viridis(100))+
  labs(size = "Species Richness", fill = "Agricultural Yield")+
  
  scale_shape_manual(name = "Shapes",
                     values = c("Compromise" = 24,"Current" = 4, "Max_HabStruct" = 25, "Max_SAR" = 21, "Max_SYM" = 22, "Max_WYLD" = 23),
                     labels = c("Compromise", "Max_HabStruct", "Max_SAR", "Max_SYM", "Max_WYLD"))+
  
  guides(shape = guide_legend(title = "Shapes", order = 1),
         size = guide_legend(title = "Species Richness", order = 2),
         color = guide_legend(title = "Agricultural Yield", order = 3)
         )
#####################################
# Plot with x-axis HabStruct, y-axis SYM, dots SAR, color SYM
ggplot(sol, aes(x=HabStruct, y=SYM)) +
  geom_point(aes(size = SAR, fill=WYLD), shape=21, col="black")+
  
  #highlight the maximised solutions
  geom_point(data = Solutions[['Max_HabStruct']], aes(x=HabStruct, y = SYM, shape = 'Max_HabStruct'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Max_SAR']], aes(x=HabStruct, y = SYM, shape = 'Max_SAR'), fill='red', col='black', size = 5)+
  geom_point(data = Solutions[['Max_SYM']], aes(x=HabStruct, y = SYM, shape = 'Max_SYM'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Max_WYLD']], aes(x=HabStruct, y = SYM, shape = 'Max_WYLD'), fill='red',col='black', size = 5)+
  geom_point(data = Solutions[['Compromise']], aes(x=HabStruct, y = SYM, shape = 'Compromise'), fill='red',col='black', size = 5)+
  #geom_point(data = Solutions[['Current']], aes(x=HabStruct, y = SAR, shape = 'Current'), fill='white', col='black', size = 5)+
  
  xlab("Habitat heterogeneity")+
  ylab("Agricultural Yield")+
  
  scale_fill_gradientn(colours = viridis(100))+
  labs(size = "Species Richness", fill = "Water Yield")+
  
  scale_shape_manual(name = "Shapes",
                     values = c("Compromise" = 24,"Current" = 4, "Max_HabStruct" = 25, "Max_SAR" = 21, "Max_SYM" = 22, "Max_WYLD" = 23),
                     labels = c("Compromise", "Max_HabStruct", "Max_SAR", "Max_SYM", "Max_WYLD"))+
  
  guides(shape = guide_legend(title = "Shapes", order = 1),
         size = guide_legend(title = "Species Richness", order = 2),
         color = guide_legend(title = "Water Yield", order = 3))
#####################################
# Plot with x-axis SYM, y-axis WYLD, dots SAR, color HabStruct
ggplot(sol, aes(x=SYM, y=WYLD)) +
  geom_point(aes(size = SAR, fill=HabStruct), shape=21, col="black")+
  
  #highlight the maximised solutions
  geom_point(data = Solutions[['Max_HabStruct']], aes(x=SYM, y = WYLD, shape = 'Max_HabStruct'), fill='red', col='black', size = 5)+
  geom_point(data = Solutions[['Max_SAR']], aes(x=SYM, y = WYLD, shape = 'Max_SAR'), fill='red', col='black', size = 5)+
  geom_point(data = Solutions[['Max_SYM']], aes(x=SYM, y = WYLD, shape = 'Max_SYM'), fill='red', col='black', size = 5)+
  geom_point(data = Solutions[['Max_WYLD']], aes(x=SYM, y = WYLD, shape = 'Max_WYLD'), fill='red', col='black', size = 5)+
  geom_point(data = Solutions[['Compromise']], aes(x=SYM, y = WYLD, shape = 'Compromise'), fill='red', col='black', size = 5)+
  #geom_point(data = Solutions[['Current']], aes(x=HabStruct, y = SAR, shape = 'Current'), fill='white', col='black', size = 5)+
  
  xlab("Agricultural Yield")+
  ylab("Water Yield") +
  
  scale_fill_gradientn(colours = viridis(100))+
  labs(size = "Species richness", fill = "Habitat heterogeneity")+
  
  scale_shape_manual(name = "Shapes",
                     values = c("Compromise" = 24,"Current" = 4, "Max_HabStruct" = 25, "Max_SAR" = 21, "Max_SYM" = 22, "Max_WYLD" = 23),
                     labels = c("Compromise", "Max_HabStruct", "Max_SAR", "Max_SYM", "Max_WYLD"))+
  
  guides(shape = guide_legend(title = "Shapes", order = 1),
         size = guide_legend(title = "Species Richness", order = 2),
         color = guide_legend(title = "Habitat Heterogeneity", order =3))
#####################################

#do a Pairwise plot
ggpairs(sol)

#plot every map in maps-list
rasterplots <- lapply(seq(maps), function(i){
  lapply(seq(maps[[i]]), function(j){
    plot(maps[[i]][[j]],
         col=pal,
         main=paste(names(maps)[i], IDs[[i]][j]),
         legend=F,
         axes=F,
         box=F)
    legend(x='topright',
           inset = c(-0.5,0),
           fill = pal,
           legend = rasternames,
           border =F,
           bty ='n',
           xpd = T)
    })  
  })

# Parallel coordinates plot
# save as web page

parcord <- overall_best_solutions %>%
  plot_ly(type = 'parcoords',
          #line = list(color = 'blue'),
          labelfont=list(size=20), # font size axes
          #width = 1500, height = 800,
          dimensions = list(
            list(range = c(min(overall_best_solutions$HabStruct),max(overall_best_solutions$HabStruct)),
                 label = 'Habitat heterogeneity', values = ~HabStruct),
            list(range = c(min(overall_best_solutions$SAR),max(overall_best_solutions$SAR)),
                 label = 'Forest Species Richness', values = ~SAR),
            list(range = c(min(overall_best_solutions$SYM),max(overall_best_solutions$SYM)),
                 label = 'Agricultural yield', values = ~SYM),
            list(range = c(min(overall_best_solutions$WYLD),max(overall_best_solutions$WYLD)),
                 label = 'Water yield', values = ~WYLD)
          )
  )
parcord