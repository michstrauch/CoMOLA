# Required packages
# install.packages("raster")
# install.packages("mco")
# install.packages("plyr")
# install.packages("ggplot2")
# install.packages("viridis")

library(mco)
library(plyr)
library(raster)
library(ggplot2)
library(viridis)
library(here)

## Define paths
# to optimization results
opt_path <- here('output')

# to some post-processing folder
post_path <- here('output_analysis')

# name of output metrics file
file_out <- "metrics.txt"

### Extract results if any results are available and feasible
setwd(opt_path)
logfile<-file(dir(pattern="_log"),"r+")
text<-readLines(logfile)
close(logfile)

if (any(grep("The optimization process needed", text) > 0)) {
  seconds <- strsplit(text[(grep("The optimization process needed", text))], "[|]") [[1]] [2]
}

if (any(grep("Best Solutions", text) > 0)) {
  bestsol <- as.matrix(text[(grep("Best Solutions", text)+2):(length(text)-1)])
  infeasible <- grep("infeasible", bestsol)
  if (length(infeasible)>0) {
    feasible <- as.matrix(c(1:length(bestsol))[-infeasible])
    bestsol <- as.matrix(bestsol[-infeasible])
  } else (feasible <- c(1:length(bestsol)))
  if (length(bestsol)>0) {
    for(k in 1:length(bestsol)){
      bestsol[k] <- strsplit(bestsol, "[|]") [[k]] [2]
    }
  }
  bestind <- bestsol
  bestfit <- bestsol
  if (length(bestsol)>0) {
    for(k in 1:length(bestsol)){
      bestind[k] <- strsplit(bestsol, "[:]") [[k]] [1]
      bestind[k] <- substr(bestind[k], 3, nchar(bestind[k])-2) 
      bestfit[k] <- strsplit(bestsol, "[:]") [[k]] [2]
      bestfit[k] <- substr(bestfit[k], 3, nchar(bestfit[k])-1) 
    }
  }
  
  if(length(bestsol)>0){
    setwd(post_path)
    write.table(file="bestind.txt",bestind, col.names=F, row.names=F, quote=F)
    write.table(file="bestfit.txt",bestfit, col.names=F, row.names=F, quote=F)
    write.table(file="seconds.txt",seconds, col.names=F, row.names=F, quote=F)
    write.table(file="feasible.txt",feasible, col.names=F, row.names=F, quote=F)
  }
}

### get maps for single maxima and compromise solution
# read in fitness values of best solutions
best.sol <- read.table("bestfit.txt", h=F, sep =",")
nobj <- length(best.sol)
best.sol$ID <- c(1:dim(best.sol)[1])

if(nobj==2){
  # extract only feasible solutions
  best.sol <- best.sol[which(best.sol$ID %in% feasible),]
  names(best.sol) <- c("obj1","obj2","ID")
  
  setwd(opt_path)
  # objective 1
  sol_max1 <- best.sol[which.max(best.sol$obj1),c(1:2)]
  ID_max1 <- best.sol[which.max(best.sol$obj1),3]
  map_max1 <- read.asc(dir(pattern=paste("map",ID_max1,".asc",sep="")), gz = FALSE)
  
  # objective 2
  sol_max2 <- best.sol[which.max(best.sol$obj2),c(1:2)]
  ID_max2 <- best.sol[which.max(best.sol$obj2),3]
  map_max2 <- read.asc(dir(pattern=paste("map",ID_max2,".asc",sep="")), gz = FALSE)
  
  # Compromise solution (sum of deviations from single mean values is minimal)
  which.mean <- function(x,y) which.min(abs(x - mean(x))+abs(y - mean(y)))
  sol_compromise <- best.sol[which.mean(best.sol$obj1,best.sol$obj2),c(1:2)]
  ID_compromise <- best.sol[which.mean(best.sol$obj1,best.sol$obj2),3]
  map_compromise <- read.asc(dir(pattern=paste("map",ID_compromise,".asc",sep="")), gz = FALSE)
  
  setwd(post_path)
  
  write.asc(map_max1, "map_max1.asc")
  write.asc(map_max2, "map_max2.asc")
  write.asc(map_compromise, "map_compromise.asc") 
}

if(nobj==3){
  # extract only feasible solutions
  best.sol <- best.sol[which(best.sol$ID %in% feasible),]
  names(best.sol) <- c("obj1","obj2","obj3","ID")
  
  ### get solutions for single maxima
  setwd(opt_path)
  # objective 1
  sol_max1 <- best.sol[which.max(best.sol$obj1),c(1:3)]
  ID_max1 <- best.sol[which.max(best.sol$obj1),4]
  map_max1 <- read.asc(dir(pattern=paste("map",ID_max1,".asc",sep="")), gz = FALSE)
  
  # objective 2
  sol_max2 <- best.sol[which.max(best.sol$obj2),c(1:3)]
  ID_max2 <- best.sol[which.max(best.sol$obj2),4]
  map_max2 <- read.asc(dir(pattern=paste("map",ID_max2,".asc",sep="")), gz = FALSE)
  
  # objective 3
  sol_max3 <- best.sol[which.max(best.sol$obj3),c(1:3)]
  ID_max3 <- best.sol[which.max(best.sol$obj3),4]
  map_max3 <- read.asc(dir(pattern=paste("map",ID_max3,".asc",sep="")), gz = FALSE)
  
  # Compromise solution (sum of deviations from single mean values is minimal)
  which.mean <- function(x,y) which.min(abs(x - mean(x))+abs(y - mean(y))+abs(z - mean(z)))
  sol_compromise <- best.sol[which.mean(best.sol$obj1,best.sol$obj2,best.sol$obj3),c(1:3)]
  ID_compromise <- best.sol[which.mean(best.sol$obj1,best.sol$obj2,best.sol$obj3),4]
  map_compromise <- read.asc(dir(pattern=paste("map",ID_compromise,".asc",sep="")), gz = FALSE)
  
  setwd(post_path)
  
  write.asc(map_max1, "map_max1.asc")
  write.asc(map_max2, "map_max2.asc")
  write.asc(map_max3, "map_max3.asc")
  write.asc(map_compromise, "map_compromise.asc") 
}

if(nobj==4){
  # extract only feasible solutions
  best.sol <- best.sol[which(best.sol$ID %in% feasible),]
  names(best.sol) <- c("obj1","obj2","obj3","obj4","ID")
  
  ### get solutions for single maxima
  setwd(opt_path)
  # objective 1
  sol_max1 <- best.sol[which.max(best.sol$obj1),c(1:4)]
  ID_max1 <- best.sol[which.max(best.sol$obj1),5]
  map_max1 <- read.asc(dir(pattern=paste("map",ID_max1,".asc",sep="")), gz = FALSE)
  
  # objective 2
  sol_max2 <- best.sol[which.max(best.sol$obj2),c(1:4)]
  ID_max2 <- best.sol[which.max(best.sol$obj2),5]
  map_max2 <- read.asc(dir(pattern=paste("map",ID_max2,".asc",sep="")), gz = FALSE)
  
  # objective 3
  sol_max3 <- best.sol[which.max(best.sol$obj3),c(1:4)]
  ID_max3 <- best.sol[which.max(best.sol$obj3),5]
  map_max3 <- read.asc(dir(pattern=paste("map",ID_max3,".asc",sep="")), gz = FALSE)
  
  # objective 4
  sol_max4 <- best.sol[which.max(best.sol$obj4),c(1:4)]
  ID_max4 <- best.sol[which.max(best.sol$obj4),5]
  map_max4 <- read.asc(dir(pattern=paste("map",ID_max4,".asc",sep="")), gz = FALSE)
  
  # Compromise solution (sum of deviations from single mean values is minimal)
  which.mean <- function(x,y) which.min(abs(x - mean(x))+abs(y - mean(y))+abs(z - mean(z))+abs(a - mean(a)))
  sol_compromise <- best.sol[which.mean(best.sol$obj1,best.sol$obj2,best.sol$obj3,best.sol$obj4),c(1:4)]
  ID_compromise <- best.sol[which.mean(best.sol$obj1,best.sol$obj2,best.sol$obj3,best.sol$obj4),5]
  map_compromise <- read.asc(dir(pattern=paste("map",ID_compromise,".asc",sep="")), gz = FALSE)
  
  setwd(post_path)
  
  write.asc(map_max1, "map_max1.asc")
  write.asc(map_max2, "map_max2.asc")
  write.asc(map_max3, "map_max3.asc")
  write.asc(map_max4, "map_max4.asc")
  write.asc(map_compromise, "map_compromise.asc") 
}
 
### Calculate performance measures
### Please note, they are only meaningful when comparing different runs

is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

calcMetrics <- function(path,file_out){
  setwd(path)
  
  # get fitness values
  sol <- read.table("bestfit.txt",sep=",",h=F,as.is=T)
  
  ## calculate performance metrics for a 2-objective optimization
  if(length(sol)==2) {
    
    ref <- matrix(data=0,nrow=1, ncol=2)
    
    metrics <- data.frame(matrix(data=NA, nrow=1, ncol=5))
    names(metrics) <- c("n","max1","max2","dHV","sec")
    
    # get number of solutions
    metrics[1,1] <- dim(sol)[1]
    
    # get single maxima
    metrics[1,2] <- round(max(sol[,1]),4)
    metrics[1,3] <- round(max(sol[,2]),4)
    
    # normalize for calculating the hypervolume
    sol.n <- sol
    sol.n[,1] <- sol[,1]/max(sol[,1])
    sol.n[,2] <- sol[,2]/max(sol[,2])
    
    # get dominated Hypervolume
    sol.val_ <- as.matrix(sol.n[,c(1:2)]*-1)
    pF <- paretoFront(sol.val_)
    pF_pos <- paretoFront(as.matrix(sol.n[,c(1:2)]))
    metrics[1,4] <- round(dominatedHypervolume(pF, ref),4)
    
    # get run time (seconds)
    sec <- read.table("seconds.txt",sep=" ",h=F,as.is=T)
    metrics[1,5] <- sec[6]
    
    write.table(metrics, file=file_out, col.names=T, row.names=F, quote=F)
  }
  
  ## calculate performance metrics for a 3-objective optimization
  if(length(sol)==3) {
    
    ref <- matrix(data=0,nrow=1, ncol=3)
    
    metrics <- data.frame(matrix(data=NA, nrow=1, ncol=6))
    names(metrics) <- c("n","max1","max2","max3","dHV","sec")
    
    # get number of solutions
    metrics[1,1] <- dim(sol)[1]
    
    # get single maxima
    metrics[1,2] <- round(max(sol[,1]),4)
    metrics[1,3] <- round(max(sol[,2]),4)
    metrics[1,4] <- round(max(sol[,3]),4)
    
    # normalize for calculating the hypervolume
    sol.n <- sol
    sol.n[,1] <- sol[,1]/max(sol[,1])
    sol.n[,2] <- sol[,2]/max(sol[,2])
    sol.n[,3] <- sol[,3]/max(sol[,3])
    
    # get dominated Hypervolume
    sol.val_ <- as.matrix(sol.n[,c(1:3)]*-1)
    pF <- paretoFront(sol.val_)
    pF_pos <- paretoFront(as.matrix(sol.n[,c(1:3)]))
    metrics[1,5] <- round(dominatedHypervolume(pF, ref),4)
    
    # get run time (seconds)
    sec <- read.table("seconds.txt",sep=" ",h=F,as.is=T)
    metrics[1,6] <- sec[6]
    
    write.table(metrics, file=file_out, col.names=T, row.names=F, quote=F)
  }
  ## calculate performance metrics for a 4-objective optimization
  if(length(sol)==4) {
    
    ref <- matrix(data=0,nrow=1, ncol=4)
    
    metrics <- data.frame(matrix(data=NA, nrow=1, ncol=7))
    names(metrics) <- c("n","max1","max2","max3","max4","dHV","sec")
    
    # get number of solutions
    metrics[1,1] <- dim(sol)[1]
    
    # get single maxima
    metrics[1,2] <- round(max(sol[,1]),4)
    metrics[1,3] <- round(max(sol[,2]),4)
    metrics[1,4] <- round(max(sol[,3]),4)
    metrics[1,5] <- round(max(sol[,4]),4)
    
    # normalize for hypervolume calculation
    sol.n <- sol
    sol.n[,1] <- sol[,1]/max(sol[,1])
    sol.n[,2] <- sol[,2]/max(sol[,2])
    sol.n[,3] <- sol[,3]/max(sol[,3])
    sol.n[,4] <- sol[,4]/max(sol[,4])
    
    # get dominated Hypervolume
    sol.val_ <- as.matrix(sol.n[,c(1:4)]*-1)
    pF <- paretoFront(sol.val_)
    pF_pos <- paretoFront(as.matrix(sol.n[,c(1:4)]))
    metrics[1,6] <- round(dominatedHypervolume(pF, ref),4)
    
    # get run time (seconds)
    sec <- read.table("seconds.txt",sep=" ",h=F,as.is=T)
    metrics[1,7] <- sec[6]
    
    write.table(metrics, file=file_out, col.names=T, row.names=F, quote=F)
  }
}

# apply function
calcMetrics(post_path, file_out)


### Plot results
# 2D version, maybe only an option if there are not too many solutions
# (if there are too many solutions you may use a filter)
setwd(post_path)
sol <- read.table("bestfit.txt",sep=",",h=F,as.is=T)

# plot for 2 objectives
if(length(sol)==2){
  names(sol) <- c("obj1","obj2")
  ggplot(sol, aes(x=obj1, y=obj2)) + 
    geom_point(aes(), shape=21, col="black") +
    xlab("Objective 1")+
    ylab("Objective 2")
}

# plot for 3 objectives
if(length(sol)==3){
  names(sol) <- c("obj1","obj2","obj3")
  ggplot(sol, aes(x=obj1, y=obj2)) + 
    geom_point(aes(fill=obj3), shape=21, col="black") +
    scale_fill_gradientn(colours=viridis(100)) +
    guides(fill = guide_legend(title="Objective 3")) +
    xlab("Objective 1")+
    ylab("Objective 2")
}

# plot for 4 objectives
if(length(sol)==4){
  names(sol) <- c("obj1","obj2","obj3","obj4")
  ggplot(sol, aes(x=obj1, y=obj2)) + 
    geom_point(aes(size=obj4, fill=obj3), shape=21, col="black") +
    scale_fill_gradientn(colours=viridis(100)) +
    guides(size = guide_legend(title="Objective 4", override.aes = list(col = "black")),
           fill = guide_legend(title="Objective 3")) +
    xlab("Objective 1")+
    ylab("Objective 2")
}

