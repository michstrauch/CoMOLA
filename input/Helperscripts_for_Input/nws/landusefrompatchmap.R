# 
# This Code generates a landusemap, patchmap and Soil-fertility-map based on
# 1.patchmap 
# 2.soilfertility-map 
# 3.landuse_grid
# given by ArcGis and makes them suitable for the optimasationprocess in CoMOLA

here::i_am('landusefrompatchmap.R')
library(here)
library(raster)


safelocation <- '/Users/victorsteffens/Documents/CoMOLA/CoMOLA_repo/COMOla/input'
location <- here()
landuse <- raster(paste(location,'nws_patch_id.asc', sep='/'),NAflag = -2, datatype = 'INT2U')
landgrid <- read.table(paste(location, 'nws_grid_landuse.csv', sep='/'), header = T, sep =';', dec=',')
soilfert <- raster(paste(location,'soil_fertility.asc', sep='/'),NAflag = -2, datatype = 'INT2U')

# Find values for landuse by using landgrid
landuseindices <- match(values(landuse), landgrid$grid_id)

# Replace patchmapvalues with corresponding values from landgrid and nullify all 8-landuseclasses. Then write new landusemap.
# Addiditionally adapt soilfertility-map
values(landuse) <- landgrid$lu_code[landuseindices]

values(landuse)[is.na(values(landuse))] <- -2

nullsoilindices <- which(values(landuse) %in% -2)
staticindices <- which(values(landuse) %in% 8)

values(soilfert)[nullsoilindices] <- -2

writeRaster(landuse, paste(location,'landusefrompatch.asc', sep='/'), overwrite = T, NAflag = -2, datatype = 'INT')
writeRaster(soilfert, paste(location,'soilfertilityfrompatch.asc', sep='/'), overwrite = T)

#Load patchmap again to create patchmap corresponding to landusemap:
patchmap <- raster(paste(location,'nws_patch_id.asc', sep='/'),NAflag = -2, datatype = 'INT2U')
#1)Nullify -2-Values and 8-Values
#2)Build consecutive order for numbering of patches
modlandgrid <- landgrid[!landgrid$lu_code=='-2',]
modlandgrid <- modlandgrid[!modlandgrid$lu_code == '8',]
modlandgrid <- modlandgrid[,-1]
rownames.modlandgrid <- NULL
#modlandgrid$FID <- c(1:nrow(modlandgrid))
orderpatchindices <- match(values(patchmap), modlandgrid$grid_id)

values(patchmap) <- orderpatchindices
values(patchmap)[staticindices] <- 0

values(patchmap)[is.na(values(patchmap))] <- 0

writeRaster(patchmap, paste(location,'newpatchfromoldpatch.asc', sep='/'), overwrite = T, NAflag = -2, datatype = 'INT')

