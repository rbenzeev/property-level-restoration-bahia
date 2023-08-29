# Load libraries
library(tidyverse)
library(raster)
library(geobr)
library(sf)
library(exactextractr)

## Read in data
forest_cover_1985 <- raster("./Data/Forest_cover/Collection_5/mapbiomas-brazil-collection-50-bahia-1985.tif") 
forest_cover_2012 <- raster("./Data/Forest_cover/Collection_5/mapbiomas-brazil-collection-50-bahia-2012.tif") 
zipfiles <- list.files(path = "./Data/Municipalities/Shapefiles", pattern = "*.zip", full.names=TRUE, recursive = TRUE)

# merge rasters to obtain study area 1985
r1 <- raster("./Data/Forest_cover/Collection_5/Restoration_1985/mapbiomas-brazil-collection-50-bahia-1985_2019-0000000000-0000000000.tif") 
r2 <- raster("./Data/Forest_cover/Collection_5/Restoration_1985/mapbiomas-brazil-collection-50-bahia-1985_2019-0000000000-0000031744.tif") # transition represents starting date until 2018
r3 <- raster("./Data/Forest_cover/Collection_5/Restoration_1985/mapbiomas-brazil-collection-50-bahia-1985_2019-0000031744-0000000000.tif") # transition represents starting date until 2018
r4 <- raster("./Data/Forest_cover/Collection_5/Restoration_1985/mapbiomas-brazil-collection-50-bahia-1985_2019-0000031744-0000031744.tif") # transition represents starting date until 2018
x <- list(r1, r2, r3, r4)
restoration_1985 <- do.call(merge, x)
plot(restoration_2012)

# merge rasters to obtain study area 2012
r5 <- raster("./Data/Forest_cover/Collection_5/Restoration_2012/mapbiomas-brazil-collection-50-bahia-2012_2019-0000000000-0000000000.tif") 
r6 <- raster("./Data/Forest_cover/Collection_5/Restoration_2012/mapbiomas-brazil-collection-50-bahia-2012_2019-0000000000-0000031744.tif") # transition represents starting date until 2018
r7 <- raster("./Data/Forest_cover/Collection_5/Restoration_2012/mapbiomas-brazil-collection-50-bahia-2012_2019-0000031744-0000000000.tif") # transition represents starting date until 2018
r8 <- raster("./Data/Forest_cover/Collection_5/Restoration_2012/mapbiomas-brazil-collection-50-bahia-2012_2019-0000031744-0000031744.tif") # transition represents starting date until 2018
y <- list(r5, r6, r7, r8)
restoration_2012 <- do.call(merge, y)
plot(restoration_2012)

rm(r1, r2, r3, r4, r5, r6, r7, r8, x, y)

### STUDY REGION 
# prep study region plot and list
meso <- read_meso_region(code_meso= 2907, year=2019) %>% 
  st_transform(crs(forest_cover_1985))
muni <- read_municipality(code_muni= 29, year=2019) %>% 
  st_transform(crs(forest_cover_1985))
muni_south_bahia <- st_intersection(muni, meso)

### PROPERTIES  
# create a scrap directory to read in zipfiles 
scrap_dir <- "./Data/Scrap"
dir.create(scrap_dir, showWarnings = F)

# create an empty list for results
result_list <- list()
# for each zipfile: unzip twice to the scrap directory, compile results in a list, and iteratively delete each set on the computer
for (i in 1:length(zipfiles)) {
  unzip(zipfiles[i], 
        exdir = scrap_dir, 
        files="RESERVA_LEGAL.zip")
  unzip("./Data/Scrap/RESERVA_LEGAL.zip", exdir = scrap_dir) 
  result_list[[i]] <- st_read(dsn = scrap_dir)
  unlink(list.files(scrap_dir,full.names = T))
} 

# bind all municipalities from the list 
result <- do.call("rbind", result_list) %>% 
  st_transform(crs(forest_cover_1985))

# create a df of all intersecting properties with the meso region 
region_intersects <- st_intersects(result, muni_south_bahia, sparse=F) %>%
  rowSums() # used to be meso

# join dfs
# select the geometry column
# then rbind into one df
properties <- result %>% 
  filter(region_intersects > 0)  %>% 
  dplyr::select(IDF, NUM_AREA, NOM_TEMA, geometry) %>% 
  rename(Name = NOM_TEMA, Size = NUM_AREA) %>% 
  mutate(Type = "Properties") %>% 
  mutate(ID = 1:nrow(.)) %>%
  dplyr::select(ID, IDF, Type, Name, Size, geometry)

### FOREST COVER
# crop rasters to study area 
forest_cover_1985 <- crop(forest_cover_1985, muni_south_bahia) 
forest_cover_2012 <- crop(forest_cover_2012, muni_south_bahia) 
restoration_1985 <- crop(restoration_1985, muni_south_bahia)
restoration_2012 <- crop(restoration_2012, muni_south_bahia)


### ANALYSIS DATAFRAME 
# forest cover reclassification 1985
# create matrix that only keeps forest class (3) as a 1
# reclassify forest data 
# 3 classes: restoration, deforestation, no change
# forest cover: forest class of 3 == restoration 
# forest transition: anything ending in 03 == restoration 
# deforestation transition: anything starting in 03 and ending in something different 
# full list of deforestation codes: 313, 315, 321, 324, 325, 332, 333 

# Reclassify the six rasters
# legend codes: https://mapbiomas-br-site.s3.amazonaws.com/_EN__C%C3%B3digos_da_legenda_Cole%C3%A7%C3%A3o_5__1_.pdf
unique(forest_cover_1985) #  [1]  0  3  4  5  9 12 13 15 21 23 24 25 29 30 31 32 33
unique(restoration_1985) # [1]    0    3   13   15   21   23   24   25   32   33  300  303  304  305  309  313
# [17]  315  320  321  323  324  325  329  330  331  332  333  336  341  403  404  409
# [33]  412  413  415  421  424  425  429  430  433  441  500  503  505  509  513  521
# [49]  523  524  525  531  532  533  903  905  909  913  915  920  921  923  924  936
# [65]  941 1203 1204 1212 1215 1221 1229 1233 1300 1303 1304 1305 1309 1313 1315 1320
# [81] 1321 1323 1324 1325 1329 1331 1332 1333 1336 1341 1500 1503 1504 1505 1509 1512
# [97] 1513 1515 1520 1521 1523 1524 1525 1529 1530 1531 1532 1533 1536 1541 2100 2103
# [113] 2104 2105 2109 2112 2113 2115 2120 2121 2123 2124 2125 2129 2130 2131 2132 2133
# [129] 2136 2141 2300 2303 2305 2309 2313 2315 2321 2323 2324 2325 2333 2423 2424 2433
# [145] 2500 2503 2504 2505 2509 2513 2515 2520 2521 2523 2524 2525 2529 2531 2532 2533
# [161] 2536 2541 2903 2904 2909 2913 2915 2921 2925 2929 2936 3025 3030 3100 3115 3123
# [177] 3125 3131 3132 3133 3200 3203 3205 3213 3215 3221 3223 3224 3225 3232 3233 3300
# [193] 3303 3304 3305 3309 3312 3313 3315 3320 3321 3323 3324 3325 3329 3331 3332 3333
# [209] 3336

forest_cover_1985[forest_cover_1985 == 3] <- 1
forest_cover_1985[forest_cover_1985 != 1] <- 0
forest_cover_2012[forest_cover_2012 == 3] <- 1
forest_cover_2012[forest_cover_2012 != 1] <- 0

# Numbers we want: 15, 19, 39, 20, 41, 36, 21, 9
# Numbers that are not in the dataset from above: 19, 39, 20, 41, 36
# This leaves us with 15, 21, and 9 
restoration_1985[restoration_1985 == 1503 |  # pasture 
                   restoration_1985 == 2103 | # mosaic of agriculture and pasture
                   restoration_1985 == 0903 ] <- 1 # forest plantation
restoration_1985[restoration_1985 != 1] <- 0 

restoration_2012[restoration_2012 == 1503 |  # pasture 
                   restoration_2012 == 2103 |  # mosaic of agriculture and pasture
                   restoration_2012 == 0903 ] <- 1 # forest plantation
restoration_2012[restoration_2012 != 1] <- 0 

### EXACTEXTRACTR 
properties$Forest_cover_1985 <- exact_extract(forest_cover_1985, properties, 'sum')
properties$Forest_cover_2012 <- exact_extract(forest_cover_2012, properties, 'sum')
properties$Percent_forest_1985 <- exact_extract(forest_cover_1985, properties, 'mean')
properties$Percent_forest_2012 <- exact_extract(forest_cover_2012, properties, 'mean')

properties$Restoration_1985 <- exact_extract(restoration_1985, properties, 'sum')
properties$Restoration_2012 <- exact_extract(restoration_2012, properties, 'sum')
properties$Percent_restoration_1985 <- exact_extract(restoration_1985, properties, 'mean')
properties$Percent_restoration_2012 <- exact_extract(restoration_2012, properties, 'mean')

# rename property size column
# create Size_original which represents the automatic sizes calculated from SICAR
# create Size which is st_area sizes 
# conversions: 1.75 pixels * 900 m^2/pixel 
# convert forest variables to m^2
# multiply all proportions by 100
# convert percent variables to %s
forest_dat_original <- properties %>% 
  dplyr::rename(Size_original = Size) %>% 
  mutate(Size = st_area(geometry)) %>% 
  mutate(Size = as.numeric(Size)) %>% 
  mutate(Size = Size * 0.0001) %>% 
  mutate(Forest_cover_1985 = Forest_cover_1985 * 900) %>% 
  mutate(Forest_cover_2012 = Forest_cover_2012 * 900) %>% 
  mutate(Restoration_1985 = Restoration_1985 * 900) %>% 
  mutate(Restoration_2012 = Restoration_2012 * 900) %>% 
  mutate(Percent_forest_1985 = Percent_forest_1985 * 100) %>% 
  mutate(Percent_forest_2012 = Percent_forest_2012 * 100) %>% 
  mutate(Percent_restoration_1985 = Percent_restoration_1985 * 100) %>% 
  mutate(Percent_restoration_2012 = Percent_restoration_2012 * 100) 

rm(muni, meso, result, result_list, restoration_1985, restoration_2012, forest_cover_1985, forest_cover_2012)

# add controls 

## PRECIP
precip <- getData('worldclim', var='prec', res=2.5)
precip_crop <- crop(precip, muni_south_bahia) 
precip_mean <- calc(precip_crop, mean, na.rm=T)
# average of the average precip per property
forest_dat_original$Avg_precip <- exact_extract(precip_mean, properties, 'mean')

# create centroid df of points data
forest_dat <- forest_dat_original %>% 
  st_centroid() %>% 
  st_intersection(muni_south_bahia)

### ROADS
# distance from property centroid
# units meters
Roads <- st_read(dsn = "./Data/Controls/Roads/GRIP4_Region2_vector_shp/GRIP4_region2.shp") %>% 
  st_intersection(muni_south_bahia) %>% 
  st_transform(crs(forest_cover_1985)) 

forest_dat <- forest_dat %>% 
  mutate(Distance_roads = st_distance(forest_dat, Roads, by_element = T)) 

# create df without restoration=0 values
forest_points <- forest_dat %>% 
  filter(Restoration_1985 !=0) 

# add distance to roads variable
forest_points <- forest_points %>% 
  mutate(Distance_roads = as.numeric(Distance_roads))

rm(forest_dat_original, precip, precip_crop, precip_mean, properties, Roads)

# sub-region analysis
muni_cacao <- muni_south_bahia %>% 
  dplyr::select(code_muni, name_muni, geom) %>% 
  filter(name_muni != "Alcobaça", name_muni != "Caravelas", name_muni != "Eunápolis",
         name_muni != "Guaratinga", name_muni != "Ibirapuã", name_muni != "Itabela",
         name_muni != "Itagimirim", name_muni != "Itamaraju", name_muni != "Itanhém", 
         name_muni != "Jucuruçu", name_muni != "Lajedão", name_muni != "Medeiros Neto",
         name_muni != "Mucuri", name_muni != "Nova Viçosa", name_muni != "Porto Seguro",
         name_muni != "Prado", name_muni != "Santa Cruz Cabrália", name_muni != "Teixeira De Freitas",
         name_muni != "Vereda")

forest_points_subset <- forest_points %>% 
  st_intersection(muni_cacao)

# create df where there are no restoration=0 values
forest_dat_subset <- forest_dat_original %>% 
  filter(Restoration_1985 !=0) 

# as.numeric
x <- t(as.data.frame(lapply(forest_dat_subset$geom, as.numeric)))
forest_dat_subset <- forest_dat_subset %>% 
  mutate(x1 = x[,1], x2 = x[,2])
# na.omit
forest_dat_subset <- distinct(forest_dat_subset, x1, x2, .keep_all=TRUE) %>% na.omit()

