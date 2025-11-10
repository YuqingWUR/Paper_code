
library(data.table)
library(terra)
library(SoilR)
library(sf)
library(parallel)
library(doParallel)
library(foreach)


# load required packages
pacman::p_load(terra, SoilR, parallel, doParallel, sf, data.table)


## defining functions  ***************************************
# function1: create a general function to read data
load_raster_data <- function(directory, pattern){
  files <- list.files(path = directory, pattern = pattern, full.names = TRUE)
  stack <- rast(files)
  return(stack)
}
fw1func <- function(P, E, S.Thick = 30, pClay = clay, pE = 1, bare) {
  M = P - E * pE
  Acc.TSMD = NULL
  for (i in 2:length(M)) {
    B = ifelse(bare[i] == FALSE, 1, 1.8)
    Max.TSMD = -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (S.Thick/23) * (1/B)
    Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
    if (Acc.TSMD[i - 1] + M[i] < 0) {
      Acc.TSMD[i] = Acc.TSMD[i - 1] + M[i]
    } else {
      Acc.TSMD[i] = 0
    }
    if (Acc.TSMD[i] <= Max.TSMD) {
      Acc.TSMD[i] = Max.TSMD
    }
  }
  b <- ifelse(Acc.TSMD > 0.444 * Max.TSMD,
              1,
              (0.2 + 0.8 * ((Max.TSMD - Acc.TSMD)/(Max.TSMD - 0.444 * Max.TSMD))))
  
  b <- terra::clamp(b, lower=0.2)  
  return(data.frame(b))
}

# function3: Vegetation Cover effects: 1 indicate bare soil, 2 indicate vegetated soil. 
RMF_SC <- function(SC){
  
  fC <- ifelse(SC < 0.2, 1, 0.6)
  
  return(fC)
}

yuqing_path <- 'W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version1_20250418/03_RData_forward'

## read data for RothC model 
Cin <- load(file.path(yuqing_path, 'region3/region3_crop_extract.RData'))
Cin_dfr1 <- crop_dfr3[1:1000000, ]
rm(crop_dfr3)

manure <- load("W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/DATA/03_Carbon_input/Manure_data_new_20251002/Manure_RData/manure_stack3.RData")
manure_dfr1 <- df_manure[1:1000000, ]
rm(df_manure)

ndvi <- load(file.path(yuqing_path, 'region3/region3_ndvi_extract.RData'))
ndvi_dfr1v <- ndvi_dfr3[1:1000000, ]
rm(ndvi_dfr3)

precip <- load(file.path(yuqing_path, 'region3/region3_precip_extract.RData'))
precip_dfr1v <- precip_dfr3[1:1000000, ]
rm(precip_dfr3)

temp <- load(file.path(yuqing_path, 'region3/region3_temp_extract.RData'))
temp_dfr1v <- temp_dfr3[1:1000000, ]
rm(temp_dfr3)

pet <- load(file.path(yuqing_path, 'region3/region3_pet_extract.RData'))
pet_dfr1v <- pet_dfr3[1:1000000, ]
rm(pet_dfr3)

clay <- load(file.path(yuqing_path, 'region3/region3_clay_extract.RData'))
clay_dfr1v <- clay_dfr3[1:1000000, ]
rm(clay_dfr3)

DR <- load(file.path(yuqing_path, 'region3/region3_DR_extract.RData'))
DR_dfr1v <- DR_dfr3[1:1000000, ]
rm( DR_dfr3)
DR_dfr1v <- DR_dfr1v/100

mr <- load("W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/DATA/mr_RData/mr3.RData")
names(df_mr) <- "mr"
df_mr1v <- df_mr[1:1000000, ]
rm(df_mr)

DR_dfr1v[is.na(DR_dfr1v)] <- 1.2995
Cin_dfr1[is.na(Cin_dfr1)] <- 0
df_mr1v[is.na(df_mr1v)] <- 1.08
manure_dfr1[is.na(manure_dfr1)] <- 0

spin_up1 <- st_read('W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version2_20251010/spinup/spinup_sf/spinup_region3.shp')
head(spin_up1)
spin_up11 <- spin_up1[1:1000000,2:6]
coords <- st_coordinates(spin_up11)
rm(spin_up1)
names(ndvi_dfr1v)[114] <- "ndvi_1995_05" 
names(ndvi_dfr1v)[429] <- "ndvi_2021_08" 

data <- cbind(clay_dfr1v, DR_dfr1v, temp_dfr1v, pet_dfr1v,precip_dfr1v,Cin_dfr1,manure_dfr1,ndvi_dfr1v,df_mr1v, spin_up11, coords)

data <- st_drop_geometry(data) 
data1 <- data[,-c(2718:2724)]

dim(data1)
data_noge <- data1
data_noge1 <- data_noge

# Initialize an empty data.table for results
results_all_GH <- data.table()

# Initialize an empty data.table for results
results_all_GH <- data.table()

# Loop over each point
for (i in 1:nrow(data_noge1)) {
  print(paste0('processing point:', i))
  
  C0 <- as.numeric(data_noge1[i, grep("DPM|RPM|BIO|HUM|IOM", names(data_noge1))])
  clay <- as.numeric(data_noge1[i, grep("clay", names(data_noge1))])
  Cin <- as.numeric(data_noge1[i, grep("Cin", names(data_noge1))])
  manure <- as.numeric(data_noge1[i, grep("manure", names(data_noge1))])
  r <- as.numeric(data_noge1[i, grep("DR", names(data_noge1))])
  mr <- as.numeric(data_noge1[i, grep("mr", names(data_noge1))])
  temp <- as.numeric(data_noge1[i, grep("temp", names(data_noge1))])
  Precip <- as.numeric(data_noge1[i, grep("precip", names(data_noge1))])
  pet <- as.numeric(data_noge1[i, grep("PET", names(data_noge1))])
  ndvi <- as.numeric(data_noge1[i, grep("ndvi", names(data_noge1))])
  
  bare <- ndvi < 0.2
  X <- as.numeric(data_noge1[i, grep("X", names(data_noge1))])
  Y <- as.numeric(data_noge1[i, grep("Y", names(data_noge1))])
  
  # Check for missing data
  if (any(is.na(c(C0, clay, Cin, manure, r, temp, Precip, pet, ndvi, X, Y)))) {
    results <- data.table(id = i, year = NA, X = X, Y = Y)
    results[, c("DPM", "RPM", "BIO", "HUM", "IOM", "total_C") := NA]
  } else {
    # RothC model calculations
    fT <- fT.RothC(as.numeric(temp))
    fW_2 <- fw1func(P = as.numeric(Precip), E = as.numeric(pet), S.Thick = 30, pClay = clay, pE = 1, bare = bare)$b 
    fC <- RMF_SC(as.numeric(ndvi))
    xi <- fT * fW_2 * fC
    timesteps <- length(Cin)
    rho <- rep(1/(r+1), each = 12)
    tau <- rep(1/(mr+1), times=timesteps) 
    # 
    # manure_vector <- c(tau, 1-tau, 0, 0.02, 0)
    
    cx <- 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
    alpha_BIO <- 0.46 / (cx + 1)
    alpha_HUM <- 0.54 / (cx + 1)
    
    ks <- c(10, 0.3, 0.66, 0.02, 0)  
    
    A <- diag(-ks)
    ai3 <- alpha_BIO * ks
    ai4 <- alpha_HUM * ks 
    A[3, ] <- A[3, ] + ai3
    A[4, ] <- A[4, ] + ai4
    A <- A / 12
    
    
    C <- matrix(0, ncol = 5, nrow = timesteps + 1)
    C[1, ] <- as.numeric(C0)
    
    for (t in 1:timesteps) {
      manure_vector <- c(0.98*(1 - tau[t]), 0.98*tau[t], 0, 0.02, 0)   # ✅ 长度固定 5
      Q <- c(1 - rho[t], rho[t], 0, 0, 0) * as.numeric(Cin[t]) + manure_vector * as.numeric(manure[t])
      C[t + 1, ] <- C[t, ] + Q + xi[t] * A %*% C[t, ]
    }
    
    year_indices <- seq(12, timesteps, by = 12)
    
    if (length(year_indices) > 0) {
      C_sub <- as.matrix(C[year_indices, ])
      results <- data.table(
        id = i,
        year = seq(1, length(year_indices)),
        X = X,
        Y = Y,
        DPM = round(C_sub[, 1], 3),
        RPM = round(C_sub[, 2], 3),
        BIO = round(C_sub[, 3], 3),
        HUM = round(C_sub[, 4], 3),
        IOM = round(C_sub[, 5], 3),
        total_C = round(rowSums(C_sub, na.rm = TRUE), 3)
      )
    } else {
      results <- data.table(id = i, year = NA, X = X, Y = Y)
      results[, c("DPM", "RPM", "BIO", "HUM", "IOM", "total_C") := NA]
    }
  }
  
  # Append results
  results_all_GH <- rbindlist(list(results_all_GH, results), fill = TRUE)
  # 
  if (i %% 10000 == 0) {
    file_path <- paste0("W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version2_20251010/forward/CSV/point1/forward_region_", i, ".csv")
    fwrite(results_all_GH, file_path)
    results_all_GH <- data.table()
    rm(results)
    gc()
  }
}


fwrite(results_all_GH, "W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version2_20251010/forward/CSV/point1/forward_region3_1final.csv")


rm(results_all_GH, data, Cin, manure, ndvi, precip, temp, pet, clay, DR, spin_up1)
gc()

print(head(results_all_GH))

