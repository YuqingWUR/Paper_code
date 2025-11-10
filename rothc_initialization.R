
  

## load required packages
pacman::p_load(terra, sf, dplyr, deSolve, SoilR, deSolve)


#### DEFINING RATE MODIFYING FACTOR FUNCTION FOR RASTER DATA ------------------

calculate_rate_modifiers <- function(Temp, Precip, Evp, soil_thick, clay, bare, Cov, pE = 1) {
  
  # # temperature
  # RMF_TMP <- function(Temp) {
  #   RM_TMP <- numeric(length(Temp))  
  #   for (m in seq_along(Temp)) {
  #     if (Temp[m] < -5) { 
  #       RM_TMP[m] <- 0
  #     } else {
  #       RM_TMP[m] <- 47.91 / (exp(106.06 / (Temp[m] + 18.27)) + 1.0)
  #     }
  #   }
  #   return(RM_TMP)
  # }
  # fT <- fT.RothC(as.numeric(Temp)) # temperature effects
  # soil moisture
  RMF_Moist <- function(P, E, S.Thick = 30, pClay, pE = 1, bare) {
    M = P - E * pE
    Acc.TSMD = numeric(length(M))
    
    for (r in 2:length(M)) {
      B = ifelse(bare[r] == FALSE, 1, 1.8)
      Max.TSMD = -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (S.Thick / 23) * (1 / B)
      Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
      if (Acc.TSMD[r - 1] + M[r] < 0) {
        Acc.TSMD[r] = Acc.TSMD[r - 1] + M[r]
      } else {
        Acc.TSMD[r] = 0
      }
      if (Acc.TSMD[r] < Max.TSMD) {
        Acc.TSMD[r] = Max.TSMD
      }
    }
    
    b = ifelse(Acc.TSMD > 0.444 * Max.TSMD, 1, (0.2 + 0.8 * ((Max.TSMD - Acc.TSMD) / (Max.TSMD - 0.444 * Max.TSMD))))
    b <- pmax(0.2, b)  # force to set the minimum b value to 0.2
    
    return(b)
  }
  
  # soil cover
  RMF_SC <- function(PC) {
    RM_SC <- ifelse(PC < 0.2, 1, 0.6)
    return(RM_SC)
  }
  
  # calculate each rate modifying factor
  # RM_TMP <- RMF_TMP(Temp[,2])
  RM_TMP <- fT.RothC(as.numeric(Temp$temp))
  RM_Moist <- RMF_Moist(P = Precip[,2], E = Evp[,2], S.Thick = soil_thick, pClay = clay, pE = 1, bare = bare)
  RM_PC <- RMF_SC(Cov)
  
  # calculate the totoal RateM
  RateM <- RM_TMP * RM_Moist * RM_PC
  # mean_RateM <- mean(RateM)
  
  # return RateM
  return(list(RateM = RateM))
}


#####  01_load point data ------------------------------------
# load point data, it may takes 4mins
sf_point1 <- st_read('W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version1_20250418/02_data_spinup/03_point_data/points_part_2.shp')

## load covariates
yuqing <- "W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version1_20250418/02_data_spinup"
setwd(paste(yuqing))

## set working directory
precip_wd <- c(paste(yuqing, "/05_mon_avg_precip_61_86", sep=""))
temp_wd <- c(paste(yuqing, "/06_mon_avg_temp_61_86", sep=""))
pet_wd <- c(paste(yuqing, "/04_mon_avg_pet_1961_1986", sep=""))
ndvi_wd <- c(paste(yuqing, "/02_monthly_ndvi_1986", sep=""))

## (1) load temperature
setwd(temp_wd)
temp_files <- list.files(pattern = '.tif')
temp_files
temp_stack <- rast(temp_files)
temp_stack
# plot(temp_stack[[1:3]])
names(temp_stack)
names(temp_stack) <- c('temp_01', 'temp_02', 'temp_03', 'temp_04', 'temp_05', 'temp_06', 'temp_07', 'temp_08', 'temp_09', 'temp_10', 'temp_11', 'temp_12')


## (2) load precipitation
setwd(precip_wd)
precip_files <- list.files(pattern = '.tif')
precip_files
precip_stack <- rast(precip_files)
precip_stack
names(precip_stack) <- c('precip_01', 'precip_02', 'precip_03', 'precip_04', 'precip_05', 'precip_06', 'precip_07', 'precip_08', 'precip_09', 'precip_10', 'precip_11', 'precip_12')


## (3) load PET
setwd(pet_wd)
pet_files <- list.files(pattern = '.tif')
pet_files
pet_stack <- rast(pet_files)
pet_stack
# plot(pet_stack)
names(pet_stack) <- c('pet_01', 'pet_02', 'pet_03', 'pet_04', 'pet_05', 'pet_06', 'pet_07', 'pet_08', 'pet_09', 'pet_10', 'pet_11', 'pet_12')


## (4) load SOCS in 1986
socs <- rast("W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version1_20250418/02_data_spinup/01_initial_carbon_stock/Initial_SOC_stock0_30cm.tif")
names(socs) <- 'socs'


## (5) load clay 
clay <- rast("W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version1_20250418/02_data_spinup/09_clay_anatol_0_30cm/clay_0_30cm_aggregated.tif")
names(clay) <- 'clay'


## (6) load ndvi 
setwd(ndvi_wd)
ndvi_files <- list.files(pattern = '.tif')
ndvi_files
ndvi_stack <- rast(ndvi_files)
ndvi_stack
names(ndvi_stack) <- c('ndvi01', 'ndvi02', 'ndvi03', 'ndvi04', 'ndvi05', 'ndvi06', 'ndvi07', 'ndvi08', 'ndvi09', 'ndvi10','ndvi11', 'ndvi12')

## (7) load DR 
DR <- rast('W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version1_20250418/02_data_spinup/07_DR/DR_1986.tif')
DR1 <- DR/100
names(DR1) <- 'DR'

# plot(DR1)
## (8) load crop residues/manure ratio
ratio_crop_manure <- rast('W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/DATA/03_Carbon_input/Manure_data_new_20251002/ratio_crop_manure.tif')

manure_dr <- rast("W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/DATA/03_Carbon_input/Manure_data_new_20251002/weighted_hc25m.tif")


## stack all variables
covariates_stack <- c(clay, socs, temp_stack, precip_stack, pet_stack, ndvi_stack, DR1, ratio_crop_manure, manure_dr)
crs(covariates_stack) <- 'epsg:28992'


st_crs(sf_point1) <- 'epsg:28992'

# extract raster values to points
Vector_variables <- terra::extract(covariates_stack, vect(sf_point1), df=TRUE)
names(Vector_variables)
head(Vector_variables)
# coords <- crds(sf_point)
coords <- st_coordinates(sf_point1)
## combine coordinates (X, Y) into the variable table
Vector_variables <- cbind(Vector_variables, coords)
names(Vector_variables)


#### -------------------03_RothC model running---------------------------

SOC_im <- Vector_variables[[3]] 

clay_im <- Vector_variables[[2]] 

DR_im <- Vector_variables[[52]]
cm_ratio_im <- Vector_variables[[53]]
manure_dr_im <- Vector_variables[[54]]

# define the temporary saving function
# here you shold change to your own saving path..
save_rothc <- function(df, iteration) {
  rothc_file <- paste0("W:/ESG/DOW_SGL/Research_PhD/YuqingLai/Paper1/RothC_new_processing/version2_20251010/spinup/point2/GH_analy_result_iter1v_", iteration, ".csv")
  write.csv(df, rothc_file, row.names = FALSE)
  message("Saved: ", rothc_file)
}

# set the number of samples to be saved at a time 
save_interval <- 30000

# 4. RothC parameter setting
ks <- c(10, 0.3, 0.66, 0.02, 0)


# initialize one dataframe to store all results 
all_results1 <- data.frame()

# loop through all points
for (i in 1:nrow(sf_point1)) {
  
  ## clean..
  message(i, 'point')
  if (i %% 1000 == 0) gc()
  
  # Extract the variables 
  Vect<-as.data.frame(Vector_variables[i,])
  
  # temperature
  Temp <- data.frame(t(Vect[,4:15]))
  colnames(Temp) <- 'temp'
  # Temp <- data.frame(Month=1:12, Temp=Temp[,1])
  
  # precipitation
  Precip <- as.data.frame(t(Vect[,16:27]))
  Precip <- data.frame(Month=1:12, Precip=Precip[,1])
  
  # pet
  Evp<-as.data.frame(t(Vect[,28:39]))
  Evp<-data.frame(Month=1:12, Evp=Evp[,1])
  
  # ndvi
  Cov<-as.data.frame(t(Vect[40:51]))
  Cov1<-data.frame(Cov=Cov[,1])
  Cov2<-data.frame(Month=1:12, Cov=Cov[,1])
  
  
  soil_thick=30  #Soil thickness (organic layer topsoil), in cm
  
  CTOT <- SOC_im[i] #Soil organic carbon in ton/ha 
  clay<-clay_im[i] #Percent clay %
  
  DR<-DR_im[i] # DPM/RPM (decomplosable vs resistant plant material.)
  bare<-(Cov1<0.2) # If the surface is bare or vegetated
  CRM_ratio <-  cm_ratio_im[i]  # ratio of crop residues to manure
  mr <- manure_dr_im [i]
  #calculate IOM using Falloon method
  FallIOM=0.049*CTOT^(1.139) 
  
  # ratio of BIO and HUM 
  x <- 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
  
  # calculate transfer coefficient
  B = 0.46/(x + 1)
  H = 0.54/(x + 1)
  
  ai3 = B * ks
  ai4 = H * ks
  A = diag(-ks)
  A[3, ] = A[3, ] + ai3
  A[4, ] = A[4, ] + ai4
  
  
  # avoid NA..
  if (any(is.na(Evp[,2])) | any(is.na(Temp)) | any(is.na(CTOT)) | 
      any(is.na(clay)) | any(is.na(Precip[,2])) | any(is.na(Cov2[,2])) | 
      any(is.na(DR)) | (CTOT < 0) | (clay < 0)) {
    output_row <- data.frame(
      ID = Vect$ID, X = Vect$X, Y = Vect$Y,
      DPM = NA, RPM = NA, BIO = NA, HUM = NA, IOM = NA, 
      totalC = NA, Cin = NA,CRP = NA,
      MP = NA
    )
  } else {
    
    # apply function to calculate RateM
    xi1 <- calculate_rate_modifiers(Temp = Temp, Precip = Precip, Evp = Evp, 
                                    soil_thick = soil_thick, clay = clay, 
                                    bare = bare, Cov = Cov)
    xi=mean(xi1$RateM)
    
    ## other parameters of RothC
    # DPM/RPM ratio
    rho <- DR / (1 + DR)
    rho
    
    tau <- mr / (1 + mr)  
    # nu <- 0.49
    
    # calculate carbon input ratio
    CR_proportion <- CRM_ratio / (CRM_ratio + 1)
    M_proportion <- 1 / (CRM_ratio + 1)
    
    # C_active and input vector of RothC model
    CTOT4 <- CTOT - FallIOM  # active carbon pool
    rho_vector <- matrix(c(rho, 1 - rho, 0, 0), nrow = 4, ncol = 1)
    tau_nu_vector <- matrix(c(0.98 * tau, 0.98*(1-tau), 0, 0.02), nrow = 4, ncol = 1)
    input_vector <- CR_proportion * rho_vector + M_proportion * tau_nu_vector
    
    # calculate carbon input
    I4 <- matrix(c(1, 1, 1, 1), nrow = 1, ncol = 4)
    A4 <- A[1:4, 1:4]
    A4inv <- solve(A4)
    
    I <- as.numeric(CTOT4 / (-(1 / xi) * I4 %*% A4inv %*% input_vector))
    
    # carbon pool calculation
    Q <- input_vector * I
    C <- -(1 / xi) * A4inv %*% Q
    
    # calculate carbon pool at equilibrium state
    Cpools <- c(C, FallIOM)
    Ctotal <- sum(Cpools)
    
    output_row <- data.frame(
      ID = Vect$ID, 
      X = Vect$X, 
      Y = Vect$Y,
      DPM = round(Cpools[1],3), 
      RPM = round(Cpools[2],3), 
      BIO = round(Cpools[3],3), 
      HUM = round(Cpools[4],3), 
      IOM = round(Cpools[5],3), 
      totalC = round(Ctotal,3), 
      Cin = round(I,3), 
      CRP = round(CR_proportion,3),
      MP = round(M_proportion,3)
    )
  }
  
  # add result to all_results
  all_results1 <- rbind(all_results1, output_row)
  
  # save result every 10000 samples
  if (i %% save_interval == 0) {
    save_rothc(all_results1, iteration = i)
    all_results1 <- data.frame()  
  }
}

# save final result after loop
if (nrow(all_results1) > 0) {
  save_rothc(all_results1, iteration = i)
}

