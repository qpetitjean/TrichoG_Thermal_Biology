######################################################################################################################
# Script Title: Retrieve and Summarize Climatic Data for Trichogramma Sampling Sites
# 
# Author: Quentin PETITJEAN
# Date Created: 22/05/2025
# Last Modified: 22/05/2025
# ============================================================================== 
# Requirements:
# - R version 4.2.3
# - Packages:
#   - easyclimate v0.2.2: For accessing daily historical climate data.
#   - data.table v1.14.8: For efficient data manipulation and file writing.
#   - progress v1.2.2: For displaying progress during data processing.
# ============================================================================== 
# Script Overview:
# This script retrieves and summarizes climatic data for each Trichogramma strain's sampling location using 
# the `easyclimate` R package. For each strain, daily minimum and maximum temperature data are extracted over 
# a 5-year period preceding the year of field sampling. The script then computes synthetic climatic indices:
#   - Number of days with Tmin < 0°C
#   - Number of days with Tmax > 25°C
#   - Average annual Tmax and Tmin across 5 years
#   - Mean, max, and min daily temperature range across 5 years
# These values are stored in a summary table, while the complete daily temperature records are saved separately.
# ============================================================================== 
# Usage:
# 1. Set the 'Path' variable to point to the Zenodo data repository.
# 2. Ensure the `FullDatasetBidime.csv` file contains the required fields: Strain, Lon, Lat, Annee.capture.
# 3. Run the script to generate:
#    - RDS file: `RawClimateData.Rds`, containing full daily temperature records per strain.
#    - CSV file: `climateIndices.csv`, containing summarized climatic indices for each strain.
# ============================================================================== 


######################################################################################################################
# Install and attach needed packages
######################################################################################################################

if (!require(easyclimate)) {
  install.packages("easyclimate")
}
if (!require(progress)) {
  install.packages("progress")
}

########################################################################################################################
# set parameters (only change the Path variable)
#
########################################################################################################################

# the path to the directory downloaded from zenodo repository
Path <- "W:/Postdoc_INRAE_SAM/Zenodo_Data_Repo"

# the path where results are stored
GlobSavingDir <- file.path(Path, "ProcessedDataset")

######################################################################################################################
# Import the dataset 
######################################################################################################################

FullDat <- data.table::fread(file.path(GlobSavingDir, "FullDatasetBidime.csv"), sep = ";", dec = ".")

######################################################################################################################
# Retrieve climatic data using easyClimate R package
######################################################################################################################

RawClimateData <- list()
climateIndices <-
  stats::setNames(
    data.frame((matrix(
      nrow = nrow(FullDat), ncol = 8
    ))),
    c("strain",
      "nDayBelow0",
      "nDayAbove25",
      "Tmax",
      "Tmin",
      "meanDailyTrange",
      "maxDailyTrange",
      "minDailyTrange"
    )
  )

# initialize progress bar
total <- nrow(FullDat)
pb <-
  progress::progress_bar$new(format = "Strain Processing [:bar] :current/:total (:percent)", total = total)
pb$tick(0)

for(i in 1:nrow(FullDat)){
  
  # retrieve spatial coordinates of the sampling point
  coords <- data.frame(lon = FullDat$Lon[i], lat = FullDat$Lat[i])
  captureYear <- FullDat$Annee.capture[i]
  if(is.na(captureYear)) captureYear = 2000
  
  # Set time range as 5 years prior sampling has been conducted in the field
  start <- paste((captureYear - 5), "01", "01", sep = "-")
  end <- paste((captureYear - 1), "12", "31", sep = "-")
  period <- paste(start, end, sep = ":")
  
  # Retrieve daily temperature data (min and max)
  dailyTemp <- easyclimate::get_daily_climate(
    coords = coords,
    period = period,
    climatic_var = c("Tmin", "Tmax"),
    version = 4
  )
  
  # compute mean temperature and range
  dailyTemp$Tmean <- (dailyTemp$Tmin + dailyTemp$Tmax) / 2
  dailyTemp$Trange <- dailyTemp$Tmax - dailyTemp$Tmin
  
  # append the strain name
  dailyTemp$strain <- FullDat$Strain[i]
  
  # Store the raw climatic data in a list
  RawClimateData[[FullDat$Strain[i]]] <- dailyTemp
  
  # Fill the climateIndices dataframe by computing synthetic climatic values
  dailyTemp$year <- format(as.Date(dailyTemp$date), "%Y")
  
  ## number of days below 0 degrees C
  nDayBelow0 <- length(which(dailyTemp$Tmin < 0))
  
  ## number of days above 25 degrees C
  nDayAbove25 <- length(which(dailyTemp$Tmax > 25))
  
  ## the average max temperature reached over 5 year
  Tmax <- mean(stats::aggregate(
    Tmax ~ year,
    data = dailyTemp,
    FUN = function(x)
      max(x, na.rm = T)
  )$Tmax)
  
  ## the average min temperature reached over 5 year
  Tmin <- mean(stats::aggregate(
    Tmin ~ year,
    data = dailyTemp,
    FUN = function(x)
      min(x, na.rm = T)
  )$Tmin)
    
  ## the mean daily temperature range over 5 year
  meanDailyTrange <- mean(dailyTemp$Trange, na.rm = T)
  
  ## the max daily temperature range over 5 year
  maxDailyTrange <- max(dailyTemp$Trange, na.rm = T)
  
  ## the min daily temperature range over 5 year
  minDailyTrange <- min(dailyTemp$Trange, na.rm = T)
  
  ## fill the dataframe
  climateIndices[i, ] <-
    c(FullDat$Strain[i],
      nDayBelow0,
      nDayAbove25,
      Tmax,
      Tmin,
      meanDailyTrange,
      maxDailyTrange,
      minDailyTrange)
  
  # update progress bar
  pb$tick(1)
}

# save the raw daily temperature as .rds
saveRDS(RawClimateData, file = file.path(GlobSavingDir, "RawClimateData.Rds"))

# save the climatic indices data
data.table::fwrite(
  climateIndices, file.path(GlobSavingDir,"climateIndices.csv"),
  sep = ";",
  dec = ".",
  na = "NA"
)
