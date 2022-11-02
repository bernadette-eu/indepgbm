###########################################################################
#
# Paths
#
###########################################################################
path_to_source <- ".../R"

###########################################################################
#
# Load libraries:
#
###########################################################################
lib <- c("ggplot2",
         "tidyverse",
         "tidyr",
         "rvest",
         "vroom",
         "readxl",
         "readr", 
         "osfr",
         "rio",
         "httr",
         "remotes", 
         'lubridate',
         "extraDistr",
         "gridExtra",
         "mgcv",
         "ggrepel",
         "patchwork",
         "rstan",
         "bayesplot",
         "Bernadette",
         "covidAgeData",
	       "loo"
)
lapply(lib, require, character.only = TRUE)

lapply(lib,
       FUN = function(x) {
         if("covidAgeData" %in% lib) {
           remotes::install_github("eshom/covid-age-data")
         } else install.packages(x, dependencies = TRUE)
         
         library(x, character.only = TRUE)
       }# End function
)

###########################################################################
#
# Source files:
#
###########################################################################

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

source(paste0(path_to_source, "/", "2_Data_Engineering", ".R"))
source(paste0(path_to_source, "/", "3_rt_ngm",           ".R"))
source(paste0(path_to_source, "/", "4_stan-utility",     ".R"))
