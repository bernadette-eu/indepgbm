###########################################################################
#
# Misc:
#
###########################################################################
`%nin%` <- Negate(`%in%`)

grid_arrange_shared_legend <- function(plots, nrow) {
  
  g       <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend  <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  
  plots_no_legend <- lapply(plots, function(x) x + theme(legend.position="none"))
  
  gridExtra::grid.arrange(
    gridExtra::arrangeGrob(grobs = plots_no_legend, nrow = nrow),
    legend,
    nrow    = 2,
    heights = grid::unit.c(unit(1, "npc") - lheight, lheight)
  )
}# End function  

# Convert prior mean and variance to alpha and beta parameters for a Beta prior:
# NOTE: These calculations will only work (alpha, beta > 0 ) if the variance is less than the mean*(1-mean).
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = c(alpha = alpha, beta = beta))
}# End function

glyph_vline <- function(data, params, size) {
  if (data$linetype == 4) {
    draw_key_vline(data, params, size)
  } else {
    zeroGrob()
  }
}# End function

###########################################################################
#
#---- Data engineering of mortality and infection counts:
#
###########################################################################
data_engineering_observations_GR <- function(start_date  = NULL,
                                             end_date    = NULL,
                                             age_mapping_deaths
){
  
  if (is.null(start_date) & 
      is.null(end_date)) message("Time series start and end date are not specified.\n The whole time series will be considered.")
  
  repo <- "https://raw.githubusercontent.com/Sandbird/covid19-Greece/master/"
  
  data <- readr::read_csv( paste0(repo, "demography_total_details", ".csv"), col_types = readr::cols() )
  
  # The data start from 2020-04-03
  data <- subset(data, date >= "2020-04-03")
  
  #---- Calculate the new deaths and cases per age group:
  data_groups_cols             <- c("Date", "Group", "Cumul_cases", "Cumul_deaths")
  dt_analysis_period           <- data.frame(matrix(ncol = length(data_groups_cols), nrow = 0))
  colnames(dt_analysis_period) <- data_groups_cols
  
  #---- Rename the columns according to the mapping:
  lookup_table <- data.frame(Initial = unique(data$category),
                             Mapping = age_mapping_deaths)
  age_bands    <- unique(unique(lookup_table$Mapping))
  age_bands_no <- length(age_bands)
  dates        <- unique(data$date)
  
  data_melt <-data %>% 
    dplyr::left_join(lookup_table, 
                     by = c(category = "Initial")) %>% 
    dplyr::select(-dplyr::one_of(c("category"))) %>% 
    dplyr::group_by(date, Mapping) %>% 
    dplyr::summarise(Cumul_cases  = sum(cases),
                     Cumul_deaths = sum(deaths)) %>% 
    dplyr::rename(Group = Mapping)
  
  for (i in 1:age_bands_no){
    
    data_subset <- data_melt[data_melt$Group == age_bands[i],]
    
    dt_age_grp <- data.frame(Date         = dates,
                             Group        = rep(age_bands[i], length(dates)),
                             Cumul_cases  = c(data_subset[data_subset$Group == age_bands[i],"Cumul_cases"]),
                             Cumul_deaths = c(data_subset[data_subset$Group == age_bands[i],"Cumul_deaths"])
    ) 
    
    # Calculate new deaths and cases from the cumulative deaths and cases, respectively:
    new_cases_age_grp    <- rep(0, nrow(dt_age_grp))
    new_cases_age_grp[1] <- data_subset$Cumul_cases[1]
    
    new_deaths_age_grp    <- rep(0, nrow(dt_age_grp))
    new_deaths_age_grp[1] <- data_subset$Cumul_deaths[1]
    
    for (t in 2:length(dates)) {
      new_cases_age_grp[t]  <- data_subset$Cumul_cases[t]  - data_subset$Cumul_cases[t-1]
      new_deaths_age_grp[t] <- data_subset$Cumul_deaths[t] - data_subset$Cumul_deaths[t-1]
      
      # Sanity check - if the cumulative cases/deaths are lower than the day before, set to zero.
      # Confirmed by EODY reports for the dates 2020-04-09 (65+), 2020-05-05 (40-64), 2020-05-21 (40-64), 2020-05-27 (65+)
      if(new_cases_age_grp[t] < 0)  new_cases_age_grp[t]  <- 0
      if(new_deaths_age_grp[t] < 0) new_deaths_age_grp[t] <- 0
      
    }# End for
    
    dt_age_grp$New_Cases  <- new_cases_age_grp
    dt_age_grp$New_Deaths <- new_deaths_age_grp
    
    dt_analysis_period <- rbind(dt_analysis_period, dt_age_grp)  
    
  }# End for  
  
  #---- Subset for the analysis period:
  dt_analysis_period <- dt_analysis_period[ dt_analysis_period$Date >= as.Date(start_date, format = "%Y-%m-%d") & 
                                              dt_analysis_period$Date <= as.Date(end_date, format = "%Y-%m-%d"),]
  
  #---- Mortality data to be processed by the .stan file:
  dt_mortality_analysis_period <- reshape(dt_analysis_period[,c(1,2,ncol(dt_analysis_period))],
                                          idvar     = "Date",
                                          timevar   = "Group",
                                          direction = "wide")
  colnames(dt_mortality_analysis_period) <- gsub("New_Deaths.", "",
                                                 colnames(dt_mortality_analysis_period))
  
  dt_mortality_analysis_period <- dt_mortality_analysis_period %>% 
    dplyr::select(-dplyr::one_of(c("Date"))) %>% 
    dplyr::mutate(New_Deaths = rowSums(.))
  
  #---- Incidence data to be processed by the .stan file:
  dt_cases_analysis_period <- reshape(dt_analysis_period[,c(1,2,ncol(dt_analysis_period)-1)],
                                      idvar     = "Date",
                                      timevar   = "Group",
                                      direction = "wide")
  colnames(dt_cases_analysis_period) <- gsub("New_Cases.", "",
                                             colnames(dt_cases_analysis_period))
  
  dt_cases_analysis_period <- dt_cases_analysis_period %>% 
    dplyr::select(-dplyr::one_of(c("Date"))) %>% 
    dplyr::mutate(New_Cases = rowSums(.))
  
  # Additional fields:
  date_subset <- unique(dt_analysis_period$Date)
  
  dt_mortality_analysis_period$Date  <- date_subset
  dt_cases_analysis_period$Date      <- dt_mortality_analysis_period$Date
  
  dt_mortality_analysis_period$Index <- 1:nrow(dt_mortality_analysis_period)
  dt_cases_analysis_period$Index     <- 1:nrow(dt_cases_analysis_period)
  
  #---- Add Week ID since first day of analysis:
  # Search: Stackoverflow "How can I group days into weeks"
  dt_mortality_analysis_period$Week_ID <- as.numeric(dt_mortality_analysis_period$Date - dt_mortality_analysis_period$Date[1]) %/% 7 + 1
  dt_cases_analysis_period$Week_ID     <- as.numeric(dt_cases_analysis_period$Date - dt_cases_analysis_period$Date[1]) %/% 7 + 1
  
  #---- Left and right time limit:
  dt_mortality_analysis_period$Right <- dt_mortality_analysis_period$Index + 1
  dt_cases_analysis_period$Right     <- dt_cases_analysis_period$Index + 1
  
  #---- Re-order the column names"
  fixed_cols        <- c("Index", "Right", "Date", "Week_ID")
  fixed_cols_cases  <- c(fixed_cols, "New_Cases")
  fixed_cols_deaths <- c(fixed_cols, "New_Deaths")
  
  reorder_colnames <- c(fixed_cols_cases, colnames(dt_cases_analysis_period)[!colnames(dt_cases_analysis_period) %in% fixed_cols_cases])
  dt_cases_analysis_period <- dt_cases_analysis_period[reorder_colnames]
  
  reorder_colnames <- c(fixed_cols_deaths, colnames(dt_mortality_analysis_period)[!colnames(dt_mortality_analysis_period) %in% fixed_cols_deaths])
  dt_mortality_analysis_period <- dt_mortality_analysis_period[reorder_colnames]
  
  #--- Table at level (Date-Group) with the cumulative cases (overall and dissaggregated) for the analysis period:
  data_cumcases_wide <- as.data.frame(data_melt[,c(1,2,3)])
  data_cumcases_wide <- reshape(data_cumcases_wide,
                                idvar     = "date",
                                timevar   = "Group",
                                direction = "wide")
  colnames(data_cumcases_wide) <- gsub("Cumul_cases.", "",
                                       colnames(data_cumcases_wide))
  colnames(data_cumcases_wide)[1] <- "Date"
  
  data_cumcases_wide <- data_cumcases_wide %>% 
    dplyr::select(-dplyr::one_of(c("Date"))) %>% 
    dplyr::mutate(Total_Cases = rowSums(.))
  
  data_cumcases_wide$Date  <- dates
  
  fixed_cols_cumcases <- c("Date", "Total_Cases")
  reorder_colnames_cumcases <- c(fixed_cols_cumcases, colnames(data_cumcases_wide)[!colnames(data_cumcases_wide) %in% fixed_cols_cumcases])
  
  data_cumcases_wide <- data_cumcases_wide[reorder_colnames_cumcases]
  
  data_cumcases_wide <- data_cumcases_wide[ data_cumcases_wide$Date >= as.Date(start_date, format = "%Y-%m-%d") & 
                                              data_cumcases_wide$Date <= as.Date(end_date,   format = "%Y-%m-%d"),]
  
  output_list <- list(Deaths    = dt_mortality_analysis_period,
                      Cases     = dt_cases_analysis_period,
                      All       = dt_analysis_period,
                      Cum_Cases = data_cumcases_wide)
  
  return(output_list)
  
}# End function

data_engineering_observations_AT <- function(start_date = NULL,
                                             end_date   = NULL,
                                             age_mapping_deaths,
                                             data_cases,
                                             age_mapping_cases){
  
  if (is.null(start_date) & 
      is.null(end_date)) message("Time series start and end date are not specified.\n The whole time series will be considered.")
  
  packages <- c("rio")
  
  lapply(packages,
         FUN = function(x) {
           if (!require(x, character.only = TRUE)) {
             install.packages(x, dependencies = TRUE)
             library(x, character.only = TRUE)
           }
         }
  )
  
  #---- Import deaths:
  data_deaths_url <- "https://www.ined.fr/fichier/rte/166/Page%20Data/Austria/Austria_2022_04_21_Deaths_Covid_19.xlsx"
  data_deaths     <- rio::import(data_deaths_url, which = 4)
  
  #---- Remove redundant rows:
  data_deaths <- data_deaths[1:16,]
  
  # Store the dates:
  store_dates <- data_deaths[5,]
  store_dates <- store_dates[!is.na(store_dates)]
  
  dates_exclude <- c("Age Group", "Population on 01/01/2019",
                     store_dates[grepl(".08.", store_dates, fixed = TRUE)])
  
  store_dates <- store_dates[store_dates %nin% dates_exclude]
  fixed_dates <- as.Date(as.numeric(store_dates), origin = "1899-12-30")
  
  data_deaths <- data_deaths[-c(1:5),]
  
  colnames(data_deaths)    <- data_deaths[1,]
  colnames(data_deaths)[1] <- "Age"
  data_deaths <- data_deaths[-1,]
  
  data_deaths <- data_deaths[,-c(2:7)]
  
  data_deaths <- data_deaths[c("Age", 
                               colnames(data_deaths)[grepl("Both", colnames(data_deaths), fixed = TRUE)])]
  
  map_dates    <- colnames(data_deaths)[-1]
  map_dates_df <- data.frame(ID      = map_dates,
                             Mapping = fixed_dates)
  
  data_deaths_long <- tidyr::pivot_longer(data_deaths, 
                                          -c(Age), 
                                          values_to = "Deaths", 
                                          names_to  = "Date") 
  
  data_deaths_long <- data_deaths_long %>% 
    as.data.frame() %>% 
    dplyr::left_join(map_dates_df, 
                     by = c(Date = "ID")) %>% 
    dplyr::select(-dplyr::one_of(c("Date"))) %>% 
    dplyr::rename(Date = Mapping,
                  Group = Age)
  
  data_cases_long <- data_cases
  
  #---- Calculate the new deaths and cases per age group:
  data_groups_cols             <- c("Date", 
                                    "Group", 
                                    "Cumul_cases", 
                                    "Cumul_deaths")
  dt_analysis_period           <- data.frame(matrix(ncol = length(data_groups_cols), nrow = 0))
  colnames(dt_analysis_period) <- data_groups_cols
  
  #---- Rename the columns according to the mapping:
  lookup_table <- data.frame(Initial = unique(data_deaths_long$Group),
                             Mapping = age_mapping_deaths)
  age_bands    <- unique(unique(lookup_table$Mapping))
  age_bands_no <- length(age_bands)
  dates        <- unique(data_deaths_long$Date)
  dates        <- dates[order(dates , decreasing = FALSE )]
  
  data_deaths_long$Deaths <- as.numeric(data_deaths_long$Deaths)
  
  data_deaths_melt <- data_deaths_long %>% 
    dplyr::left_join(lookup_table, 
                     by = c(Group = "Initial")) %>% 
    dplyr::select(-dplyr::one_of(c("Group"))) %>% 
    dplyr::group_by(Date, Mapping) %>% 
    dplyr::summarise(Cumul_deaths = sum(Deaths)) %>% 
    dplyr::rename(Group = Mapping)
  
  lookup_table_cases <- data.frame(Initial = unique(data_cases_long$Group),
                                   Mapping = age_mapping_cases)
  
  data_cases_melt <-data_cases_long %>% 
    dplyr::left_join(lookup_table_cases, 
                     by = c(Group = "Initial")) %>% 
    dplyr::select(-dplyr::one_of(c("Group"))) %>% 
    dplyr::group_by(Date, Mapping) %>% 
    dplyr::summarise(Cumul_cases = sum(Cases)) %>% 
    dplyr::rename(Group = Mapping) %>% 
    dplyr::filter(Date >= min(dates) & Date <= max(dates)) %>% 
    dplyr::mutate(Cumul_cases = round(Cumul_cases,0))
  
  for (i in 1:age_bands_no){
    
    data_deaths_subset <- data_deaths_melt[data_deaths_melt$Group == age_bands[i],]
    data_cases_subset  <- data_cases_melt[data_cases_melt$Group   == age_bands[i],]
    
    dt_age_grp <- data.frame(Date         = dates,
                             Group        = rep(age_bands[i], length(dates)),
                             Cumul_cases  = c(data_cases_subset[data_cases_subset$Group == age_bands[i],"Cumul_cases"]),
                             Cumul_deaths = c(data_deaths_subset[data_deaths_subset$Group == age_bands[i], "Cumul_deaths"])
    ) 
    
    # Calculate new deaths and cases from the cumulative deaths and cases, respectively:
    new_cases_age_grp    <- rep(0, nrow(dt_age_grp))
    new_cases_age_grp[1] <- data_cases_subset$Cumul_cases[1]
    
    new_deaths_age_grp    <- rep(0, nrow(dt_age_grp))
    new_deaths_age_grp[1] <- data_deaths_subset$Cumul_deaths[1]
    
    for (t in 2:length(dates)) {
      
      new_cases_age_grp[t]  <- data_cases_subset$Cumul_cases[t]   - data_cases_subset$Cumul_cases[t-1]
      new_deaths_age_grp[t] <- data_deaths_subset$Cumul_deaths[t] - data_deaths_subset$Cumul_deaths[t-1]
      
      # Sanity check - if the cumulative cases/deaths are lower than the day before, set to zero.
      # Confirmed by EODY reports for the dates 2020-04-09 (65+), 2020-05-05 (40-64), 2020-05-21 (40-64), 2020-05-27 (65+)
      
      if(new_cases_age_grp[t] < 0)  new_cases_age_grp[t]  <- 0
      if(new_deaths_age_grp[t] < 0) new_deaths_age_grp[t] <- 0
      
    }# End for
    
    dt_age_grp$New_Cases  <- new_cases_age_grp
    dt_age_grp$New_Deaths <- new_deaths_age_grp
    
    dt_analysis_period <- rbind(dt_analysis_period, dt_age_grp)  
    
  }# End for  
  
  #---- Subset for the analysis period:
  dt_analysis_period <- dt_analysis_period[ dt_analysis_period$Date >= as.Date(start_date, format = "%Y-%m-%d") & 
                                              dt_analysis_period$Date <= as.Date(end_date, format = "%Y-%m-%d"),]
  
  #---- Mortality data to be processed by the .stan file:
  dt_mortality_analysis_period <- reshape(dt_analysis_period[,c(1,2,ncol(dt_analysis_period))],
                                          idvar     = "Date",
                                          timevar   = "Group",
                                          direction = "wide")
  colnames(dt_mortality_analysis_period) <- gsub("New_Deaths.", "",
                                                 colnames(dt_mortality_analysis_period))
  
  dt_mortality_analysis_period <- dt_mortality_analysis_period %>% 
    dplyr::select(-dplyr::one_of(c("Date"))) %>% 
    dplyr::mutate(New_Deaths = rowSums(.))
  
  #---- Incidence data to be processed by the .stan file:
  dt_cases_analysis_period <- reshape(dt_analysis_period[,c(1,2,ncol(dt_analysis_period)-1)],
                                      idvar     = "Date",
                                      timevar   = "Group",
                                      direction = "wide")
  colnames(dt_cases_analysis_period) <- gsub("New_Cases.", "",
                                             colnames(dt_cases_analysis_period))
  
  dt_cases_analysis_period <- dt_cases_analysis_period %>% 
    dplyr::select(-dplyr::one_of(c("Date"))) %>% 
    dplyr::mutate(New_Cases = rowSums(.))
  
  # Additional fields:
  date_subset <- unique(dt_analysis_period$Date)
  
  dt_mortality_analysis_period$Date  <- date_subset
  dt_cases_analysis_period$Date      <- dt_mortality_analysis_period$Date
  
  dt_mortality_analysis_period$Index <- 1:nrow(dt_mortality_analysis_period)
  dt_cases_analysis_period$Index     <- 1:nrow(dt_cases_analysis_period)
  
  #---- Add Week ID since first day of analysis:
  dt_mortality_analysis_period$Week_ID <- as.numeric(dt_mortality_analysis_period$Date - dt_mortality_analysis_period$Date[1]) %/% 7 + 1
  dt_cases_analysis_period$Week_ID     <- as.numeric(dt_cases_analysis_period$Date - dt_cases_analysis_period$Date[1]) %/% 7 + 1
  
  #---- Left and right time limit:
  dt_mortality_analysis_period$Right <- dt_mortality_analysis_period$Index + 1
  dt_cases_analysis_period$Right     <- dt_cases_analysis_period$Index + 1
  
  #---- Re-order the column names"
  fixed_cols        <- c("Index", "Right", "Date", "Week_ID")
  fixed_cols_cases  <- c(fixed_cols, "New_Cases")
  fixed_cols_deaths <- c(fixed_cols, "New_Deaths")
  
  reorder_colnames <- c(fixed_cols_cases, colnames(dt_cases_analysis_period)[!colnames(dt_cases_analysis_period) %in% fixed_cols_cases])
  dt_cases_analysis_period <- dt_cases_analysis_period[reorder_colnames]
  
  reorder_colnames <- c(fixed_cols_deaths, colnames(dt_mortality_analysis_period)[!colnames(dt_mortality_analysis_period) %in% fixed_cols_deaths])
  dt_mortality_analysis_period <- dt_mortality_analysis_period[reorder_colnames]
  
  #--- Table at level (Date-Group) with the cumulative cases (overall and dissaggregated) for the analysis period:
  data_cumcases_wide <- as.data.frame(data_cases_melt)
  data_cumcases_wide <- reshape(data_cumcases_wide,
                                idvar     = "Date",
                                timevar   = "Group",
                                direction = "wide")
  colnames(data_cumcases_wide) <- gsub("Cumul_cases.", "",
                                       colnames(data_cumcases_wide))
  
  data_cumcases_wide <- data_cumcases_wide[ data_cumcases_wide$Date >= as.Date(start_date, format = "%Y-%m-%d") & 
                                              data_cumcases_wide$Date <= as.Date(end_date,   format = "%Y-%m-%d"),]
  
  data_cumcases_wide <- data_cumcases_wide %>% 
    dplyr::select(-dplyr::one_of(c("Date"))) %>% 
    dplyr::mutate(Total_Cases = rowSums(.))
  
  data_cumcases_wide$Date  <- date_subset
  
  fixed_cols_cumcases <- c("Date", "Total_Cases")
  reorder_colnames_cumcases <- c(fixed_cols_cumcases, colnames(data_cumcases_wide)[!colnames(data_cumcases_wide) %in% fixed_cols_cumcases])
  
  data_cumcases_wide <- data_cumcases_wide[reorder_colnames_cumcases]
  
  output_list <- list(Deaths    = dt_mortality_analysis_period,
                      Cases     = dt_cases_analysis_period,
                      All       = dt_analysis_period,
                      Cum_Cases = data_cumcases_wide)
  
  return(output_list)
  
}# End function


data_engineering_observations_ENG <- function(start_date = NULL,
                                              end_date   = NULL,
                                              age_mapping_deaths,
                                              age_mapping_cases){
  
  if (is.null(start_date) & 
      is.null(end_date)) message("Time series start and end date are not specified.\n The whole time series will be considered.")
  
  packages <- c("rio", "httr")#), "remotes", "osfr", "covidAgeData")
  
  lapply(packages,
         FUN = function(x) {
           if (!require(x, character.only = TRUE)) {
             install.packages(x, dependencies = TRUE)
             library(x, character.only = TRUE)
           }
         }
  )
  
  #---- Import deaths:
  data_deaths_url <- "https://www.ined.fr/fichier/rte/166/Page%20Data/England%20and%20Wales/EnglandWales_2022_04_22_Deaths_Covid-19.xlsx"
  data_deaths     <- rio::import(data_deaths_url, which = 6)
  
  #---- Import cases:
  urlcontent      <- httr::GET("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDateAgeDemographics&format=csv&release=2022-05-06")
  data_cases      <- read.csv(text=content(urlcontent, type = "text", encoding="UTF-8"))
  data_cases      <- data_cases[c("date", "age", "cases" )]
  data_cases_long <- subset(data_cases, age %nin% c("00_59", "60+", "unassigned"))
  
  #---- Age distribution (based on the Office for National Statistics, 2020):
  url_ons <- "https://www.ons.gov.uk/file?uri=%2fpeoplepopulationandcommunity%2fpopulationandmigration%2fpopulationestimates%2fdatasets%2fpopulationestimatesforukenglandandwalesscotlandandnorthernireland%2fmid2020/ukpopestimatesmid2020on2021geography.xls"
  age_distribution <- rio::import(url_ons, which = 6)
  age_distribution <- age_distribution[-c(1:11),c(1,5)]
  
  colnames(age_distribution) <- c("AgeGrp", "PopTotal")
  age_distribution$PopTotal  <- as.numeric(age_distribution$PopTotal)
  
  lookup_table_cases <- data.frame(ID      = unique(age_distribution$AgeGrp),
                                   Mapping = age_mapping_cases)
  
  age_distribution_aggr <-age_distribution %>% 
                          dplyr::left_join(lookup_table_cases, 
                                           by = c(AgeGrp = "ID")) %>% 
                          dplyr::select(-dplyr::one_of(c("AgeGrp"))) %>% 
                          dplyr::group_by(Mapping) %>% 
                          dplyr::summarise(PopTotal = sum(PopTotal)) %>% 
                          dplyr::rename(AgeGrp = Mapping)
  
  age_distribution$AgeGrp <- gsub(" to ", "-", age_distribution$AgeGrp)
  age_distribution$AgeGrp <- gsub(" and over", "+", age_distribution$AgeGrp)

  #---- Remove redundant rows:
  data_deaths <- data_deaths[1:11,]
  
  # Store the dates - Option 1:
  store_dates <- data_deaths[5,]
  store_dates <- store_dates[!is.na(store_dates)]
  
  dates_exclude <- c("Age Group", 
                     "Population on 30/06/2018",
                     "Total", 
                     "Awaiting verification", 
                     "Population on 30/06/2020") 
  
  store_dates <- store_dates[store_dates %nin% dates_exclude]
  fixed_dates <- as.Date(as.numeric(store_dates), origin = "1899-12-30")
  
  data_deaths <- data_deaths[,-c(2:9)]
  data_deaths <- data_deaths[-c(1:4,6),]
  
  colnames(data_deaths)    <- data_deaths[1,]
  colnames(data_deaths)[1] <- "Age"
  data_deaths <- data_deaths[-1,]
  
  lookup_table_deaths <- data.frame(ID      = unique(data_deaths$Age),
                                    Mapping = age_mapping_deaths)
  
  map_dates    <- colnames(data_deaths)[-1]
  map_dates_df <- data.frame(ID      = map_dates,
                             Mapping = fixed_dates)
  
  data_deaths_long <- tidyr::pivot_longer(data_deaths, 
                                          -c(Age), 
                                          values_to = "Deaths", 
                                          names_to  = "Date") 
  
  data_deaths_long <- data_deaths_long %>% 
                      as.data.frame() %>% 
                      dplyr::left_join(map_dates_df, 
                                       by = c(Date = "ID")) %>% 
                      dplyr::select(-dplyr::one_of(c("Date"))) %>% 
                      dplyr::rename(Date = Mapping,
                                    Group = Age)
  
  #---- Calculate the new deaths and cases per age group:
  data_groups_cols             <- c("Date", 
                                    "Group", 
                                    "Cumul_cases", 
                                    "Cumul_deaths")
  dt_analysis_period           <- data.frame(matrix(ncol = length(data_groups_cols), nrow = 0))
  colnames(dt_analysis_period) <- data_groups_cols
  
  #---- Rename the columns according to the mapping:
  age_bands    <- unique(unique(lookup_table_deaths$Mapping))
  age_bands_no <- length(age_bands)
  dates        <- unique(data_deaths_long$Date)
  dates        <- dates[order(dates , decreasing = FALSE )]
  
  data_deaths_long$Deaths <- as.numeric(data_deaths_long$Deaths)
  
  data_deaths_melt <- data_deaths_long %>% 
                      dplyr::left_join(lookup_table_deaths, 
                                       by = c(Group = "ID")) %>% 
                      dplyr::select(-dplyr::one_of(c("Group"))) %>% 
                      dplyr::group_by(Date, Mapping) %>% 
                      dplyr::summarise(Cumul_deaths = sum(Deaths)) %>% 
                      dplyr::rename(Group = Mapping)
  
  lookup_table_cases <- data.frame(Initial = unique(data_cases_long$age),
                                   Mapping = age_mapping_cases)
  
  data_cases_melt <- data_cases_long %>% 
                     dplyr::left_join(lookup_table_cases, 
                                      by = c(age = "Initial")) %>% 
                     dplyr::select(-dplyr::one_of(c("age"))) %>% 
                     dplyr::group_by(date, Mapping) %>% 
                     dplyr::summarise(New_Cases = sum(cases)) %>% 
                     dplyr::rename(Group = Mapping,
                                   Date  = date) %>% 
                     dplyr::filter(Date >= min(dates) & Date <= max(dates))
                     # Match with the dates of available mortality records
  
  data_cumcases_long <- data_cases_long %>% 
                        dplyr::group_by(age) %>% 
                        dplyr::arrange(date) %>% 
                        dplyr::mutate(Cumul_cases = cumsum(cases)) %>% 
                        dplyr::select(-dplyr::one_of(c("cases"))) %>% 
                        dplyr::left_join(lookup_table_cases, 
                                         by = c(age = "Initial"))
  
  data_cumcases_long$age <- NULL
  
  data_cumcases_long <- data_cumcases_long %>% 
                        dplyr::group_by(date, Mapping) %>% 
                        dplyr::summarise(Cumul_cases = sum(Cumul_cases)) %>% 
                        dplyr::rename(Group = Mapping,
                                      Date  = date) %>% 
                        dplyr::filter(Date >= min(dates) & Date <= max(dates))
  
  for (i in 1:age_bands_no){
    
    data_deaths_subset <- data_deaths_melt[data_deaths_melt$Group == age_bands[i],]
    data_cases_subset  <- data_cases_melt[data_cases_melt$Group   == age_bands[i],]
    
    dt_age_grp <- data.frame(Date         = dates,
                             Group        = rep(age_bands[i], length(dates)),
                             Cumul_deaths = c(data_deaths_subset[data_deaths_subset$Group == age_bands[i], "Cumul_deaths"])
    ) 
    
    # Calculate new deaths and cases from the cumulative deaths and cases, respectively:
    new_deaths_age_grp    <- rep(0, nrow(dt_age_grp))
    new_deaths_age_grp[1] <- data_deaths_subset$Cumul_deaths[1]
    
    for (t in 2:length(dates)) {
      
      new_deaths_age_grp[t] <- data_deaths_subset$Cumul_deaths[t] - data_deaths_subset$Cumul_deaths[t-1]
      
      # Sanity check - if the cumulative cases/deaths are lower than the day before, set to zero.
      # Confirmed by EODY reports for the dates 2020-04-09 (65+), 2020-05-05 (40-64), 2020-05-21 (40-64), 2020-05-27 (65+)
      if(new_deaths_age_grp[t] < 0) new_deaths_age_grp[t] <- 0
      
    }# End for
    
    dt_age_grp$New_Cases  <- data_cases_subset$New_Cases
    dt_age_grp$New_Deaths <- new_deaths_age_grp
    
    dt_analysis_period <- rbind(dt_analysis_period, dt_age_grp)  
    
  }# End for  
  
  #---- Subset for the analysis period:
  dt_analysis_period <- dt_analysis_period[ dt_analysis_period$Date >= as.Date(start_date, format = "%Y-%m-%d") & 
                                              dt_analysis_period$Date <= as.Date(end_date, format = "%Y-%m-%d"),]
  
  #---- Mortality data to be processed by the .stan file:
  dt_mortality_analysis_period <- reshape(dt_analysis_period[,c(1,2,ncol(dt_analysis_period))],
                                          idvar     = "Date",
                                          timevar   = "Group",
                                          direction = "wide")
  colnames(dt_mortality_analysis_period) <- gsub("New_Deaths.", "",
                                                 colnames(dt_mortality_analysis_period))
  
  dt_mortality_analysis_period <- dt_mortality_analysis_period %>% 
                                  dplyr::select(-dplyr::one_of(c("Date"))) %>% 
                                  dplyr::mutate(New_Deaths = rowSums(.))
  
  #---- Incidence data to be processed by the .stan file:
  dt_cases_analysis_period <- reshape(dt_analysis_period[,c(1,2,ncol(dt_analysis_period)-1)],
                                      idvar     = "Date",
                                      timevar   = "Group",
                                      direction = "wide")
  colnames(dt_cases_analysis_period) <- gsub("New_Cases.", "",
                                             colnames(dt_cases_analysis_period))
  
  dt_cases_analysis_period <- dt_cases_analysis_period %>% 
                              dplyr::select(-dplyr::one_of(c("Date"))) %>% 
                              dplyr::mutate(New_Cases = rowSums(.))
  
  # Additional fields:
  date_subset <- unique(dt_analysis_period$Date)
  
  dt_mortality_analysis_period$Date  <- date_subset
  dt_cases_analysis_period$Date      <- dt_mortality_analysis_period$Date
  
  dt_mortality_analysis_period$Index <- 1:nrow(dt_mortality_analysis_period)
  dt_cases_analysis_period$Index     <- 1:nrow(dt_cases_analysis_period)
  
  #---- Add Week ID since first day of analysis:
  dt_mortality_analysis_period$Week_ID <- as.numeric(dt_mortality_analysis_period$Date - dt_mortality_analysis_period$Date[1]) %/% 7 + 1
  dt_cases_analysis_period$Week_ID     <- as.numeric(dt_cases_analysis_period$Date - dt_cases_analysis_period$Date[1]) %/% 7 + 1
  
  #---- Left and right time limit:
  dt_mortality_analysis_period$Right <- dt_mortality_analysis_period$Index + 1
  dt_cases_analysis_period$Right     <- dt_cases_analysis_period$Index + 1
  
  #---- Re-order the column names"
  fixed_cols        <- c("Index", "Right", "Date", "Week_ID")
  fixed_cols_cases  <- c(fixed_cols, "New_Cases")
  fixed_cols_deaths <- c(fixed_cols, "New_Deaths")
  
  reorder_colnames <- c(fixed_cols_cases, colnames(dt_cases_analysis_period)[!colnames(dt_cases_analysis_period) %in% fixed_cols_cases])
  dt_cases_analysis_period <- dt_cases_analysis_period[reorder_colnames]
  
  reorder_colnames <- c(fixed_cols_deaths, colnames(dt_mortality_analysis_period)[!colnames(dt_mortality_analysis_period) %in% fixed_cols_deaths])
  dt_mortality_analysis_period <- dt_mortality_analysis_period[reorder_colnames]
  
  #--- Table at level (Date-Group) with the cumulative cases (overall and dissaggregated) for the analysis period:
  data_cumcases_wide <- as.data.frame(data_cumcases_long)
  data_cumcases_wide <- reshape(data_cumcases_wide,
                                idvar     = "Date",
                                timevar   = "Group",
                                direction = "wide")
  colnames(data_cumcases_wide) <- gsub("Cumul_cases.", "",
                                       colnames(data_cumcases_wide))

  data_cumcases_wide <- data_cumcases_wide[ data_cumcases_wide$Date >= as.Date(start_date, format = "%Y-%m-%d") & 
                                            data_cumcases_wide$Date <= as.Date(end_date,   format = "%Y-%m-%d"),]
  
  data_cumcases_wide <- data_cumcases_wide %>% 
                        dplyr::select(-dplyr::one_of(c("Date"))) %>% 
                        dplyr::mutate(Total_Cases = rowSums(.))
  
  data_cumcases_wide$Date  <- date_subset
  
  fixed_cols_cumcases <- c("Date", "Total_Cases")
  reorder_colnames_cumcases <- c(fixed_cols_cumcases, colnames(data_cumcases_wide)[!colnames(data_cumcases_wide) %in% fixed_cols_cumcases])
  
  data_cumcases_wide <- data_cumcases_wide[reorder_colnames_cumcases]
  
  output_list <- list(Deaths                = dt_mortality_analysis_period,
                      Cases                 = dt_cases_analysis_period,
                      All                   = dt_analysis_period,
                      Age_distribution      = age_distribution,
                      Age_distribution_aggr = age_distribution_aggr,
                      Cum_Cases             = data_cumcases_wide)
  
  return(output_list)
  
}# End function

data_engineering_observations <- function(country            = c("Greece", "Austria", "England"),
                                          start_date         = NULL,
                                          end_date           = NULL,
                                          data_cases_AT      = NULL, # Set to "Austria_5yo_bands_data",
                                          age_mapping_deaths = NULL,
                                          age_mapping_cases  = NULL)
{
  
  if (country == "Greece"){
    
    data_all <- data_engineering_observations_GR(start_date         = start_date,
                                                 end_date           = end_date,
                                                 age_mapping_deaths = age_mapping_deaths)
    
  } else if (country == "Austria"){
   
    data_all <- data_engineering_observations_AT(start_date         = start_date,
                                                 end_date           = end_date,
                                                 data_cases         = data_cases_AT,
                                                 age_mapping_deaths = age_mapping_deaths,
                                                 age_mapping_cases  = age_mapping_cases)
  
  } else if (country == "England"){
    
    data_all <- data_engineering_observations_ENG(start_date         = start_date,
                                                  end_date           = end_date,
                                                  age_mapping_deaths = age_mapping_deaths,
                                                  age_mapping_cases  = age_mapping_cases)
  
  }# End if
  
  if ( country %in% c("Greece", "Austria") ){
    
  output <- list(Deaths    = data_all$Deaths,
                 Cases     = data_all$Cases,
                 All       = data_all$All)
  
  } else if ( country == "England" ){
    
    output <- list(Deaths                = data_all$Deaths,
                   Cases                 = data_all$Cases,
                   All                   = data_all$All,
                   Age_distribution      = data_all$Age_distribution,
                   Age_distribution_aggr = data_all$Age_distribution_aggr,
                   Cusum_cases           = data_all$Cum_Cases
                   )
    
  }# End if
  
  return(output)
  
}# End function

data_engineering_npi <- function(npi_file_path,
                                 country,
                                 start_date,
                                 end_date
){
  
  npi_data_dt            <- readxl::read_xlsx(npi_file_path, col_types = c("text", "text", "date", "date"))
  npi_data_dt$date_start <- gsub(npi_data_dt$date_start, pattern = " 00:00:00", replacement = "", fixed = T)
  npi_data_dt$date_start <- as.Date(npi_data_dt$date_start)
  
  npi_data_dt$date_end <- gsub(npi_data_dt$date_end, pattern = " 00:00:00", replacement = "", fixed = T)
  npi_data_dt$date_end <- as.Date(npi_data_dt$date_end)
  
  filter1 <- npi_data_dt$Country == country
  
  npi_data_dt <- npi_data_dt[filter1,]
  
  filter2 <- npi_data_dt$date_start >= start_date & npi_data_dt$date_start <= end_date
  
  npi_data_dt <- npi_data_dt[filter2,]
  
  npi_data_dt <- npi_data_dt[order(npi_data_dt$date_start),]
  
  return(npi_data_dt)
  
}# End function

###########################################################################
#
# Estimated aggregate mortality and infection counts:
#
###########################################################################
posterior_aggregated_counts <- function(mortality_data,
                                        infections_data,
                                        nuts_fit){
  
  #---- Checks:
  if(class(nuts_fit)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(nuts_fit)
  
  #---- Deaths (posterior median, 50% and 95% credible intervals):
  fit_dead        <- posterior_draws$E_deaths
  median_fit_dead <- apply(fit_dead, 2, median)
  low_fit_dead    <- apply(fit_dead, 2, quantile, probs = c(0.025)) 
  high_fit_dead   <- apply(fit_dead, 2, quantile, probs = c(0.975))
  low25_fit_dead  <- apply(fit_dead, 2, quantile, probs = c(0.25))
  high75_fit_dead <- apply(fit_dead, 2, quantile, probs = c(0.75))
  
  deaths_output        <- mortality_data 
  deaths_output$median <- median_fit_dead
  deaths_output$low    <- low_fit_dead
  deaths_output$high   <- high_fit_dead
  deaths_output$low25  <- low25_fit_dead
  deaths_output$high75 <- high75_fit_dead
  
  #---- Infections (posterior median, 50% and 95% credible intervals):
  fit_cases        <- posterior_draws$E_cases
  median_fit_cases <- apply(fit_cases, 2, median)
  low_fit_cases    <- apply(fit_cases, 2, quantile, probs = c(0.025)) 
  high_fit_cases   <- apply(fit_cases, 2, quantile, probs = c(0.975)) 
  low25_fit_cases  <- apply(fit_cases, 2, quantile, probs = c(0.25))
  high75_fit_cases <- apply(fit_cases, 2, quantile, probs = c(0.75))
  
  infections_output        <- infections_data
  infections_output$median <- median_fit_cases
  infections_output$low    <- low_fit_cases
  infections_output$high   <- high_fit_cases
  infections_output$low25  <- low25_fit_cases
  infections_output$high75 <- high75_fit_cases
  
  output <- list(Deaths     = deaths_output,
                 Infections = infections_output)
  
  return(output)
  
}# End function

###########################################################################
#
# Estimated age-specific mortality and infection counts:
#
###########################################################################
posterior_age_specific_counts <- function(mortality_data, 
                                          infections_data,
                                          nuts_fit,
                                          cov_data,
                                          aggr_age){
  
  #---- Checks:
  if(class(nuts_fit)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(nuts_fit)
  
  #---- Placeholder for output tables:
  data_cases_cols      <- c("Date", "Week_ID", "Group", "New_Cases", "median", "low", "high")
  data_cases           <- data.frame(matrix(ncol = length(data_cases_cols), nrow = 0))
  colnames(data_cases) <- data_cases_cols
  
  data_deaths_cols      <- c("Date", "Week_ID", "Group", "New_Deaths", "median", "low", "high")
  data_deaths           <- data.frame(matrix(ncol = length(data_deaths_cols), nrow = 0))
  colnames(data_deaths) <- data_deaths_cols
  
  # Loop over each age group:
  for (i in 1:cov_data$A){
    
    fit_cases  <- posterior_draws$E_casesByAge[,,i] 
    fit_deaths <- posterior_draws$E_deathsByAge[,,i]
    
    dt_cases_age_grp  <- data.frame(Date  = infections_data$Date,
                                    Group = rep(aggr_age$AgeGrp[i], nrow(infections_data)))
    dt_deaths_age_grp <- dt_cases_age_grp
    
    # Bring the Week_ID field:
    dt_cases_age_grp  <- dt_cases_age_grp %>% 
                         left_join(infections_data[,c("Date", "Week_ID")], by = "Date")
    dt_deaths_age_grp <- dt_deaths_age_grp %>% 
                         left_join(mortality_data[,c("Date", "Week_ID")], by = "Date")
    
    dt_cases_analysis_period_melt <- infections_data %>% 
                                     dplyr::select(-one_of("Index","Week_ID", "Right")) %>% 
                                     tidyr::gather(Group, New_cases, -Date)
    
    dt_mortality_analysis_period_melt <- mortality_data %>% 
                                         dplyr::select(-one_of("Index","Week_ID", "Right")) %>%
                                         tidyr::gather(Group, New_deaths, -Date)
    
    dt_cases_age_grp <- dt_cases_age_grp %>%
                        left_join(dt_cases_analysis_period_melt,
                                  by = c("Date"  = "Date",
                                         "Group" = "Group"))
    
    dt_deaths_age_grp <- dt_deaths_age_grp %>%
                         left_join(dt_mortality_analysis_period_melt,
                                   by = c("Date"  = "Date",
                                          "Group" = "Group"))
    
    # Add quantiles from the model outputs:
    dt_cases_age_grp$median <- apply(fit_cases, 2, median)
    dt_cases_age_grp$low    <- apply(fit_cases, 2, quantile, probs = c(0.025))
    dt_cases_age_grp$high   <- apply(fit_cases, 2, quantile, probs = c(0.975)) 
    dt_cases_age_grp$low25  <- apply(fit_cases, 2, quantile, probs = c(0.25)) 
    dt_cases_age_grp$high75 <- apply(fit_cases, 2, quantile, probs = c(0.75))
    
    dt_deaths_age_grp$median <- apply(fit_deaths, 2, median)
    dt_deaths_age_grp$low    <- apply(fit_deaths, 2, quantile, probs = c(0.025))
    dt_deaths_age_grp$high   <- apply(fit_deaths, 2, quantile, probs = c(0.975))
    dt_deaths_age_grp$low25  <- apply(fit_deaths, 2, quantile, probs = c(0.25))
    dt_deaths_age_grp$high75 <- apply(fit_deaths, 2, quantile, probs = c(0.75))
    
    data_cases  <- rbind(data_cases, dt_cases_age_grp)
    data_deaths <- rbind(data_deaths, dt_deaths_age_grp)
    
  }# End for
  
  output <- list(Deaths     = data_deaths,
                 Infections = data_cases)
  
  return(output)
  
}# End function

###########################################################################
#
# Posterior random draws of age-specific mortality counts:
#
###########################################################################
deaths_random_draws <- function(dispersion_type = "time-independent",
                                cov_data,
                                model_out,
                                age_specific_deaths,
                                aggregated_deaths){
  
  set.seed(1)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  E_deathsByAge   <- posterior_draws[["E_deathsByAge"]]
  E_deaths        <- posterior_draws[["E_deaths"]]
  phi             <- posterior_draws[["phiD"]]
  
  mcmc_length     <- dim(E_deathsByAge)[1]
  ts_length       <- dim(E_deathsByAge)[2]
  age_grps        <- dim(E_deathsByAge)[3]
  dates           <- age_specific_deaths$Date
  
  death_draws           <- array(NA, c(mcmc_length, ts_length, age_grps))
  data_deaths_cols      <- c("Date",
                             "Group",
                             "rng_low",
                             "rng_low25",
                             "rng_high75",
                             "rng_high")
  age_rng_draws           <- data.frame(matrix(ncol = length(data_deaths_cols),
                                               nrow = 0))
  colnames(age_rng_draws) <- data_deaths_cols
  aggregated_rng_draws    <- age_rng_draws
  
  message(" > Estimation for age-spefic deaths")
  
  for (k in 1:age_grps) {
    message(paste0(" > Estimation in group ", k))
    
    for (i in 1:mcmc_length) {
      for (j in 1:ts_length) {
        
        if(dispersion_type == "time-independent"){
          death_draws[i,j,k] <- rnbinom(1,
                                        mu   = E_deathsByAge[i,j,k],
                                        size = phi[i])
          
        } else {
          death_draws[i,j,k] <- rnbinom(1,
                                        mu   = E_deathsByAge[i,j,k],
                                        size = (E_deathsByAge[i,j,k] / phi[i]))
          
        }# End if
      }# End for
    }# End for
  }# End for
  
  for (k in 1:age_grps) {
    rng_draws_age_grp  <- data.frame(Date  = dates,
                                     Group = rep( colnames(cov_data$y_deaths)[k],
                                                  length(dates)
                                     ))
    rng_draws_age_grp$rng_low    <- apply(death_draws[,,k], 2, quantile, probs = c(0.025))
    rng_draws_age_grp$rng_low25  <- apply(death_draws[,,k], 2, quantile, probs = c(0.25))
    rng_draws_age_grp$rng_high75 <- apply(death_draws[,,k], 2, quantile, probs = c(0.75))
    rng_draws_age_grp$rng_high   <- apply(death_draws[,,k], 2, quantile, probs = c(0.975))
    
    age_rng_draws <- rbind(age_rng_draws, rng_draws_age_grp)
  }# End for
  
  message(" > Estimation for aggregated deaths")
  
  aggregated_death_draws <- apply(death_draws, MARGIN = c(1,2), sum)
  aggregated_rng_draws   <- data.frame(Date  = dates)
  
  aggregated_rng_draws$rng_low    <- apply(aggregated_death_draws, 2, quantile, probs = c(0.025)) 
  aggregated_rng_draws$rng_low25  <- apply(aggregated_death_draws, 2, quantile, probs = c(0.25))  
  aggregated_rng_draws$rng_high75 <- apply(aggregated_death_draws, 2, quantile, probs = c(0.75))  
  aggregated_rng_draws$rng_high   <- apply(aggregated_death_draws, 2, quantile, probs = c(0.975))
  
  age_rng_draws        <- age_rng_draws %>% distinct(Date, Group, .keep_all = TRUE)
  aggregated_rng_draws <- aggregated_rng_draws %>% distinct(Date, .keep_all = TRUE)
  
  output <- list(age_specific_deaths = age_rng_draws,
                 aggregated_deaths   = aggregated_rng_draws)
  
  return(output)
  
}# End function
deaths_random_draws <- compiler::cmpfun(deaths_random_draws)

###########################################################################
#
# Age-specific daily reporting ratio:
#
###########################################################################
posterior_reporting_ratio <- function(infections_age_data,
                                      infections_plot_data,
                                      nuts_fit,
                                      cov_data,
                                      age_mapping_deaths,
                                      reporting_delay_days = 6){
  
  #---- Checks:
  if(class(nuts_fit)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(nuts_fit)
  
  #---- Placeholder for output tables:
  output_cols     <- c("Date", "Group", "median", "low", "high", "low_25", "high_75")
  output          <- data.frame(matrix(ncol = length(output_cols), nrow = 0))
  colnames(output)<- output_cols
  
  group_labels <- unique(age_mapping_deaths)
  
  for (i in 1:cov_data$A){
    
    fit_cases  <- posterior_draws$E_casesByAge[,,i]
    
    # Initiate the data frame with zeros:
    reporting_ratio_group           <- data.frame( matrix(0, nrow = nrow(fit_cases), ncol = ncol(fit_cases) ) )
    colnames(reporting_ratio_group) <- paste0("Week_", 1:ncol(fit_cases))
    
    data_cases_filt <- subset(infections_age_data, Group == group_labels[i])
    
    for (t in (reporting_delay_days + 1):ncol(fit_cases) ) reporting_ratio_group[,t] <- data_cases_filt$New_cases[t] / fit_cases[,(t - reporting_delay_days)]
    
    reporting_ratio_group <- reporting_ratio_group[,(reporting_delay_days + 1):ncol(fit_cases)]
    
    output_group <- data.frame(Date    = infections_plot_data$Date[(reporting_delay_days + 1):length(infections_plot_data$Date)],
                               Group   = rep(group_labels[i], ncol(reporting_ratio_group)),
                               median  = apply(reporting_ratio_group, 2, median),
                               low     = apply(reporting_ratio_group, 2, quantile, probs = c(0.025)),
                               high    = apply(reporting_ratio_group, 2, quantile, probs = c(0.975)),
                               low_25  = apply(reporting_ratio_group, 2, quantile, probs = c(0.25)),
                               high_75 = apply(reporting_ratio_group, 2, quantile, probs = c(0.75)),
                               stringsAsFactors = FALSE)
    
    output <- rbind(output, output_group)
    rm(output_group)
  }# End for
  
  #---- Set "Date" column to class "Date":
  output$Date <- as.Date(output$Date)
  
  return(output)
  
}# End function  

###########################################################################
#
# Transmissibility, age-specific transmission rate, effective 
# reproduction number:
#
###########################################################################
posterior_infection_dynamics <- function(cov_data,
                                         model_out,
                                         data_cases,
                                         data_rt){

  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  
  #---- Checks:
  if(ncol(posterior_draws$cm_sample) != cov_data$A) stop( paste0("The number of rows in the age distribution table must be equal to ", cov_data$A) )
  
  age_grps <- cov_data$A
  
  #---- Identify between weekly and daily GBM:
  if ( "beta_weekly" %in% names(posterior_draws) ){ 

    beta_draws   <- posterior_draws$beta_weekly
    chain_length <- nrow(beta_draws)
    ts_length    <- dim(beta_draws)[2]
    dates        <- data_cases %>%
                    dplyr::select(Date, Week_ID) %>%
                    dplyr::group_by((Week_ID)) %>%
                    dplyr::summarise(Date = min(Date)) %>% 
                    dplyr::pull(Date)
    
  } else if ( "beta_N" %in% names(posterior_draws) ) {
    
    beta_draws   <- posterior_draws$beta_N
    chain_length <- nrow(beta_draws)
    
    if( !is.na(dim(beta_draws)[3]) ){
    #-- Case of multi-BM

      # Keep every 7 days if we assume weekly GBMs:
      if("n_weeks" %in% names(cov_data) ){
        beta_draws2 <- array(NA, c(chain_length, cov_data$n_weeks, cov_data$A))
        
        for (k in 1:cov_data$A) beta_draws2[,,k] <- t( t(beta_draws[,,k]) %>% 
                                                       as.data.frame %>% 
                                                       slice(which(row_number() %% 7 == 1)) )
        
        ts_length    <- dim(beta_draws2)[2]
        dates        <- data_cases %>%
                        dplyr::select(Date, Week_ID) %>%
                        dplyr::group_by((Week_ID)) %>%
                        dplyr::summarise(Date = min(Date)) %>% 
                        dplyr::pull(Date)
        
        beta_draws <- beta_draws2
        rm(beta_draws2)
        
      } else {
        ts_length    <- dim(beta_draws)[2]
        dates        <- data_cases$Date
      }# End if 
      
    } else {
    #-- Case of single-BM
      
      # Keep every 7 days if we assume weekly GBMs:
      if("n_weeks" %in% names(cov_data) ){
        
        beta_draws   <- t( t(beta_draws) %>% as.data.frame %>% slice(which(row_number() %% 7 == 1)))
        
        ts_length    <- dim(beta_draws)[2]
        dates        <- data_cases %>%
                        dplyr::select(Date, Week_ID) %>%
                        dplyr::group_by((Week_ID)) %>%
                        dplyr::summarise(Date = min(Date)) %>% 
                        dplyr::pull(Date)
        
      } else {
        ts_length    <- dim(beta_draws)[2]
        dates        <- data_cases$Date
      }# End if 
      
    }# End if
    
  }# End if

  #---- Initiate output table for the transmission rate
  #     per age group over time:
  data_transmission_rate_cols      <- c("Date", 
                                        "Group", 
                                        "median", 
                                        "low0025", 
                                        "low25", 
                                        "high75", 
                                        "high975")
  data_transmission_rate           <- data.frame(matrix(ncol = length(data_transmission_rate_cols), nrow = 0))
  colnames(data_transmission_rate) <- data_transmission_rate_cols
  
  message(" > Estimation of age-specific transmission rate")
  
  for (k in 1:age_grps){
    
    trans_rate_temp          <- matrix(0L, nrow = chain_length, ncol = ts_length)
    data_trans_rate_age_grp  <- data.frame(Date  = dates,
                                           Group = rep( colnames(cov_data$y_deaths)[k],
                                                        length(dates)
                                                        ))
    
    for (j in 1:ts_length){
      # Common GBM across groups:
      if( is.na(dim(beta_draws)[3]) ){
        
        trans_rate_temp[,j] <- beta_draws[,j] * posterior_draws$cm_sample[,k,k]
        
        # Multiple GBMs: 
      } else if ( !is.na(dim(beta_draws)[3]) ){
        
        trans_rate_temp[,j] <- beta_draws[,j,k] * posterior_draws$cm_sample[,k,k]
        
      }# End if
    }# End for 
      
    data_trans_rate_age_grp$median  <- apply(trans_rate_temp, 2, median)
    data_trans_rate_age_grp$low0025 <- apply(trans_rate_temp, 2, quantile, probs = c(0.025)) # c(0.025)
    data_trans_rate_age_grp$low25   <- apply(trans_rate_temp, 2, quantile, probs = c(0.25))  # c(0.025)
    data_trans_rate_age_grp$high75  <- apply(trans_rate_temp, 2, quantile, probs = c(0.75))  # c(0.975)
    data_trans_rate_age_grp$high975 <- apply(trans_rate_temp, 2, quantile, probs = c(0.975)) # c(0.975)
    
    data_transmission_rate <- rbind(data_transmission_rate,
                                    data_trans_rate_age_grp)
    
  }# End for
  
  message(" > Estimation of effective contact rate")
  
  #---- Effective contact rate
  # Identify between common and multiple GBM:
  if( is.na(dim(beta_draws)[3]) ){
   
    data_eff_contact_rate  <- data.frame(Date = dates)
    
    # Add quantiles from the model outputs:
    data_eff_contact_rate$Beta_median  <- apply(beta_draws, 2, median)
    data_eff_contact_rate$Beta_low0025 <- apply(beta_draws, 2, quantile, probs = c(0.025)) # c(0.025)
    data_eff_contact_rate$Beta_low25   <- apply(beta_draws, 2, quantile, probs = c(0.25))  # c(0.025)
    data_eff_contact_rate$Beta_high75  <- apply(beta_draws, 2, quantile, probs = c(0.75))  # c(0.975)
    data_eff_contact_rate$Beta_high975 <- apply(beta_draws, 2, quantile, probs = c(0.975)) # c(0.975)
    
  } else if ( !is.na(dim(beta_draws)[3]) ){

    data_eff_contact_rate_cols     <- c("Date", 
                                        "Group", 
                                        paste("Beta", c("median", "low0025", "low25", "high75", "high975")))
    data_eff_contact_rate           <- data.frame(matrix(ncol = length(data_eff_contact_rate_cols), nrow = 0))
    colnames(data_eff_contact_rate) <- data_eff_contact_rate_cols
    
    for (i in 1:cov_data$A){
      fit_beta <- beta_draws[,,i]
      
      data_eff_contact_rate_age_grp  <- data.frame(Date  = dates,
                                                   Group = rep( colnames(cov_data$y_deaths)[i], 
                                                                length(dates)))
      
      # Add quantiles from the model outputs:
      data_eff_contact_rate_age_grp$Beta_median  <- apply(fit_beta, 2, median)
      data_eff_contact_rate_age_grp$Beta_low0025 <- apply(fit_beta, 2, quantile, probs = c(0.025)) # c(0.025)
      data_eff_contact_rate_age_grp$Beta_low25   <- apply(fit_beta, 2, quantile, probs = c(0.25)) # c(0.025)
      data_eff_contact_rate_age_grp$Beta_high75  <- apply(fit_beta, 2, quantile, probs = c(0.75)) # c(0.975)
      data_eff_contact_rate_age_grp$Beta_high975 <- apply(fit_beta, 2, quantile, probs = c(0.975)) # c(0.975)
      
      data_eff_contact_rate <- rbind(data_eff_contact_rate, 
                                     data_eff_contact_rate_age_grp)
    
    }# End for
  }# End if
    
  message(" > Estimation of reproduction number")
  
  data_repnumber                <- data.frame(Date  = dates)
  data_repnumber$eff_rt_median  <- apply(post_Rt$R_t, 2, median)
  data_repnumber$eff_rt_low0025 <- apply(post_Rt$R_t, 2, quantile, probs = c(0.025)) # c(0.025)
  data_repnumber$eff_rt_low25   <- apply(post_Rt$R_t, 2, quantile, probs = c(0.25))  # c(0.025)
  data_repnumber$eff_rt_high75  <- apply(post_Rt$R_t, 2, quantile, probs = c(0.75))  # c(0.975)
  data_repnumber$eff_rt_high975 <- apply(post_Rt$R_t, 2, quantile, probs = c(0.975)) # c(0.975)
  
  message(" > Estimation of effective reproduction number")

  data_eff_repnumber                <- data.frame(Date  = dates)
  data_eff_repnumber$eff_rt_median  <- apply(post_Rt$R_eff_t, 2, median)
  data_eff_repnumber$eff_rt_low0025 <- apply(post_Rt$R_eff_t, 2, quantile, probs = c(0.025)) # c(0.025)
  data_eff_repnumber$eff_rt_low25   <- apply(post_Rt$R_eff_t, 2, quantile, probs = c(0.25))  # c(0.025)
  data_eff_repnumber$eff_rt_high75  <- apply(post_Rt$R_eff_t, 2, quantile, probs = c(0.75))  # c(0.975)
  data_eff_repnumber$eff_rt_high975 <- apply(post_Rt$R_eff_t, 2, quantile, probs = c(0.975)) # c(0.975)
  
  #---- Output:
  output <- list(R_t               = data_repnumber,
                 R_eff_t           = data_eff_repnumber,
                 GBM               = data_eff_contact_rate,
                 Transmission_rate = data_transmission_rate)
  
  return(output)
  
}# End function

###########################################################################
#
# Age distribution - adjusted figures of estimated infections 
# (Table 2, Ward 2021):
#
###########################################################################
aggregate_react_cumcases <- function(x, 
                                     user_AgeGrp) 
{
  options(dplyr.summarise.inform = FALSE)
  if (length(user_AgeGrp) != nrow(x)) 
    stop("The mapped age group labels do not correspond to the age group labels of the aggregated age distribution matrix.\n")
  
  REACT2_cumul_cases <- data.frame(AgeGrp = c("0-14", "15-44", "45-64", "65-74", "74-100"),
                                   Mean   = c(0, 1536 * 1e+03, 895 * 1e+03 ,181 * 1e+03, 166 * 1e+03),
                                   Lower  = c(0, 1437 * 1e+03, 837 * 1e+03 ,153 * 1e+03, 131 * 1e+03),
                                   Upper  = c(0, 1635 * 1e+03, 953 * 1e+03 ,209 * 1e+03, 201 * 1e+03) )
  
  REACT2_cumul_cases$AgeGrpStart <- sapply(1:nrow(REACT2_cumul_cases), function(x) {
    min(as.numeric(strsplit(REACT2_cumul_cases$AgeGrp, "-")[[x]]))
  })
  
  REACT2_cumul_cases$AgeGrpEnd <- sapply(1:nrow(REACT2_cumul_cases), function(x) {
    max(as.numeric(strsplit(REACT2_cumul_cases$AgeGrp, "-")[[x]]))
  })
  
  temp_x <- x
  temp_x$AgeGrp <- gsub("\\+", "-100", temp_x$AgeGrp)
  temp_x$AgeGrpEnd <- sapply(1:nrow(temp_x), function(x) {
    max(as.numeric(strsplit(temp_x$AgeGrp, "-")[[x]]))
  })
  
  temp_x$Group_mapping <- user_AgeGrp
  
  if ("AgeGrpStart" %nin% colnames(temp_x)) 
    temp_x$AgeGrpStart <- sapply(1:nrow(temp_x), function(x) {
      min(as.numeric(strsplit(temp_x$AgeGrp, "-")[[x]]))
    })
  
  for (i in 1:nrow(temp_x)) {
    for (j in 1:nrow(REACT2_cumul_cases)) {
      if ( (temp_x$AgeGrpStart[i] >= REACT2_cumul_cases$AgeGrpStart[j]) & 
           (temp_x$AgeGrpEnd[i]   <= REACT2_cumul_cases$AgeGrpEnd[j]) ) 
        
        temp_x$Mean[i]  <- REACT2_cumul_cases$Mean[j]
      
    }
  }
  
  for (i in 1:nrow(temp_x)) {
    for (j in 1:nrow(REACT2_cumul_cases)) {
      if ( (temp_x$AgeGrpStart[i] >= REACT2_cumul_cases$AgeGrpStart[j]) & 
           (temp_x$AgeGrpEnd[i]   <= REACT2_cumul_cases$AgeGrpEnd[j]) ) 
        
        temp_x$Lower[i] <- REACT2_cumul_cases$Lower[j]
      
    }
  }
  
  for (i in 1:nrow(temp_x)) {
    for (j in 1:nrow(REACT2_cumul_cases)) {
      if ( (temp_x$AgeGrpStart[i] >= REACT2_cumul_cases$AgeGrpStart[j]) & 
           (temp_x$AgeGrpEnd[i]   <= REACT2_cumul_cases$AgeGrpEnd[j]) ) 
        
        temp_x$Upper[i] <- REACT2_cumul_cases$Upper[j]
      
    }
  }
  
  output <- temp_x %>% 
    as.data.frame() %>%
    dplyr::group_by(Group_mapping) %>% 
    dplyr::mutate(PopPerc = prop.table(PopTotal), 
                  Mean    = sum(Mean  * PopPerc),
                  Lower   = sum(Lower * PopPerc),
                  Upper   = sum(Upper * PopPerc)
    ) %>% 
    dplyr::select(dplyr::one_of(c("Group_mapping", 
                                  "Mean", "Lower", "Upper"))) %>% 
    dplyr::group_by(Group_mapping) %>% 
    dplyr::slice(1) %>% 
    as.data.frame()
  
  return(output)
}
