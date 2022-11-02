###########################################################################
#
# Paths
#
###########################################################################
path_to_source      <- "C:/Users/bouranis/OneDrive - aueb.gr/BERNADETTE/10_Stan_Project_code/11_SIR_GR_Age_model4/Github_files"
path_to_stan_output <- store_model_out #"..."

###########################################################################
#
# Load model output (stored in the object nuts_fit_1):
#
###########################################################################
load(file = path_to_stan_output)

###########################################################################
#
# Data preparation:
#
###########################################################################
if (country == "Austria"){
  data_all <- data_engineering_observations(country = country,
                                            start_date,
                                            end_date,
                                            data_cases_AT = Output_5_clean, 
                                            age_mapping_deaths,
                                            age_mapping_cases)

} else if (country == "Greece"){
  data_all <- data_engineering_observations(country = country,
                                            start_date,
                                            end_date,
                                            data_cases_AT = NULL, 
                                            age_mapping_deaths,
                                            age_mapping_deaths)

  
} else if (country == "England"){
  
  data_all <- data_engineering_observations(country = country,
                                            start_date,
                                            end_date,
                                            data_cases_AT = NULL, 
                                            age_mapping_deaths,
                                            age_mapping_cases)
  
}# End if

dt_analysis_period <- data_all$All
mortality_data     <- data_all$Deaths
infections_data    <- data_all$Cases
cusum_infections   <- data_all$Cusum_cases

age_specific_fits  <- posterior_age_specific_counts(mortality_data, 
                                                    infections_data,
                                                    nuts_fit = nuts_fit_1,
                                                    cov_data,
                                                    aggr_age)

aggregated_fits    <- posterior_aggregated_counts(mortality_data,
                                                  infections_data,
                                                  nuts_fit = nuts_fit_1)

aggregated_deaths_data   <- aggregated_fits$Deaths
aggregated_cases_data    <- aggregated_fits$Infections
age_specific_cases_data  <- age_specific_fits$Infections
age_specific_deaths_data <- age_specific_fits$Deaths

#---- Generate random samples of age-stratified deaths from the model:
deaths_rng      <- deaths_random_draws(dispersion_type     = "time-dependent",
                                       cov_data,
                                       model_out           = nuts_fit_1,
                                       age_specific_deaths = age_specific_deaths_data,
                                       aggregated_deaths   = aggregated_deaths_data)

age_deaths_rng  <- deaths_rng$age_specific_deaths
aggr_deaths_rng <- deaths_rng$aggregated_deaths

age_specific_deaths_data <- age_specific_deaths_data %>% left_join(age_deaths_rng,  by = c("Date", "Group"))
aggregated_deaths_data   <- aggregated_deaths_data   %>% left_join(aggr_deaths_rng, by = c("Date"))

#---- Estimation of the effective reproduction number:
time_start_Rt <- Sys.time()
post_Rt <- effective_rt_seeiir(cov_data     = cov_data,
                               model_out    = nuts_fit_1,
                               progress_bar = TRUE)
time_end_Rt <- Sys.time()
duration_Rt <- time_end_Rt - time_start_Rt

save(post_Rt,
     duration_Rt,
     file = rt_path)

# Load the estimates of the effective reproduction number:
# load(rt_path)

#---- Datasets related to the age-stratified transmissibility and transmission rate, 
#     and the effective reproduction number: 
data_eff       <- posterior_infection_dynamics(cov_data,
                                               model_out = nuts_fit_1,
                                               age_specific_cases_data,
                                               data_rt = post_Rt)

#---- Dataset related to the age-stratified reporting ratio:
data_rep_ratio <- posterior_reporting_ratio(infections_age_data  = age_specific_cases_data,
                                            infections_plot_data = aggregated_cases_data,
                                            nuts_fit             = nuts_fit_1,
                                            cov_data,
                                            age_mapping_deaths,
                                            reporting_delay_days = 6)

###########################################################################
#
# Data on non-pharmaceutical interventions:
#
###########################################################################
npi_data_path <- paste0(path_to_source,
                        "/ECDC_Country_NPI.xlsx")

#--- Import country-specific NPIs. Source: ECDC.
data_npi <- data_engineering_npi(npi_data_path,
                                 country,
                                 start_date,
                                 end_date)

#--- Select a few NPIs for demonstration:

# Greece:
data_npi_filter <- c(4, 18, 21, 22, 26)

# Austria:
#data_npi_filter <- c(2, 5, 6, 22, 26, 27)

data_npi        <- data_npi[data_npi_filter,]

data_npi_long <- data_npi %>% 
  tidyr::gather(date_start, date_end, -c(Country, Response_measure)) %>% 
  dplyr::mutate(NPI = ifelse(date_start == "date_start", 
                             Response_measure, 
                             paste("End of ", Response_measure, sep = ""))) %>% 
  dplyr::rename(Date_ID = date_start,
                Date    = date_end) %>% 
  dplyr::filter(Date <= end_date) %>% 
  dplyr::mutate(NPI_ID    = 1:n(),
                NPI_label = paste(NPI_ID, NPI, sep = " - "))

###########################################################################
#
# Age-stratified new observed deaths:
#
###########################################################################
country_observed_deaths <- 
  ggplot(dt_analysis_period,
       aes(Date       = Date,
           New_Deaths = New_Deaths)) +
  geom_point(aes(x   = Date,
                 y   = New_Deaths,
                 fill= "Reported deaths")) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right",
             shrink = T) +  
  scale_fill_manual(values = c('Reported deaths' = 'black'), 
                    guide  = "none") +
  labs(x = "Epidemiological Date",
       y = "Number of new deaths") +
  scale_x_date(date_breaks = "1 month", 
               date_labels =  "%b %Y") +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
		    axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank() )

###########################################################################
#
# Goodness of fit to new overall deaths:
#
###########################################################################
fit_deaths_all <- 
  ggplot(aggregated_deaths_data,
         aes(Date       = Date,
             New_Deaths = New_Deaths)) +
  geom_point(aes(x   = Date, 
                 y   = New_Deaths,
                 fill= "Reported")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = rng_low,
                  ymax = rng_high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "New daily deaths") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(
    limits = c(0,   max(aggregated_deaths_data$high)*1.4),    
    breaks = seq(0, max(aggregated_deaths_data$high)*1.4, 20) 
  ) +
  scale_fill_manual(values = c('Reported' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black", 
                                 "95% CrI" = "gray40"       
                      )
  ) +
  theme_bw() +
  theme(strip.placement  = "outside",
        strip.background = element_rect(fill = NA, colour="grey50"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        panel.spacing    = unit(0,"cm"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin()) 

###########################################################################
#
# Goodness of fit to new age-stratified deaths:
#
###########################################################################
fit_deaths_age <- 
  ggplot(age_specific_deaths_data,
         aes(Date       = Date,
             New_Deaths = New_deaths)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  geom_point(aes(x   = Date, 
                 y   = New_deaths,
                 fill= "Reported")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = rng_low,
                  ymax = rng_high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "New daily deaths") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_fill_manual(values = c('Reported' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black", 
                                 "95% CrI" = "gray40"       
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())
  
###########################################################################
#
# Aggregated infections:
#
###########################################################################

if (country == "Austria") y_NPI <- 24000 else if (country == "Greece") y_NPI <- 17000

fit_cases_all <- 
  ggplot(aggregated_cases_data,
         aes(Date       = Date,
             New_Cases = New_Cases)) +
  geom_point(aes(x   = Date, 
                 y   = New_Cases,
                 fill= "Reported")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "New daily infections") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(
    limits = c(0, max(aggregated_cases_data$high75)*1.3)  , 
    breaks = seq(0, max(aggregated_cases_data$high75)*13, 3000) 
  ) +
  scale_fill_manual(values = c('Reported' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black", 
                                 "50% CrI" = "gray40"
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank() ) +
  geom_vline(data    = data_npi_long,
             mapping = aes(xintercept = Date,
                           color   = "black"), 
             key_glyph = glyph_vline,
             linetype  = 2) +
  ggrepel::geom_text_repel(data = data_npi_long,
                           aes(x     = Date,
                               y     = y_NPI,
                               label = NPI),
                           color        = "black",
                           force_pull   = 0,
                           nudge_y      = 0.1,
                           direction    = "x",
                           angle        = 90,
                           hjust        = 0.5,
                           segment.size = 0.2,
                           max.iter     = 1e4,
                           max.time     = 1,
                           inherit.aes  = FALSE)

###########################################################################
#
# Age-stratified infections:
#
###########################################################################

#---- Austria. Plot age-stratified infections vs NPIs together and the aggregated infections independently:
if (country == "Austria"){
  age_specific_cases_data_grp1 <- subset(age_specific_cases_data, Group == "<65")
  age_specific_cases_data_grp2 <- subset(age_specific_cases_data, Group == "65+")
  
  age_specific_cases_data_grp1_plot <- 
    ggplot(age_specific_cases_data_grp1,
           aes(Date      = Date,
               New_Cases = New_cases)) +
    geom_point(aes(x   = Date, 
                   y   = New_cases,
                   fill= "Reported")) +
    geom_line(aes(x     = Date,
                  y     = median,
                  color = "Median"),
              size = 1.3) +
    geom_ribbon(aes(x    = Date,
                    ymin = low25,
                    ymax = high75,
                    fill = "50% CrI"),
                alpha = 0.5) +
    labs(x = "Epidemiological Date",
         y = "New daily infections") +
    scale_x_date(date_breaks       = "1 month",
                 date_minor_breaks = "1 month") +
    scale_y_continuous(
      limits = c(0, max(age_specific_cases_data_grp1$high75)*1.8)  , 
      breaks = seq(0, max(age_specific_cases_data_grp1$high75)*1.8, 2000) 
    ) +
    scale_fill_manual(values = c('Reported' = "black")) +
    scale_colour_manual(name = '', 
                        values = c('Median'  = "black",
                                   "50% CrI" = "gray40")) +
    theme_bw() +
    theme(panel.spacing    = unit(0.2,"cm"),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          axis.title.x     = element_text(size = 14, face = "bold"),
          axis.title.y     = element_text(size = 14, face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_blank() ) +
    geom_vline(data    = data_npi_long,
               mapping = aes(xintercept = Date,
                             color   = "black"),
               key_glyph = glyph_vline,
               linetype  = 2) +
    ggrepel::geom_text_repel(data = data_npi_long,
                             aes(x     = Date,
                                 y     = 22500,
                                 label = NPI),
                             color        = "black",
                             force_pull   = 0,
                             nudge_y      = 0.1,
                             direction    = "x",
                             angle        = 90,
                             hjust        = 0.5,
                             segment.size = 0.2,
                             max.iter     = 1e4,
                             max.time     = 1,
                             inherit.aes  = FALSE)
  
  age_specific_cases_data_grp2_plot <- 
    ggplot(age_specific_cases_data_grp2,
           aes(Date      = Date,
               New_Cases = New_cases)) +
    geom_point(aes(x   = Date, 
                   y   = New_cases,
                   fill= "Reported")) +
    geom_line(aes(x     = Date,
                  y     = median,
                  color = "Median"),
              size = 1.3) +
    geom_ribbon(aes(x    = Date,
                    ymin = low25,
                    ymax = high75,
                    fill = "50% CrI"),
                alpha = 0.5) +
    labs(x = "Epidemiological Date",
         y = "New daily infections") +
    scale_y_continuous(
      limits = c(0, max(age_specific_cases_data_grp2$high75)*1.9)  , 
      breaks = seq(0, max(age_specific_cases_data_grp2$high75)*1.9, 500) 
    ) +
    scale_x_date(date_breaks       = "1 month",
                 date_minor_breaks = "1 month") +
    scale_fill_manual(values = c('Reported' = "black")) +
    scale_colour_manual(name = '', 
                        values = c('Median'  = "black",
                                   "50% CrI" = "gray40")) +
    theme_bw() +
    theme(panel.spacing    = unit(0.2,"cm"),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          axis.title.x     = element_text(size = 14, face = "bold"),
          axis.title.y     = element_text(size = 14, face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_blank() ) +
    geom_vline(data    = data_npi_long,
               mapping = aes(xintercept = Date,
                             color   = "black"),
               key_glyph = glyph_vline,
               linetype  = 2) +
    ggrepel::geom_text_repel(data = data_npi_long,
                             aes(x     = Date,
                                 y     = 3100,
                                 label = NPI),
                             color        = "black",
                             force_pull   = 0,
                             nudge_y      = 0.1,
                             direction    = "x",
                             angle        = 90,
                             hjust        = 0.5,
                             segment.size = 0.2,
                             max.iter     = 1e4,
                             max.time     = 1,
                             inherit.aes  = FALSE)
  
  } else if (country == "Greece") {
    
    age_specific_cases_data_grp1 <- subset(age_specific_cases_data, Group == "0-39")
    age_specific_cases_data_grp2 <- subset(age_specific_cases_data, Group == "40-64")
    age_specific_cases_data_grp3 <- subset(age_specific_cases_data, Group == "65+")
    
    age_specific_cases_data_grp1_plot <- 
      ggplot(age_specific_cases_data_grp1,
             aes(Date      = Date,
                 New_Cases = New_cases)) +
      geom_point(aes(x   = Date, 
                     y   = New_cases,
                     fill= "Reported")) +
      geom_line(aes(x     = Date,
                    y     = median,
                    color = "Median"),
                size = 1.3) +
      geom_ribbon(aes(x    = Date,
                      ymin = low25,
                      ymax = high75,
                      fill = "50% CrI"),
                  alpha = 0.5) +
      labs(x = "Epidemiological Date",
           y = "New daily infections") +
      scale_x_date(date_breaks       = "1 month",
                   date_minor_breaks = "1 month") +
      scale_y_continuous(
        limits = c(0, max(age_specific_cases_data_grp1$high75)*1.3)  , 
        breaks = seq(0, max(age_specific_cases_data_grp1$high75)*1.3, 2000) 
      ) +
      scale_fill_manual(values = c('Reported' = "black")) +
      scale_colour_manual(name = '', 
                          values = c('Median'  = "black",
                                     "50% CrI" = "gray40")) +
      theme_bw() +
      theme(panel.spacing    = unit(0.2,"cm"),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            axis.title.x     = element_text(size = 14, face = "bold"),
            axis.title.y     = element_text(size = 14, face = "bold"),
            legend.position  = "bottom",
            legend.title     = element_blank() ) +
      geom_vline(data    = data_npi_long,
                 mapping = aes(xintercept = Date,
                               color   = "black"),
                 key_glyph = glyph_vline,
                 linetype  = 2) +
      ggrepel::geom_text_repel(data = data_npi_long,
                               aes(x     = Date,
                                   y     = 7000,
                                   label = NPI),
                               color        = "black",
                               force_pull   = 0,
                               nudge_y      = 0.1,
                               direction    = "x",
                               angle        = 90,
                               hjust        = 0.5,
                               segment.size = 0.2,
                               max.iter     = 1e4,
                               max.time     = 1,
                               inherit.aes  = FALSE)
    
    age_specific_cases_data_grp2_plot <- 
      ggplot(age_specific_cases_data_grp2,
             aes(Date      = Date,
                 New_Cases = New_cases)) +
      geom_point(aes(x   = Date, 
                     y   = New_cases,
                     fill= "Reported")) +
      geom_line(aes(x     = Date,
                    y     = median,
                    color = "Median"),
                size = 1.3) +
      geom_ribbon(aes(x    = Date,
                      ymin = low25,
                      ymax = high75,
                      fill = "50% CrI"),
                  alpha = 0.5) +
      labs(x = "Epidemiological Date",
           y = "New daily infections") +
      scale_y_continuous(
        limits = c(0, max(age_specific_cases_data_grp2$high75)*1.3)  , 
        breaks = seq(0, max(age_specific_cases_data_grp2$high75)*1.3, 1000) 
      ) +
      scale_x_date(date_breaks       = "1 month",
                   date_minor_breaks = "1 month") +
      scale_fill_manual(values = c('Reported' = "black")) +
      scale_colour_manual(name = '', 
                          values = c('Median'  = "black",
                                     "50% CrI" = "gray40")) +
      theme_bw() +
      theme(panel.spacing    = unit(0.2,"cm"),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            axis.title.x     = element_text(size = 14, face = "bold"),
            axis.title.y     = element_text(size = 14, face = "bold"),
            legend.position  = "bottom",
            legend.title     = element_blank() ) +
      geom_vline(data    = data_npi_long,
                 mapping = aes(xintercept = Date,
                               color   = "black"),
                 key_glyph = glyph_vline,
                 linetype  = 2) +
      ggrepel::geom_text_repel(data = data_npi_long,
                               aes(x     = Date,
                                   y     = 7000,
                                   label = NPI),
                               color        = "black",
                               force_pull   = 0,
                               nudge_y      = 0.1,
                               direction    = "x",
                               angle        = 90,
                               hjust        = 0.5,
                               segment.size = 0.2,
                               max.iter     = 1e4,
                               max.time     = 1,
                               inherit.aes  = FALSE)
    
    age_specific_cases_data_grp3_plot <- 
      ggplot(age_specific_cases_data_grp3,
             aes(Date      = Date,
                 New_Cases = New_cases)) +
      geom_point(aes(x   = Date, 
                     y   = New_cases,
                     fill= "Reported")) +
      geom_line(aes(x     = Date,
                    y     = median,
                    color = "Median"),
                size = 1.3) +
      geom_ribbon(aes(x    = Date,
                      ymin = low25,
                      ymax = high75,
                      fill = "50% CrI"),
                  alpha = 0.5) +
      labs(x = "Epidemiological Date",
           y = "New daily infections") +
      scale_y_continuous(
        limits = c(0, max(age_specific_cases_data_grp3$high75)*1.4)  , 
        breaks = seq(0, max(age_specific_cases_data_grp3$high75)*1.4, 500) 
      ) +
      scale_x_date(date_breaks       = "1 month",
                   date_minor_breaks = "1 month") +
      scale_fill_manual(values = c('Reported' = "black")) +
      scale_colour_manual(name = '', 
                          values = c('Median'  = "black",
                                     "50% CrI" = "gray40")) +
      theme_bw() +
      theme(panel.spacing    = unit(0.2,"cm"),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            axis.title.x     = element_text(size = 14, face = "bold"),
            axis.title.y     = element_text(size = 14, face = "bold"),
            legend.position  = "bottom",
            legend.title     = element_blank() ) +
      geom_vline(data    = data_npi_long,
                 mapping = aes(xintercept = Date,
                               color   = "black"),
                 key_glyph = glyph_vline,
                 linetype  = 2) +
      ggrepel::geom_text_repel(data = data_npi_long,
                               aes(x     = Date,
                                   y     = 2000,
                                   label = NPI),
                               color        = "black",
                               force_pull   = 0,
                               nudge_y      = 0.1,
                               direction    = "x",
                               angle        = 90,
                               hjust        = 0.5,
                               segment.size = 0.2,
                               max.iter     = 1e4,
                               max.time     = 1,
                               inherit.aes  = FALSE)
    
}# End if

###########################################################################
#
# Cumulative age-stratified infections (available for England only):
#
###########################################################################
unique_dates <- unique(age_specific_cases_data$Date)

#---- Locate duplicates, after rounding to 0 dps:
cusum_data_plot <- age_specific_cases_data %>% 
                   mutate_if(is.numeric, round, 0) %>% 
                   dplyr::select(-dplyr::one_of(c("Date"))) %>% 
                   dplyr::group_by(Group) %>% 
                   dplyr::mutate(Date      = unique_dates, 
                                 Obs_cases = cumsum(New_cases * !duplicated(New_cases)),
                                 median    = cumsum(median * !duplicated(median)),
                                 low       = cumsum(low    * !duplicated(low)),
                                 high      = cumsum(high   * !duplicated(high)),
                                 low25     = cumsum(low25  * !duplicated(low25)),
                                 high75    = cumsum(high75 * !duplicated(high75)) )

REACT2_adjusted_cusum_infections <- aggregate_react_cumcases(x           = data_all$Age_distribution, 
                                                             user_AgeGrp = age_mapping_UN) 
colnames(REACT2_adjusted_cusum_infections)[1] <- "Group"

#---- Plot the age-stratified cumulative infections:
ggplot(cusum_data_plot,
       aes(Date      = Date,
           New_Cases = median)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right",
             shrink = T) +
  geom_point(aes(x   = Date,
                 y   = Obs_cases,
                 fill= "Cumulative reported infections"
  )) +
  geom_line(aes(x   = Date, 
                y   = Obs_cases,
                color= "Cumulative reported infections")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median Estimated cumulative infections"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI Estimated cumulative infections"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI Estimated cumulative infections"),
              alpha = 0.5)  +
  labs(x = "Epidemiological Date",
       y = "Cumulative infections") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  scale_fill_manual(values = c('Cumulative reported infections' = "black")) +
  scale_colour_manual(name = '', 
                      values = c("Median Estimated cumulative infections"  = "black", 
                                 "50% CrI Estimated cumulative infections" = "gray20",
                                 "95% CrI Estimated cumulative infections" = "gray40" 
                      )
  ) +
  geom_hline(data = REACT2_adjusted_cusum_infections,
             aes(yintercept = Lower),
             colour   = "red", 
             linetype = "dashed") +
  geom_hline(data = REACT2_adjusted_cusum_infections,
             aes(yintercept = Upper),
             colour   = "red", 
             linetype = "dashed") +
  geom_hline(data = REACT2_adjusted_cusum_infections,
             aes(yintercept = Mean),
             colour   = "red", 
             linetype = "solid") +
  geom_vline(xintercept = as.Date("2020-07-17"), colour = "blue" )  +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.box       = "vertical",
        legend.title     = element_blank() ) 

###########################################################################
#
# Prior of the contact matrix vs the posterior - Feedback from infection 
# data to the contacts.
#
###########################################################################
posts_1      <- rstan::extract(nuts_fit_1)          
niters       <- nrow(posts_1[["cm_sample"]])
prior_cm     <- array(numeric(),c(niters, cov_data$A, cov_data$A)) 
plot_cm_list <- list()
plot_cm_indx <- 1

for (i in 1:cov_data$A){
  for (j in 1:cov_data$A){
    
    # Generate draws from the prior distribution:
    for(k in 1:niters){
      
      L_raw    <- matrix(0, nrow = age_bands, ncol = age_bands)
      L_sample <- L_raw
      
      for(col in 1:nrow(L_raw)){
        for(row in col:nrow(L_raw)){
          L_raw[row,col]    <- rnorm(1, 0, 1)
          L_sample[row,col] <- L_cm[row,col] + (cov_data$p_sigmaCM*L_cm[row,col]) * L_raw[row,col]
        }
      }
      
      # Revert to cm_sample:
      cm_sample <- solve( diag(aggr_age$PopTotal) ) %*% L_sample %*% t(L_sample)
      
      prior_cm[k, i, j] <- cm_sample[i, j]
      
    }# End for
    
    # Create the appropriate dataset in a format acceptable by ggplot:
    dt_cm_post <- data.frame(Prior     = prior_cm[,i, j],
                             Posterior = posts_1[["cm_sample"]][,i, j])
    dt_cm_post2 <- reshape2::melt(dt_cm_post)
    
    p <- ggplot(dt_cm_post2, 
                aes(x    = value, 
                    fill = variable)) +
      geom_density(alpha = 0.8) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
      
      
      # scale_x_continuous(
      #   limits = c(min(dt_cm_post2$value)*0.95, max(dt_cm_post2$value)*1.05)
      # ) +
      labs(x = "",
           y = "Density",
           title = paste0("C[",i,",",j,"]")) +
      theme(strip.placement  = "outside",
            strip.background = element_rect(fill = NA, colour="grey50"),
            panel.spacing    = unit(0,"cm"),
            legend.position  = "bottom",
            legend.title     = element_blank() )
    
    plot_cm_list[[plot_cm_indx]] <- p
    
    # Increment the index of stored graph:
    plot_cm_indx <- plot_cm_indx + 1
    
  }# End for
}# End for

#--- Remove redundant objects:
rm(i, dt_cm_post, dt_cm_post2)

cm_posterior_store_graph <- grid_arrange_shared_legend(plot_cm_list[1:(cov_data$A^2)], 
                                                       nrow = cov_data$A)
  
###########################################################################
#
# Effective contact rate:
#
###########################################################################
model <- "Multi-BM"

if (model == "Multi-BM"){
  plot_betat <- 
    ggplot(data_eff$GBM) +
    facet_wrap(. ~ Group, 
               scales = "free_y", 
               ncol   = 1,
               strip.position = "right") +
    geom_line(aes(x     = Date,
                  y     = Beta_median,
                  color = "Median"),
              size  = 1.3) +
    geom_ribbon(aes(x    = Date,
                    ymin = Beta_low25,
                    ymax = Beta_high75,
                    fill = "50% CrI"),
                alpha = 0.6) +
    labs(x = "Epidemiological Date",
         y = "Effective contact rate") +
    scale_x_date(date_breaks       = "1 month", 
                 date_minor_breaks = "1 month") + 
    scale_fill_manual(values = c("50% CrI" = "gray40")) +
    scale_colour_manual(name   = '',
                        values = c('Median' = "black")) +
    theme_bw() +
    theme(panel.spacing    = unit(0.2,"cm"),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          axis.title.x     = element_text(size = 14, face = "bold"),
          axis.title.y     = element_text(size = 14, face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_blank() ) 
  
} else {

  plot_betat <- 
    ggplot(data_eff$GBM) +
    geom_line(aes(x     = Date,
                  y     = Beta_median,
                  color = "Median"),
              size  = 1.3) +
    geom_ribbon(aes(x    = Date,
                    ymin = Beta_low25,
                    ymax = Beta_high75,
                    fill = "50% CrI"),
                alpha = 0.6) +
    labs(x = "Epidemiological Date",
         y = "Effective contact rate") +
    scale_x_date(date_breaks       = "1 month", 
                 date_minor_breaks = "1 month") + 
    scale_y_continuous(
      limits = c(0,   max(data_eff$GBM$Beta_high75)*1.2),    
      breaks = seq(0, max(data_eff$GBM$Beta_high75)*1.2, 0.05) 
    ) +
    scale_fill_manual(values = c("50% CrI" = "gray40"
    )) +
    scale_colour_manual(name   = '',
                        values = c('Median' = "black")) + 
    theme_bw() +
    theme(panel.spacing    = unit(0.2,"cm"),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          axis.title.x     = element_text(size = 14, face = "bold"),
          axis.title.y     = element_text(size = 14, face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_blank() )  

}# End if

###########################################################################
#
# Age-stratified transmission rate:
#
###########################################################################
plot_transrate <- 
  ggplot(data_eff$Transmission_rate) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size  = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.6) +
  labs(x = "Epidemiological Date",
       y = "Transmission rate") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_fill_manual(values = c("50% CrI" = "gray40"
  )) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +
  theme_bw() +
  theme(panel.spacing    = unit(0.25,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank() )   
  
###########################################################################
#
# Effective reproduction number 
#
###########################################################################
ggplot(data_eff$R_eff_t) +
  geom_line(aes(x     = Date,
  			y     = eff_rt_median,
  			color = "Median"),
  		size  = 1.3) +
  geom_ribbon(aes(x    = Date,
  			  ymin = eff_rt_low25,
  			  ymax = eff_rt_high75,
  			  fill = "50% CrI"),
  		  alpha = 0.6) +
  geom_hline(yintercept = 1, color = "black") +
  labs(x = "Epidemiological Date",
     y = "Effective reproduction number") +
  scale_x_date(date_breaks = "1 month", 
  		   date_minor_breaks = "1 month") + 
  scale_y_continuous(limits = c(0,   max(data_eff$R_eff_t$eff_rt_high75)*1.1), 
  				 breaks = c(seq(0, max(data_eff$R_eff_t$eff_rt_high75)*1.1, 0.2)) ) +
  scale_fill_manual(values = c("50% CrI" = "gray40"
  ) ) +
  scale_colour_manual(name   = '',
  				  values = c('Median' = "black")) +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
      	axis.text.x      = element_text(angle = 45, hjust = 1),
      	axis.title.x     = element_text(size = 14, face = "bold"),
      	axis.title.y     = element_text(size = 14, face = "bold"),
      	legend.position  = "bottom",
      	legend.title     = element_blank() ) +
  geom_vline(data    = data_npi_long,
  		 mapping = aes(xintercept = Date,
  					   color   = "black"), 
  		 key_glyph = glyph_vline,
  		 linetype  = 2) +
  ggrepel::geom_text_repel(data = data_npi_long,
  					   aes(x     = Date,
  						   y     = 1.8,
  						   label = NPI),
  					   color        = "black",
  					   force_pull   = 0,
  					   nudge_y      = 0.1,
  					   direction    = "x",
  					   angle        = 90,
  					   hjust        = 0.5,
  					   segment.size = 0.2,
  					   max.iter     = 1e4,
  					   max.time     = 1,
  					   inherit.aes  = FALSE)

###########################################################################
#
# Age-stratified reporting ratio:
#
###########################################################################
ggplot(data_rep_ratio,
       aes(x = Date,
           y = median)) +
  stat_smooth(method  = "gam",
              formula = y ~ s(x),
              n       = nrow(data_rep_ratio)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  labs(x = "Epidemiological Date",
       y = "Estimated daily reporting ratio") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L)) +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank() )

###########################################################################
#
# Store graphs
#
###########################################################################

# Size: 7x12

#---- Model fit to observed deaths:
(fit_deaths_age | fit_deaths_all) + plot_annotation(tag_levels =  "A")

#---- Estimated infections:

if (country == "Austria"){
  # Results for Austria - Age-stratified plots together and aggregated infections individually:
  (age_specific_cases_data_grp1_plot | age_specific_cases_data_grp2_plot) + plot_annotation(tag_levels =  "A")

  # Save the plot of aggregated infections separately:
  fit_cases_all
  
} else if (country == "Greece"){

  ( (age_specific_cases_data_grp1_plot / age_specific_cases_data_grp2_plot) | 
    (age_specific_cases_data_grp3_plot / fit_cases_all) )+ 
    plot_annotation(tag_levels =  "A")
  
}# End if

#---- Transmissibility, transmission rate, contact matrix:
( (plot_betat / plot_transrate) | cm_posterior_store_graph ) + plot_annotation(tag_levels =  "A")
