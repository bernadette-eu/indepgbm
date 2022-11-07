#################################################################################
#
# Guidelines: https://timriffe.github.io/covid_age/GettingStarted.html
#
#################################################################################

#---- This downloads and reads in one line:
dest      <- "...//"
data_path <- "...//"

cdb_repo <- osf_retrieve_node("mpwjq")

setwd(dest)
osf_ls_files(cdb_repo, path = "Data") %>%
    dplyr::filter(name == "Output_5.zip") %>%
    osf_download(conflicts = "overwrite")

Output_5 <- read_csv("Output_5.zip",
                     skip      = 3,
                     col_types = "ccccciiddd")

#---- Convert to date class:
Austria_5yo_bands_data <- Output_5 %>%
                          dplyr::filter(Country == "Austria",
                                        Region == "All",
                                        Sex    == "b",
                                        !is.na(Cases)
                                        ) %>%
                          dplyr::mutate(Date = dmy(Date),
                                 AgeInt = Age + 4) %>%
                          dplyr::rename(Group = Age) %>%
                          dplyr::select(c(Date, Group, AgeInt, Cases)) %>%
                          dplyr::mutate(Group = paste(Group, AgeInt, sep = "-") )%>%
                          dplyr::select(c(Date, Group, Cases))

c(min(Austria_5yo_bands_data$Date), max(Austria_5yo_bands_data$Date))

save(Austria_5yo_bands_data,
     file = paste0(data_path, "Austria_5yo_bands_data", ".RData"))

#---- Load the dataset:
load(file = paste0(data_path, "Austria_5yo_bands_data", ".RData"))