# Author: Ram Siwakoti
# Date: 6/23/22
# Objective: This program will merge phthalates data with cord blood metabolomics data

library(readxl)
library(dplyr)
library(tidyselect)
library(tidyr)
library(stringr)
library(data.table)

#Set up working directory
setwd("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/cord_blood_data")

#Import the data
phthalates <- read.csv(file = "phthalates_cordmetab.csv", check.names = FALSE)
metabolites_neg <- read.csv(file = "negative-mode-data.csv", check.names = FALSE, na.strings='.')
metabolites_pos <- read.csv(file = "positive-mode-data.csv", check.names = FALSE, na.strings='.')
metabolites <- read.csv(file = "data-processed-named-288.csv", check.names = FALSE)
idfile <- read_excel(path = "C:/Users/rsiwa/Dropbox (University of Michigan)/Ram_shared/data/cord blood metabolomics/EX01061 - Watkins - Neg mode Final Report.xlsx", 
                     sheet = "Sample Information")
idfile$`Subject ID` <- substr(idfile$`Subject ID`, 1,4)
idfile2 <- idfile %>% select(`Sample ID`,  `Sample Name`, `Subject ID`)

covars2a <- read.csv(file = "C:/Users/rsiwa/Dropbox (University of Michigan)/Ram_shared/data/2_v54_006_deb_hsd_20220914093338_updated.csv") 
covars2b <- read_xlsx(path = "C:/Users/rsiwa/Dropbox (University of Michigan)/Ram_shared/data/Copy of IDs_with_missing_values_cordblood_IAOrev ZRP10.13.22.xlsx", sheet = "sheet1") %>%
  rename("Subject ID" = `Study ID`, "ISAGE_1" = ISAGE, "PPSEX_1" = PPSEX, "PPWEIGHTP_1" = PPWEIGHTP, "PPWEIGHTO_1" = PPWEIGHTO)

covars2a <- covars2a %>% mutate(`Subject ID` = PROTECTID)  %>% 
  select(`Subject ID`, PPSEX, isage, PPWEIGHTP, PPWEIGHTO, WTPREPREG, FVCURRHT_FOOT, FVCURRHT_INCH) %>%
  rename("ISAGE" = isage) %>%
  mutate(height_inch = (as.numeric(FVCURRHT_FOOT))*12 + as.numeric(FVCURRHT_INCH)) %>%
  mutate(prebmi = (WTPREPREG/(height_inch)^2)*703)

covars2 <- covars2a %>% left_join(covars2b, by = "Subject ID")
covars2 <- covars2 %>% mutate(ISAGE_2 = ifelse(is.na(ISAGE), ISAGE_1, ISAGE),
                              PPSEX_2 = ifelse(is.na(PPSEX), PPSEX_1, PPSEX),
                              PPWEIGHTP_2 = ifelse(is.na(PPWEIGHTP), PPWEIGHTP_1, PPWEIGHTP),
                              PPWEIGHTO_2 = ifelse(is.na(PPWEIGHTO), PPWEIGHTO_1, PPWEIGHTO))

chk <- covars2 %>% filter(is.na(PPSEX))

covars2 <- covars2 %>% select(`Subject ID`, PPSEX_2, ISAGE_2, PPWEIGHTP_2, PPWEIGHTO_2, prebmi)
names(covars2) <- c("Subject ID", "PPSEX", "ISAGE", "PPWEIGHTP", "PPWEIGHTO", "prebmi")

visdat::vis_miss(covars2)

#Transpose metabolite files so that I can link them to the column names from the processed file 
t_metabolites_neg <- data.table::transpose(metabolites_neg, keep.names = "Metabolite", make.names = "Sample ID")
t_metabolites_pos <- data.table::transpose(metabolites_pos, keep.names = "Metabolite", make.names = "Sample ID")
metab_name <- as.data.frame(names(metabolites))
colnames(metab_name) <- "Metabolite"

#Since there are six metabolites that are identified via both negative and positive mode, I am going to change their names in each file like Gayatri did. 
t_metabolites_neg <- t_metabolites_neg %>% mutate(Metabolite2 = case_when(Metabolite == "Glutamine" ~ "Glutamine-N",
                                                                          Metabolite == "Isoleucine" ~ "Isoleucine-N",
                                                                          Metabolite == "Leucine/Isoleucine" ~ "Leucine/Isoleucine-N",
                                                                          Metabolite == "Tryptophan" ~ "Tryptophan-N",
                                                                          Metabolite == "Tyrosine" ~ "Tyrosine-N",
                                                                          Metabolite == "Phenylalanine" ~ "Phenylalanine-N"),
                                                  Metabolite3 = ifelse(is.na(Metabolite2), Metabolite, Metabolite2)) %>%
  select(-c(Metabolite, Metabolite2)) %>% rename("Metabolite" = "Metabolite3")

t_metabolites_pos <- t_metabolites_pos %>% mutate(Metabolite2 = case_when(Metabolite == "Glutamine" ~ "Glutamine-P",
                                                                          Metabolite == "Isoleucine" ~ "Isoleucine-P",
                                                                          Metabolite == "Leucine/Isoleucine" ~ "Leucine/Isoleucine-P",
                                                                          Metabolite == "Tryptophan" ~ "Tryptophan-P",
                                                                          Metabolite == "Tyrosine" ~ "Tyrosine-P",
                                                                          Metabolite == "Phenylalanine" ~ "Phenylalanine-P"),
                                                  Metabolite3 = ifelse(is.na(Metabolite2), Metabolite, Metabolite2)) %>%
  select(-c(Metabolite, Metabolite2)) %>% rename("Metabolite" = "Metabolite3")

metabolites_neg2 <- metab_name %>% inner_join(t_metabolites_neg, by = "Metabolite")
metabolites_pos2 <- metab_name %>% inner_join(t_metabolites_pos, by = "Metabolite")

#make sure the name of columns in both tables are the same 
names_long_neg <- as.data.frame(names(metabolites_neg2))
colnames(names_long_neg) <- "Sample ID2"
names_long_neg <- names_long_neg %>% filter(`Sample ID2` != "Metabolite") %>% mutate(`Sample ID` = substr(`Sample ID2`, 30, 38), 
                                                                                     `Sample Name` = substr(`Sample ID2`, 40, 46)) %>% select(-`Sample ID2`)
names_long_neg <- names_long_neg %>% distinct()

names_long_pos<- as.data.frame(names(metabolites_pos2))
colnames(names_long_pos) <- "Sample ID2"
names_long_pos <- names_long_pos %>% filter(`Sample ID2` != "Metabolite") %>% mutate(`Sample ID` = substr(`Sample ID2`, 30, 38), 
                                                                                     `Sample Name` = substr(`Sample ID2`, 40, 46)) %>%  select(-`Sample ID2`)
names_long_pos <- names_long_pos %>% distinct()

#Check again
sum(names_long_pos$`Sample ID` != names_long_neg$`Sample ID`)
sum(names_long_pos$`Sample Name` != names_long_neg$`Sample Name`)

#change column names of metabolites_pos2 and metabolites_neg2 
samp_ids <- c("Metabolite", names_long_pos$`Sample ID`)
names(metabolites_neg2) <- samp_ids
names(metabolites_pos2) <- samp_ids

#Combine neg and pos files 
metabolites2 <- rbind.data.frame(metabolites_neg2, metabolites_pos2)

metabolites_t <- data.table::transpose(metabolites2, keep.names = "Sample ID", make.names = "Metabolite")


#Link metabolites file with the subject ID file
metabolites_t <- metabolites_t %>% inner_join(idfile2, by = "Sample ID") 
metabolites_t <- metabolites_t %>% select(`Sample ID`, `Sample Name`, `Subject ID`, everything()) %>% rename("Study ID" = "Subject ID")

#uniq phthalates
phthalates_names <- names(phthalates)

phthalates <- phthalates %>% rename("Study ID" = "studyid")
str(phthalates)

#Combine covariates, phthalates, and metabolomics data
covars2 <- covars2 %>% rename("Study ID" = "Subject ID")
str(covars2)

#Add covariates to the phthalates data
phthalates_covars <- phthalates %>% inner_join(covars2, by = c("Study ID")) %>% distinct()

covars_chk <- phthalates_covars %>% select(`Study ID`, PPSEX, PPWEIGHTP, PPWEIGHTO, ISAGE, prebmi)
covars_chk_NA <- covars_chk %>% filter(is.na(PPSEX))
write.csv(covars_chk_NA, file = "IDs_with_missing_values_cordblood.csv")

#Add metabolomics
phthalates_covars$`Study ID` <- as.character(phthalates_covars$`Study ID`)
phthalates_met <- metabolites_t %>% inner_join(phthalates_covars, by = "Study ID") %>% distinct()

#Save the files 
saveRDS(metabolites_t, "metabolites_t.rds")                              
saveRDS(phthalates, "phthalates.rds")
saveRDS(phthalates_met, "phthalates_met.rds")
saveRDS(phthalates_covars, "phthalates_covars.rds")