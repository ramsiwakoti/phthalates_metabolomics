#######################################################
#######################################################
## This script is used to publish some of the tables ##
## and figures for metabolomics paper                ##
#######################################################
#######################################################



##Load the libraries 
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(tidyverse)



##Import the datasets
setwd("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/")
load("phthalates_met2_b.RData")
load("phthalates_met2.RData")
blod_summary_prenatal <- readRDS("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/blod_summary.rds") #prenatal
phthalates4_wide2_allphthalates <- readRDS("phthalates4_wide2_allphthalates.rds")


phthalates_met2_b <- phthalates_met2_b %>% select(StudyID, DINCH, DEHTP, MNP, MONP) %>% mutate(across(2:5, exp))

phthalates_met2 <- phthalates_met2 %>% inner_join(phthalates_met2_b, by = c("StudyID"))


##Make table 1

phthalates_met2_covars <- phthalates_met2 %>% select(isage, prebmi, PPSEX)


table1::table1(~isage + prebmi + PPSEX, data = phthalates_met2_covars)
table(phthalates_met2_covars$PPSEX)


# table1::table1(~ISAGE + PPSEX + BW, phthalates_met2_cord)



## Make table 2 - exposure distributions 

#Prenatal data

prenatal_dat_long <- phthalates_met2 %>% 
  select(StudyID, MCOP, MBZP, MCNP, MCPP, MEP, MNP, MONP, DBP, DEHP, DIBP, DINCH, DEHTP) %>%
  pivot_longer(cols = MCOP:DEHTP, names_to = "Phthalates", values_to = "Phthalates_conc") %>%
  mutate(Phthalates_conc = exp(Phthalates_conc))
  
geomean <- function(x){
  gm  = exp(mean(log(x), na.rm = TRUE))
  return(gm)
}

sum_prenatal_dat <- prenatal_dat_long %>% group_by(Phthalates) %>%
  summarise(GM = round(geomean(Phthalates_conc), 2), 
            median = round(quantile(Phthalates_conc, 0.5, na.rm = TRUE), 2),
            median_form = formatC(median, format = "f", digits = 2),
            q1 = round(quantile(Phthalates_conc, 0.25, na.rm = TRUE), 2), 
            q1_form = formatC(q1, format = "f", digits = 2),
            q3 = round(quantile(Phthalates_conc, 0.75, na.rm = TRUE), 2), 
            q3_form = formatC(q3, format = "f", digits = 2),
            min = round(min(Phthalates_conc, na.rm = TRUE), 2), 
            min_form = formatC(min, format = "f", digits = 2), 
            max = round(max(Phthalates_conc, na.rm = TRUE), 2),
            max_form = formatC(max, format = "f", digits = 2), 
            med_minmax = paste0(median, " ", "(", min_form, ",", " ", max_form, ")")) 

max_length <- max(nchar(sum_prenatal_dat$med_minmax))
table(phthalates_met2$PPSEX)


sum_prenatal_dat$med_minmax2 = str_pad(sum_prenatal_dat$med_minmax, max_length, side = "right")
write.csv(sum_prenatal_dat, "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/manuscript_writing/prenatal_exp_summary.csv")


################################################################################################################################################
#Prenatal data with all phthalates

#First get the gestational average exposure 

phthalates4_wide2_allphthalates2<- phthalates4_wide2_allphthalates[,c(1, 4:22)]

conflicted::conflicts_prefer(dplyr::summarize)

#NAs per phthalates stratified by study visit

na_count_by_visit <- phthalates4_wide2_allphthalates %>%
  select(VISITID, 4:22) %>%
  group_by(VISITID) %>%
  summarize(across(everything(), ~sum(is.na(.)), .names = "na_count_{.col}"))

# geomean  <- function(x) {
#   # Exclude NA values
#   x <- x[!is.na(x)]
#   if(length(x) == 0) return(NA)
#   # Compute geometric mean
#   return(exp(mean(log(x))))
# }

  
geomean <- function(x){
  gm  = exp(mean(log(x), na.rm = TRUE))
  return(gm)
}

phthalates4_wide2_allphthalates_gm <- phthalates4_wide2_allphthalates2 %>%
  group_by(StudyID) %>%
  summarise(across(MNP:MONP, geomean, .names = "GM_{.col}"))



prenatal_dat_long2 <- phthalates4_wide2_allphthalates_gm %>% 
  pivot_longer(cols = GM_MNP:GM_MONP, names_to = "Phthalates", values_to = "Phthalates_conc") 


sum_prenatal_dat2 <- prenatal_dat_long2 %>% group_by(Phthalates) %>%
  summarise(GM = round(geomean(Phthalates_conc), 2), 
            median = round(quantile(Phthalates_conc, 0.5, na.rm = TRUE), 2),
            median_form = formatC(median, format = "f", digits = 2),
            q1 = round(quantile(Phthalates_conc, 0.25, na.rm = TRUE), 2), 
            q1_form = formatC(q1, format = "f", digits = 2),
            q3 = round(quantile(Phthalates_conc, 0.75, na.rm = TRUE), 2), 
            q3_form = formatC(q3, format = "f", digits = 2),
            min = round(min(Phthalates_conc, na.rm = TRUE), 2), 
            min_form = formatC(min, format = "f", digits = 2), 
            max = round(max(Phthalates_conc, na.rm = TRUE), 2),
            max_form = formatC(max, format = "f", digits = 2), 
            med_minmax = paste0(median, " ", "(", min_form, ",", " ", max_form, ")"),
            med_q1q3 = paste0(median, " ", "(", q1_form, ",", " ", q3_form, ")"),
            
            sum_na = sum(is.na(Phthalates_conc))) 

max_length2 <- max(nchar(sum_prenatal_dat2$med_minmax))


sum_prenatal_dat2$med_minmax2 = str_pad(sum_prenatal_dat2$med_minmax, max_length2, side = "right")
sum_prenatal_dat2$med_q1q3_2 = str_pad(sum_prenatal_dat2$med_q1q3, max_length2, side = "right")
sum_prenatal_dat2 <- sum_prenatal_dat2 %>% mutate(Phthalates2 = substr(Phthalates, 4, nchar(Phthalates)))


#Look at proportion below LOD

blod_summary_prenatal_wide <- blod_summary_prenatal %>% pivot_wider(values_from = nlt_bld, names_from = VISITID)
names(blod_summary_prenatal_wide) <- c("Phthalates2", "Visit1", "Visit2", "Visit3")

sum_prenatal_dat2 <- sum_prenatal_dat2 %>% inner_join(blod_summary_prenatal_wide, by = "Phthalates2")
write.csv(sum_prenatal_dat2, "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/manuscript_writing/prenatal_exp_summary_allphthalates.csv")

####################################################################################################################################################
#Cord blood data
load("./cord_blood_data/phthalates_met2_cordblood.RData")
phthalates_met2_cord <- phthalates_met2
# 
# phthalates_met2_cord <- phthalates_met2_cord %>%
#   select(-c(SG_MEHP_GM, SG_MEHHP_GM, SG_MEOHP_GM, SG_MECPP_GM, SG_MBP_GM, SG_MHBP_GM, SG_MIBP_GM, SG_MHIBP_GM, SG_MECPTP_GM))

cord_dat_long <- phthalates_met2_cord %>% 
  select(`Study ID`, 298:317) %>%
  pivot_longer(cols = SG_MBP_GM:SG_DEHTP_GM, names_to = "Phthalates", values_to = "Phthalates_conc") 


sum_cord_dat <- cord_dat_long %>% group_by(Phthalates) %>%
  summarise(GM = round(geomean(Phthalates_conc), 2), 
            median = round(quantile(Phthalates_conc, 0.5, na.rm = TRUE), 2),
            median_form = formatC(median, format = "f", digits = 2),
            q1 = round(quantile(Phthalates_conc, 0.25, na.rm = TRUE), 2), 
            q1_form = formatC(q1, format = "f", digits = 2),
            q3 = round(quantile(Phthalates_conc, 0.75, na.rm = TRUE), 2), 
            q3_form = formatC(q3, format = "f", digits = 2),
            min = round(min(Phthalates_conc, na.rm = TRUE), 2), 
            min_form = formatC(min, format = "f", digits = 2), 
            max = round(max(Phthalates_conc, na.rm = TRUE), 2),
            max_form = formatC(max, format = "f", digits = 2), 
            med_minmax = paste0(median, " ", "(", min_form, ",", " ", max_form, ")"),
            med_q1q3 = paste0(median, " ", "(", q1_form, ",", " ", q3_form, ")"),
            
            sum_na = sum(is.na(Phthalates_conc))) %>%
  mutate(phthalates = substr(Phthalates, 1, nchar(Phthalates)-3))

max_length <- max(nchar(sum_cord_dat$med_minmax))
sum_cord_dat$med_minmax2 = str_pad(sum_cord_dat$med_minmax, max_length, side = "right")
write.csv(sum_cord_dat, "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/manuscript_writing/cord_exp_summary.csv")


  

#Figure 1






