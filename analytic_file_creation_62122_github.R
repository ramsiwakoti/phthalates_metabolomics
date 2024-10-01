# Author: Ram Siwakoti
# Date: 1/6/2022
# Objective: This program will merge phthalates data with metabolomics data
# Date: 1/6/2022
# Revised: 6/21/22 to include DEHP metabolites in the table. Also changed lines 317, 324, 331

library(readxl)
library(dplyr)
library(tidyselect)
library(tidyr)
library(stringr)
conflicted::conflicts_prefer(dplyr::filter)
library(ggplot2)

#Set up working directory
setwd("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/data_62122")

#Import the data
sg <- read.csv(file = "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/original_data/-1_v54_006_specific_gravity_20200616100219.csv" )
phthalates <- read.csv(file = "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/original_data/2_v54_006_Phthalates_no_dup_03092021_20210625153132.csv")
metabolites <- read_excel(path = "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/original_data/EX00864_INTEGRATED_REPORT_20190108_165901.xlsx", sheet = "Integrated report")
covars <- read.csv(file = "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/original_data/covars.csv")
metabolites_untar_pos <- read_excel(path = "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/original_data/EX00864_INTEGRATED_REPORT_20190108_165901.xlsx", sheet = "Untargeted positive ANNOTATIONS")
metabolites_untar_neg <- read_excel(path = "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/original_data/EX00864_INTEGRATED_REPORT_20190108_165901.xlsx", sheet = "Untargeted negative ANNOTATIONS")

#Put sample IDs from metabolites file in rows as opposed to column
metabolites2 <- metabolites %>% dplyr::select(starts_with("POOLED") | ends_with("v3")|starts_with("Compound Name"))

#Transpose the data file
metabolites_t <- as.data.frame(t(metabolites2[, -ncol(metabolites2)]))

#Get the column names and assign them to the transposed data file
colname <- as.vector(metabolites2$`Compound Name`)
colnames(metabolites_t) <- colname

#Now add sample IDs to the transposed data file.
metabolites_t <- cbind.data.frame(samp = t(colnames(metabolites2))[-110], metabolites_t)

#Create StudyID and VisitID
metabolites_t$StudyID <- substr(metabolites_t$samp, 1,4)
metabolites_t$VISITID <- substr(metabolites_t$samp, 7,7)

metabolites_t <- metabolites_t %>% relocate(StudyID, .after = samp) %>% relocate(VISITID, .after = StudyID)

#Check the confidence level of individual metabolites 
metabolites_untar_pos2 <- metabolites_untar_pos %>% select(`Compound Name`, `ID confidence`) %>% 
  filter(!is.na(`Compound Name`)) %>% distinct() %>%
  mutate(assay = "Untargeted positive")

metabolites_untar_neg2 <- metabolites_untar_neg %>% select(`Compound Name`, `ID confidence`) %>% 
  filter(!is.na(`Compound Name`)) %>% distinct() %>%
  mutate(assay = "Untargeted negative")

metabolites_untar_both <- rbind.data.frame(metabolites_untar_pos2, metabolites_untar_neg2) %>% distinct()

# Remove StudyID = "POOL" from the metabolites_t table

#Add .1 to duplicate columns (added on 8/23/23 because I was getting an error)

# Function to add ".1" to the second occurrence of a column name
make_unique <- function(x) {
  # Find out which names are duplicated
  dupes <- duplicated(x)
  # Add ".1" to the second occurrence
  x[dupes] <- paste0(x[dupes], ".1")
  return(x)
}

# Apply the function to the column names of your data frame
names(metabolites_t) <- make_unique(names(metabolites_t))

metabolites_t <- metabolites_t %>% dplyr::filter(StudyID != "POOL")

#Put phthalates metabolites from phthalates data in a wide format

#uniq phthalates
phthalates_names <- phthalates %>% dplyr::select(ANALYTE, ANALYTE_FULL) %>% distinct()

phthalates2 <- phthalates %>% dplyr::select(StudyID, VISITID, CONC, CONCUNITS, BLOD, LOD, ANALYTE)

#Only select samples present in metabolomic data
phthalates2$StudyID <- as.character(phthalates2$StudyID)

#There is a single instance where studyid is just 1 number (studyid = 2). Lets add zero before it to be consistent with the length
#of studyid in metabolites file 
phthalates2$StudyID <- str_pad(phthalates2$StudyID, 4, pad = "0")

phthalates3 <- phthalates2 %>% inner_join(metabolites_t, by = c("StudyID")) %>% 
  dplyr::select(StudyID, VISITID.x, CONC, CONCUNITS, BLOD, LOD, ANALYTE) %>% rename("VISITID" = "VISITID.x")

#Count number of samples below LOD for each chemical for visit 3
blod_summary1 <- phthalates3 %>% group_by(ANALYTE, VISITID) %>% dplyr::summarize(detrate = 100-sum(BLOD)/99*100) %>% dplyr::filter(VISITID == 1)
blod_summary2 <- phthalates3 %>% group_by(ANALYTE, VISITID) %>% dplyr::summarize(detrate = 100-sum(BLOD)/99*100) %>% dplyr::filter(VISITID == 2)
blod_summary3 <- phthalates3 %>% group_by(ANALYTE, VISITID) %>% dplyr::summarize(detrate = 100-sum(BLOD)/99*100) %>% dplyr::filter(VISITID == 3)

blod_summary <- rbind.data.frame(blod_summary1, blod_summary2, blod_summary3) 

blod_summary_chk <- phthalates3 %>% group_by(ANALYTE) %>% dplyr::summarize(detrate = 100-sum(BLOD)/297*100, nasum = sum(is.na(BLOD))/297*100) 

ggplot(blod_summary, aes(y = ANALYTE, x = detrate/3, fill = ANALYTE)) + 
  geom_bar(stat = "identity") + theme(legend.position = "none") + xlab("det rate") + ylab("Phthalates") +
  ggtitle("Number of detects for each Phthalate metoabolite")

#So, MHINCH, MNP, MCOCH, and MEHP have more than 20% of values below the LOD. 

#Number of obs with unique phthalates samples 
temp <- phthalates2 %>% dplyr::select(StudyID) %>% distinct()

#Number of phthalates per sample per visit
phthalates_sum <- phthalates3 %>% group_by(StudyID, VISITID) %>% summarise(numofanalyte = n())
table(phthalates_sum$numofanalyte)

#Convert data to wide format
phthalates3_wide <- phthalates3 %>% select (StudyID, VISITID, CONC, CONCUNITS, ANALYTE) %>% 
  pivot_wider(names_from = ANALYTE, values_from = CONC)

#Calculate the total molar amount of parent phthalates - ng/ml * mol/gm = nmol/ml 

#MW of phthalates
MEHPmw <- 278.348
MEHHPmw <- 294.34
MEOHPmw <- 292.33
MECPPmw <- 308.3264
MBPmw <- 222.24
MHBPmw <- 238.24
MIBPmw <- 222.24
MHIBPmw <- 238.24
MCOCHmw <- 328.4 #https://pubchem.ncbi.nlm.nih.gov/compound/145459122
MHINCHmw <- 314.4 #https://pubchem.ncbi.nlm.nih.gov/compound/129853628
MECPTPmw <- 308.33 #https://pubchem.ncbi.nlm.nih.gov/compound/145459123
MEHHTPmw <- 294.34 #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6917904/
MNPmw <- 292.4 #https://pubchem.ncbi.nlm.nih.gov/compound/Monoisononyl-phthalate
MCOPmw <- 322.357 #http://exposome-explorer.iarc.fr/compounds/1455
DEHPmw <- (MEHPmw + MEHHPmw +MEOHPmw +MECPPmw)/4
DBPmw <- (MBPmw +MHBPmw)/2
DIBPmw <- (MIBPmw + MHIBPmw)/2
DINCHmw <- (MCOCHmw + MHINCHmw)/2
DEHTPmw <- (MECPTPmw + MEHHTPmw)/2

phthalates3_wide <- phthalates3_wide %>% mutate(molMEHP = MEHP/MEHPmw, molMEOHP = MEOHP/MEOHPmw, molMECPP = MECPP/MECPPmw, 
                                                molMEHHP = MEHHP/MEHHPmw, molMBP = MBP/MBPmw, molMHBP = MHBP/MHBPmw,
                                                molMIBP = MIBP/MIBPmw, molMHIBP = MHIBP/MHIBPmw,
                                                molMCOCH = MCOCH/MCOCHmw, molMHINCH = MHINCH/MHINCHmw,
                                                molMECPTP = MECPTP/MECPTPmw, molMEHHTP = MEHHTP/MEHHTPmw)

phthalates3_wide <- phthalates3_wide %>% rowwise %>% mutate(DEHP = ifelse(sum(is.na(c(MEHP, MEHHP, MEOHP, MECPP))) == 4, NA, 
                                                                          sum(molMEHP, molMEHHP, molMEOHP, molMECPP, na.rm = TRUE)*DEHPmw),
                                                            DBP = ifelse(sum(is.na(c(MBP, MHBP))) == 2, NA, sum(molMBP, molMHBP, na.rm = TRUE)*DBPmw),
                                                            DIBP = ifelse(sum(is.na(c(MIBP, MHIBP))) == 2, NA, sum(molMIBP, molMHIBP, na.rm = TRUE)*DIBPmw),
                                                            DINCH = ifelse(sum(is.na(c(MCOCH, MHINCH))) == 2, NA, sum(molMCOCH, molMHINCH, na.rm = TRUE)*DINCHmw),
                                                            DEHTP = ifelse(sum(is.na(c(MECPTP, MEHHTP))) == 2, NA, sum(molMECPTP,  molMEHHTP, na.rm = TRUE)*DEHTPmw))

#Now only keep the variables we want 
phthalates3_wide_chk <- phthalates3_wide %>% select(StudyID, VISITID, CONCUNITS, DEHP, DBP, DIBP, DINCH, DEHTP, MEHP, 
                                                    MEHHP, MEOHP, MECPP, MBP, MHBP, MIBP, MHIBP, MCOCH, MHINCH, MECPTP, MEHHTP, molMHINCH, molMCOCH)

phthalates3_wide_chk2 <- phthalates3_wide
phthalates3_wide_chk2$VISITID <- as.character(phthalates3_wide_chk2$VISITID)

#To include individual metabolites of DEHP
phthalates3_wide <- phthalates3_wide %>% select(StudyID, VISITID, CONCUNITS, DEHP, MEHP, MEHHP, MEOHP, MECPP, DBP, MBP, MHBP, DIBP, MHIBP, MIBP, DINCH, DEHTP, MBZP, MEP, MCPP, MCNP, MNP, MCOP, MONP)
phthalates3_wide$VISITID <- as.character(phthalates3_wide$VISITID)

#Merge sg and metabolites
sg <- sg %>% rename("StudyID" = "studyid") 
sg <- sg %>% mutate(StudyID = as.character(StudyID))

#Similar to phthalates, lets add zeros infront of studyid = 2
sg$StudyID = str_pad(sg$StudyID, 4, pad = "0")

sg2 <- sg %>% inner_join(metabolites_t, by = c("StudyID")) %>% select(StudyID, specificgravity_v1, specificgravity_v2, specificgravity_v3)

#There are some missing values for sg - replace missing values of sg with the mean of available sgs for that sample 
sg_t <- t(sg2)

#replace NAs with the mean
for(i in 1:ncol(sg_t)) {
  sg_t[ , i][is.na(sg_t[ , i])] <- mean(as.numeric(sg_t[2:nrow(sg_t) , i]), na.rm=TRUE)
}

#transpose sg_t now
sg3 <- as.data.frame(t(sg_t))
sg3$StudyID <- as.character(sg3$StudyID)

#Convert all sgs from character to numeric
for(i in 2:ncol(sg3)){
  sg3[,i] = as.numeric(sg3[,i])
}

#Now convert sg data to long format so that it can be merged with phthalates data by StudyID and visit
sg3_long <- sg3 %>% pivot_longer(specificgravity_v1:specificgravity_v3, names_to = "visit", values_to = "sg") %>% 
            mutate(VISITID = substr(visit,18, 18)) %>% select(!visit)

#Merge sg data with phthalates data
phthalates4_wide <- phthalates3_wide %>% left_join(sg3_long, by = c("StudyID", "VISITID"))
phthalates4_wide_chk <- phthalates3_wide_chk2 %>% left_join(sg3_long, by = c("StudyID", "VISITID"))

#Now I want to adjust individual phthalates for specific gravity using the formula used in Rodriguez-Carmona et al. 
#Pc = P(Median SG-1)/(SG-1)

phthalates4_wide_copy = phthalates4_wide

#Divide phthalates4_wide by VISITID, apply the above function, and then recombine them again
phthalates4_wide_v1 <- phthalates4_wide %>% filter(VISITID == 1)
phthalates4_wide_v2 <- phthalates4_wide %>% filter(VISITID == 2)
phthalates4_wide_v3 <- phthalates4_wide %>% filter(VISITID == 3)

#For visit = 1
for(i in 4:19) {
  phthalates4_wide_v1[ , i] <- (phthalates4_wide_v1[ , i])*(1.019-1)/(phthalates4_wide_v1$sg-1)
}

#For visit = 2
for(i in 4:19) {
  phthalates4_wide_v2[ , i] <- (phthalates4_wide_v2[ , i])*(1.019-1)/(phthalates4_wide_v2$sg-1)
}

#For visit = 3
for(i in 4:19) {
  phthalates4_wide_v3[ , i] <- (phthalates4_wide_v3[ , i])*(1.02-1)/(phthalates4_wide_v3$sg-1)
}

#Now combine all three visits together
phthalates4_wide2 <- rbind.data.frame(phthalates4_wide_v1, phthalates4_wide_v2, phthalates4_wide_v3)
phthalates4_wide2_chk <- phthalates4_wide2 %>% select(StudyID) %>% distinct()

#ICC calculation
dehp <- phthalates4_wide2 %>% 
  select(StudyID, VISITID, DEHP) %>% mutate(DEHP = log(DEHP)) %>%
  pivot_wider(values_from = "DEHP", names_from = "VISITID")
psych::ICC(dehp[,-1])

mcpp <- phthalates4_wide2 %>% 
  select(StudyID, VISITID, MCPP) %>% mutate(MCPP = log(MCPP)) %>% 
  pivot_wider(values_from = "MCPP", names_from = "VISITID")
psych::ICC(mcpp[,-1])

mcpp2 <- lme4::lmer(log(MCPP) ~ (1|StudyID), data = phthalates4_wide2)
performance::icc(mcpp2)

mep <- phthalates4_wide2 %>% 
  select(StudyID, VISITID, MEP) %>% 
  pivot_wider(values_from = "MEP", names_from = "VISITID")
psych::ICC(mep[,-1])

#Sum of NAs
na_counts_by_visit <- with(phthalates4_wide2, tapply(DEHP, VISITID, function(x) sum(is.na(x))))
for (colname in c("MEHP", "MEHHP", "MEOHP", "MECPP", "DBP", "MBP", "MHBP", "DIBP", "MHIBP", "MIBP", "DINCH", "DEHTP", "MBZP", "MEP", "MCPP", "MCNP", "MNP", "MCOP", "MONP")) {
  na_counts_by_visit <- cbind(na_counts_by_visit, with(phthalates4_wide2, tapply(get(colname), VISITID, function(x) sum(is.na(x)))))
}
colnames(na_counts_by_visit) <- c("DEHP", "MEHP", "MEHHP", "MEOHP", "MECPP", "DBP", "MBP", "MHBP", "DIBP", "MHIBP", "MIBP", "DINCH", "DEHTP", "MBZP", "MEP", "MCPP", "MCNP", "MNP", "MCOP", "MONP")

print(na_counts_by_visit)

#Lot of missing values for:
#MHBP
#MHIBP
#DINCH
#DEHTP
#MNP
#MONP

#These phthalates will be excluded from analyses. 

#Calculate GM per person per phthalate

# Define a function to compute geometric mean
geometric_mean <- function(x) {
  valid_x <- x[!is.na(x) & x > 0]  # Excluding NAs and non-positive values
  exp(mean(log(valid_x)))
}

# Calculate geometric mean per person
gm_list <- lapply(phthalates4_wide2[, -c(1:3, ncol(phthalates4_wide2))], function(col) {
  tapply(col, phthalates4_wide2$StudyID, geometric_mean)
})

# Convert list to data.frame
gm_per_person_df <- as.data.frame(gm_list)
gm_per_person_df$StudyID <- rownames(gm_per_person_df)
rownames(gm_per_person_df) <- NULL

# If you want 'StudyID' as the first column
gm_per_person_df <- gm_per_person_df[, c(ncol(gm_per_person_df), 1:(ncol(gm_per_person_df)-1))]

#Count number of NAs or NaNs for each phthalate 
na_nan_counts <- sapply(gm_per_person_df, function(col) sum(is.na(col) | is.nan(col)))
print(na_nan_counts)


#Consequently, we will exclude these phthalates from analysis
phthalates4_wide2 <- phthalates4_wide2 %>% select(-c(MHBP, MHIBP, MNP, DINCH, DEHTP, MONP, DBP, DIBP))

#Now combine metabolites data with phthalates data
phthalates_met <- metabolites_t %>% inner_join(phthalates4_wide2, by = c("StudyID"))

#Add covariates to the metabolomics and phthlalates data
covars <- covars %>% rename("StudyID" = "studyid")
covars$StudyID = str_pad(covars$StudyID, 4, pad = "0")
covars$StudyID <- as.character(covars$StudyID)

studyid <- phthalates_met %>% select(StudyID)

covars2 <- covars %>% inner_join(studyid, by = c("StudyID")) %>% distinct()

covars2_chk <- covars %>% anti_join(studyid, by = "StudyID")

phthalates_met <- phthalates_met %>% inner_join(covars2, by = c("StudyID"))

#Save the files 
# saveRDS(metabolites_t, "metabolites_t.rds")
# saveRDS(phthalates4_wide2, "phthalates4_wide2.rds")
# saveRDS(sg3, "sg3.rds")
# saveRDS(phthalates_met, "phthalates_met.rds")
# saveRDS(blod_summary, "blod_summary.rds")
# saveRDS(covars2, "covars2.rds")

rm(list = ls())