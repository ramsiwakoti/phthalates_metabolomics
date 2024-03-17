
library(dplyr)

setwd("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/manuscript_writing")

#cordblood metabolites

load("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/cord_blood_data/phthalates_met2_cordblood.RData")
cordmet_data <- phthalates_met2

cordmet <- data.frame(met = names(cordmet_data[,4:276]))

#prenatal metabolites
load("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/phthalates_met2.RData")
prenat_data <- phthalates_met2
prenatmet <- data.frame(met = names(prenat_data[,4:326]))

#remove white space
cordmet$met2 <- stringr::str_trim(toupper(cordmet$met))
prenatmet$met2 <- stringr::str_trim(toupper(prenatmet$met))

metboth <- cordmet %>% inner_join(prenatmet, by = "met2")
met_prenatonly <- prenatmet %>% anti_join(cordmet, by = "met2")
met_cordonly <- cordmet %>% anti_join(prenatmet, by = "met2")

nrow(metboth)
nrow(met_prenatonly)
nrow(met_cordonly)


#bring in a file with p values from prenatal analyiss
prenatal_pval <- read.csv("fit_lin34_adjust_long_prenatal.csv") %>% mutate(met2 = stringr::str_trim(toupper(samp)))


#bring in a file with p values from cord blood analyiss

cordmet_pval <- read.csv("fit_lin34_adjust_long_cordblood.csv") %>% mutate(met2 = stringr::str_trim(toupper(Metabolite)), Phthalates = Phthalate) 

#merge common file with prenatal p values file

prenatal_pval_common <- prenatal_pval %>% inner_join(metboth, by = "met2")

#merge prenatal_pval_common file with cord blood file

prenatal_cord_pval_common2 <- prenatal_pval_common %>% inner_join(cordmet_pval, by = c("met2", "Phthalates")) %>%
  select(met2, Phthalates, betas.x, adj_p.x, betas.y, adj_p.y) 
names(prenatal_cord_pval_common2) <- c("Metabolites", "Phthalates", "Prenatal_beta", "Prenatal_pval", "Cord_beta", "Cord_pval")



write.csv(prenatal_cord_pval_common2, "metabolites_common_prenatal_cord.csv")

library(ggplot2)

prenatal_cord_pval_common2 <- prenatal_cord_pval_common2 %>% mutate(prenatal_pind = ifelse(Prenatal_beta < 0, "Neg", "Pos"),
                                                          cord_pind = ifelse(Cord_beta < 0, "Neg", "Pos"))

ggplot(prenatal_cord_pval_common2, aes(x = Prenatal_beta, y = Cord_beta)) + geom_point() + facet_wrap(~Phthalates)
  

#Make a confusion matrix

mcpp_chk <- prenatal_cord_pval_common2 %>% filter(Phthalates == "MCPP")
gmodels::CrossTable(mcpp_chk$prenatal_pind, mcpp_chk$cord_pind)



#Add classes of common metabolic features 

classes <- read.csv("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/cord_blood_data/list_of_cordblood_metabolites2_csv.csv") %>% distinct()
classes <- classes %>% mutate(class = trimws(class), Class_new2 = class) %>% rename("Metabolites" = "name")


metboth2 <- metboth %>% rename("Metabolites" = "met2") %>% select(Metabolites) %>% inner_join(classes, by = "Metabolites")

#Distribution of classes

class_sum_both <- metboth2 %>% group_by(Class_new2) %>% tally()





