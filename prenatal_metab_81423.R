#Date: 8/14/2023

#In this file, I will run main analysis for prenatal metabolomics data for manuscript production. Most of the code 
#will come from phthalates_metabolites_62023.

#revised on 10/3/23 to remove phthalate metabolites with >30% missing data from the analysis 


########################################################################################################################################################################
########################################################################################################################################################################

#Set up a working directory

setwd("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/manuscript_writing/")


#Load some packages 


library(readxl)
library(plyr)
library(tidyselect)
library(visdat) #missing values visualization
library(missForest)
library(ggplot2)
library(tidyr) #converting data from wide to long and vice-versa 
library(directlabels) #to add labels to line plots
library(ggcorrplot) #correlation plots
library(pca3d)
library(gplots)
library(knitr)
library(Hmisc)
library(purrr)
library(dplyr)
library(corrplot)

conflicted::conflict_prefer("summarise", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("group_by", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("mutate", "dplyr")
conflicted::conflict_prefer("arrange", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

########################################################################################################################################################################

#Obtain the analytic file created in phthalates_metabolites_62023. 

#This file has phthalates data, metabolites data, as well as covariates. Phthalates have been log-transformed
#whereas metabolites data have been auto-scaled. 

load("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/phthalates_met2.RData")
classes <- readxl::read_excel("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/significant_metabolites_with_class_annotated.xlsx", sheet = "classes_81523")

phthalates_met2_short <- phthalates_met2 %>% filter(StudyID != "3034")

########################################################################################################################################################################

#Correlation plot

corph <- corrplot(corr = cor(phthalates_met2[,327:334], method = "spearman", use = "complete"), method = "number", type = "lower")




####################################################################################################################################################
#Run MWAS 


fit_lin3 <- matrix(0, nrow = 323, ncol = 8)
fit_lin4 <- matrix(0, nrow = 323, ncol = 8)
fit_lin5 <- matrix(0, nrow = 323, ncol = 8)


phthalates_met2[, 327:334] <- scale(phthalates_met2[,327:334]) #scaling does not change the final results



for (i in 4:326)
  for (j in 327:334){
    
    # fit_lin3[i-3,j-330] <- summary(lm(phthalates_met2[,i] ~ phthalates_met2[, j] + sg + isage + prebmi, data = phthalates_met2))$coefficients[17] #p value
    # 
    # fit_lin4[i-3,j-330] <- summary(lm(phthalates_met2[,i] ~ phthalates_met2[, j] + sg + isage  + prebmi, data = phthalates_met2))$coefficients[2] #beta
    
    fit_lin3[i-3,j-326] <- summary(lm(phthalates_met2[,i] ~ phthalates_met2[, j] +isage + prebmi, data = phthalates_met2))$coefficients[2,4] #p value
    
    fit_lin4[i-3,j-326] <- summary(lm(phthalates_met2[,i] ~ phthalates_met2[, j] + isage  + prebmi, data = phthalates_met2))$coefficients[2,1] #beta
    
    fit_lin5[i-3,j-326] <- summary(lm(phthalates_met2[,i] ~ phthalates_met2[, j] + isage  + prebmi, data = phthalates_met2))$coefficients[2,2] #beta
    
  }


fit_lin3 <- as.data.frame(fit_lin3)
fit_lin4 <- as.data.frame(fit_lin4)
fit_lin5<- as.data.frame(fit_lin5)

# fit_lin3 <- as.data.frame(fit_lin3)

colnames(fit_lin3) <- names(phthalates_met2[,327:334])
colnames(fit_lin4) <- names(phthalates_met2[,327:334])
colnames(fit_lin5) <- names(phthalates_met2[,327:334])


fit_lin3 <- cbind.data.frame(fit_lin3, samp = names(phthalates_met2[,4:326]))
fit_lin4 <- cbind.data.frame(fit_lin4, samp = names(phthalates_met2[,4:326]))
fit_lin5 <- cbind.data.frame(fit_lin5, samp = names(phthalates_met2[,4:326]))


# Now, I will derive FDR adjusted p-values again. Using Q value approach would result in same conclusions regarding the significant metabolites.

# ```{r checking for significant metabolites after adding covariates}

#Define FDR threshold value

pthres = 0.050
pthres = 0.100

fit_lin3_adjust <- cbind.data.frame(samp = fit_lin3[,9], 
                                    apply(fit_lin3[,1:8], 2, function(x) p.adjust(x, method = "BH")))

fit_lin3_adjust_long <- fit_lin3_adjust %>% pivot_longer(!samp, names_to = "Phthalates", values_to = "adj_p") %>% 
  mutate(pval_ind = ifelse(adj_p > pthres, 0, 1))


fit_lin4_adjust_long <- fit_lin4 %>% pivot_longer(!samp, names_to = "Phthalates", values_to = "betas") 
fit_lin5_adjust_long <- fit_lin5 %>% pivot_longer(!samp, names_to = "Phthalates", values_to = "se") 


fit_lin34_adjust_long <- cbind.data.frame(fit_lin3_adjust_long, betas = fit_lin4_adjust_long$betas, se = fit_lin5_adjust_long$se)

fit_lin34_adjust_long <- fit_lin34_adjust_long %>% mutate(beta_sign = ifelse(betas < 0, "Negative", "Positive"))


########################################################################################################################################################################

# Now draw a Manhattan plot of p-values


group.colors <- c(Negative = "skyblue4", Positive = "tomato1")


fit_lin34_adjust_long <- fit_lin34_adjust_long %>% 
  mutate(adj_p = round(adj_p, 3)) %>%
  mutate(`Beta coefficient` = beta_sign, Label = ifelse(adj_p <= pthres, "1", "")) 

fit_lin34_adjust_long <- fit_lin34_adjust_long %>% mutate(Phthalates = ifelse(Phthalates == "MBZP", "MBzP", Phthalates))

fit_lin34_adjust_long2 <- fit_lin34_adjust_long %>% filter(adj_p <= pthres)
fit_lin34_adjust_long2$`Beta coefficient` <- factor(fit_lin34_adjust_long2$`Beta coefficient`, levels = c("Negative", "Positive"))


man_plot2 <- ggplot(fit_lin34_adjust_long, aes(x = samp, y = -log10(adj_p), color = `Beta coefficient`)) +    
  geom_point(shape = 15) +     
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),    
        panel.background = element_rect(fill = 'white'),    
        strip.background = element_rect(fill = 'white'), 
        strip.placement = "outside", 
        axis.line = element_line(),    
        legend.position = "bottom") +    
  scale_y_continuous(expression(paste(-log10('p-value'))), expand = c(0, 0)) + 
  scale_color_manual(values = c("Positive" = "skyblue4", "Negative" = "red")) +# Modified line
  facet_wrap(~Phthalates, strip.position = "bottom", nrow = 1)

man_plot2 <- man_plot2  + theme(panel.spacing.x = unit(0.1, "lines")) +    
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 2) +    
  geom_hline(yintercept = -log10(0.1), col = "blue", linetype = 2) +   

  scale_fill_manual(values = group.colors)

man_plot2


########################################################################################################################################################################

# Now, I will tabulate and plot only the significant metabolites


all_sig2 <- fit_lin34_adjust_long %>% filter(adj_p <= pthres) %>% arrange(Phthalates, adj_p) %>% 
  select(samp, Phthalates, adj_p, betas, se, beta_sign) %>% 
  rename("Metabolites" = "samp")

# all_sig2 %>% group_by(Phthalates) %>% dplyr::summarize(Nsig = n())

ggplot(all_sig2, aes(x = Phthalates, group = Phthalates, fill = Phthalates)) + geom_bar() + 
  theme(legend.position = "none", axis.title.x  = element_blank()) + 
  ggtitle("Number of significant metabolic markers per phthalate metabolite")  + 
  facet_wrap(~beta_sign)


########################################################################################################################################################################
# Plot significant metabolites using a heat plot. Here, I am only plotting metabolites with 
# adjusted p-values \< 0.1 to make it easier to visuaulize the result.


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 50)

all_sig2_plot <- all_sig2 %>% filter(adj_p <= pthres) %>% mutate(`Beta coefficient` = betas) 

# Create the heatmap
all_sig2_plot <- all_sig2_plot %>%  
  mutate(Metabolites2 = case_when(Metabolites == "AMINO FATTY ACIDS (C8H17NO2)" ~ "2-Aminooctanoic acid",
                                  Metabolites == "DICARBOXYLIC ACIDS (C16H30O4)" ~ "Hexadecanedioic acid",
                                  Metabolites == "DICARBOXYLIC ACIDS (C12H22O4)" ~ "Dodecanedioic acid",
                                  Metabolites == "DICARBOXYLIC ACIDS (C18H32O4)" ~ "Octadecenedioic acid",
                                  Metabolites == "FA(18:1(Ke)).1" ~ "FA(18:1(Ke))"))

all_sig2_plot <- all_sig2_plot %>% mutate(Metabolites = ifelse(is.na(Metabolites2), Metabolites, Metabolites2)) %>%
  mutate(Metabolites = toupper(Metabolites))



#Check # of significantly associated phthalates per metabolite


sum_all_sig2 <- all_sig2 %>% group_by(Metabolites) %>% summarise(NumPhthalates = n()) %>% 
  arrange(desc(NumPhthalates))

########################################################################################################################################################################

#Add classes of metabolites


###Add classes of significant metabolites

all_sig2_plot$Metabolites <- trimws(all_sig2_plot$Metabolites, which = c("both"))
all_sig2_plot$Metabolites <- toupper(all_sig2_plot$Metabolites)
all_sig2_plot <- all_sig2_plot %>% select(-Metabolites2)


##Change names of some metabolites in Classes file too

classes <- classes %>% mutate(Metabolites2 = case_when(Metabolites == "AMINO FATTY ACIDS (C8H17NO2)" ~ "2-Aminooctanoic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C16H30O4)" ~ "Hexadecanedioic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C12H22O4)" ~ "Dodecanedioic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C18H32O4)" ~ "Octadecenedioic acid",
                                                       Metabolites == "FA(18:1(Ke)).1" ~ "FA(18:1(Ke))"))
classes <- classes %>% mutate(Metabolites = ifelse(is.na(Metabolites2), Metabolites, Metabolites2)) %>%
  mutate(Metabolites = toupper(Metabolites))

all_sig2_class <- all_sig2_plot %>% left_join(classes, by = c("Metabolites"))

all_sig2_class2 <- all_sig2_class %>% filter(adj_p <= pthres)

all_sig2_class3 <- all_sig2_class2 %>% select(Metabolites, Phthalates, beta_sign, Class_new2)
all_sig2_class3_wide <- pivot_wider(all_sig2_class3, names_from = Phthalates, values_from = beta_sign)


# Now create a plot showing metabolites classes vs. phthalates


all_sig2_class2 <- all_sig2_class2 %>% mutate(`Beta coefficient2` = `Beta coefficient`, `Beta coefficient` = beta_sign)
all_sig2_class2$`Beta coefficient` <- factor(all_sig2_class2$`Beta coefficient`, levels = c("Negative", "Positive"))


classplot1 <- ggplot(all_sig2_class2, aes(y = Class_new2,  fill = `Beta coefficient`)) + geom_bar() + theme_classic() +
  facet_grid(~Phthalates, scales = "free_y") + theme(legend.position = "none", 
                                                     strip.background = element_rect(fill = 'white', linetype = 0), 
                                                     strip.placement = "outside", axis.line = element_line(), 
                                                     panel.spacing.x = unit(0.5, "cm"), axis.ticks.y = element_blank())  +
  # ggtitle("Figure 1. Class distribution of significantly associated metabolic markers from MWAS") +
  xlab("Count")  + ylab("Metabolite classes")+
  scale_fill_manual(values = group.colors)

classplot1


#Also make the heatplot

all_sig2_class2 <- all_sig2_class2 %>%
  arrange(Class_new2)

ggplot(all_sig2_class2, aes(Phthalates, Metabolites, fill = `Beta coefficient2`)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "skyblue4", mid = "white", 
                       midpoint = 0, limit = c(-0.5,0.5), 
                       name = "Beta coefficient") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1), legend.position = "bottom") +
  labs(x = "Phthalates", y = "Plasma metabolic features")



#Also export significant file 

sigfile_to_export <- all_sig2_class2 %>%
  mutate(
    lower_ci = round(betas - 1.96 * se, 2),
    upper_ci = round(betas + 1.96 * se, 2),
    beta_ci = sprintf("%.2f (%.2f, %.2f)", betas, lower_ci, upper_ci), # Using sprintf to format the string
    adj_p = round(adj_p, 3)) %>%  # Round adjusted p-values to 3 decimal points
  select(Phthalates, Metabolites, Class_new2, beta_ci, adj_p) %>% arrange(Phthalates, Class_new2) 

library(flextable)
sigfile_to_export_flex <- flextable::flextable(sigfile_to_export)
sigfile_to_export_flex <- flextable::merge_v(sigfile_to_export_flex, j = ~Phthalates)
sigfile_to_export_flex <- flextable::merge_v(sigfile_to_export_flex, j = ~Class_new2)
flextable::save_as_docx(sigfile_to_export_flex, path = paste0(getwd(), "/significant_metabolites_by_phthalates_prenatal.docx"),
                        align = "left")

write.csv(sigfile_to_export, file = "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/manuscript_writing/significant_metabolites_by_phthalates_prenatal.csv")


#####################################################################################################################################################
#Look at original p-value - added on 1/28/2024

origpvalue_prenat <- cbind.data.frame(Metabolite = fit_lin3[,9], fit_lin3[,1:8]) %>%
   mutate(Metabolites2 = case_when(Metabolite == "AMINO FATTY ACIDS (C8H17NO2)" ~ "2-Aminooctanoic acid",
                                                         Metabolite == "DICARBOXYLIC ACIDS (C16H30O4)" ~ "Hexadecanedioic acid",
                                                         Metabolite == "DICARBOXYLIC ACIDS (C12H22O4)" ~ "Dodecanedioic acid",
                                                         Metabolite == "DICARBOXYLIC ACIDS (C18H32O4)" ~ "Octadecenedioic acid",
                                                         Metabolite == "FA(18:1(Ke)).1" ~ "FA(18:1(Ke))",
                                   TRUE ~ Metabolite))

origpvalue_prenat2 <- origpvalue_prenat %>% 
  rename("Metabolites" = "Metabolites2") %>% 
  pivot_longer(cols = DEHP:MCOP, names_to = "Phthalates", values_to = "Unadj_p") %>%
  mutate(Metabolites = toupper(Metabolites), Phthalates = toupper(Phthalates)) %>%  
  inner_join(sigfile_to_export, by = c("Metabolites", "Phthalates"))

sigfile_to_export %>% anti_join(origpvalue_prenat2, by = c("Metabolites", "Phthalates"))



# Define the desired order of Metabolites2
origpvalue_prenat2_mcpp <- origpvalue_prenat2 %>% filter(Phthalates == "MCPP")
origpvalue_prenat2_dehp <- origpvalue_prenat2 %>% filter(Phthalates == "DEHP")
origpvalue_prenat2_mep <- origpvalue_prenat2 %>% filter(Phthalates == "MEP")

mcpp_order <- c("L-GAMMA-GLUTAMYL-L-ISOLEUCINE", "L-GAMMA-GLUTAMYL-L-LEUCINE", "GLUTAMYLPHENYLALANINE", 
                "N-ACETYLLEUCINE", "N-ACETYL-L-PHENYLALANINE", "UROCANIC ACID", "1-AMINOCYCLOPROPANECARBOXYLIC ACID", 
                "L-PHENYLALANINE", "L-TRYPTOPHAN", "2-AMINOOCTANOIC ACID", "L-ISOLEUCINE", "L-METHIONINE", 
                "L-TYROSINE", "ISOLEUCYL-ISOLEUCINE", "4-HYDROXY-L-PROLINE", "L-LYSINE", "HYDROXYPROPIONIC ACID",
                "BRANCHED FATTY ACIDS (C23H46O2)", "ISOVALERYLCARNITINE", "CAR(5:1)", "CAR(10:3)", "CER(D40:1)",
                "CER(D42:2)", "CER(D41:2)", "DG(16:0/16:0/0:0)", "DODECANEDIOIC ACID", "METHYLMALONIC ACID", 
                "2-HYDROXY-3-METHYLBUTYRIC ACID", "3-METHYL-2-OXOVALERIC ACID", "LPC(18:2)", "FA(6:0)",
                "FA(10:0(OH))", "PIMELIC ACID", "VANILLYLMANDELIC ACID", "4-NITROPHENOL", "PC(32:0)", 
                "DI(2-ETHYLHEXYL)PHTHALATE", "ADENINE", "GUANINE", "1-METHYLADENOSINE", "2-ACETYLPYRROLIDINE", 
                "FA(5:0(OH))", "SM(D18:1/24:0)", "SM(D34:1)", "PROGESTERONE", "PREGNENOLONE")

origpvalue_prenat2_mcpp <- origpvalue_prenat2_mcpp %>%
  arrange(factor(Metabolites, levels = mcpp_order))


dehp_order <- c("N-ACETYLLEUCINE", "HEXADECANEDIOIC ACID", "OCTADECENEDIOIC ACID",
                   "FA(26:0)", "BEHENIC ACID", "FA(6:0(OH))", "FA(11:1)", "4-NITROPHENOL",
                   "DI(2-ETHYLHEXYL)PHTHALATE", "1,7-DIMETHYLGUANOSINE", "PREGNENOLONE",
                   "PROGESTERONE")

origpvalue_prenat2_dehp <- origpvalue_prenat2_dehp %>%
  arrange(factor(Metabolites, levels = dehp_order))



mep_order <- c("1-AMINOCYCLOPROPANECARBOXYLIC ACID", "N-ACETYLLEUCINE", "L-METHIONINE", "L-ISOLEUCINE",
               "N-ACETYL-L-PHENYLALANINE", "L-GAMMA-GLUTAMYL-L-LEUCINE", "PROTOCATECHUIC ACID", "CER(D40:1)", 
               "FERULIC ACID", "3-METHYL-2-OXOVALERIC ACID", "FA(12:0(KE))", "DI(2-ETHYLHEXYL)PHTHALATE", "GUANINE")


origpvalue_prenat2_mep <- origpvalue_prenat2_mep %>%
  arrange(factor(Metabolites, levels = mep_order))

origpvalue_prenat3 <- rbind.data.frame(origpvalue_prenat2_dehp, origpvalue_prenat2_mcpp, origpvalue_prenat2_mep) %>%
  mutate(Unadj_p = round(Unadj_p, 6))

write.csv(origpvalue_prenat3, "origpvalue_prenat.csv")

#########################################################################################################################################################


#Add class to a file with all metabolites
fit_lin34_adjust_long_chk <- fit_lin34_adjust_long %>% mutate(Metabolites = toupper(samp)) %>%
  inner_join(classes, by = c("Metabolites")) %>%
  mutate(
  lower_ci = round(betas - 1.96 * se, 2),
  upper_ci = round(betas + 1.96 * se, 2),
  beta_ci = sprintf("%.2f (%.2f, %.2f)", betas, lower_ci, upper_ci), # Using sprintf to format the string
  adj_p = round(adj_p, 3)) %>%  # Round adjusted p-values to 3 decimal points
  select(Phthalates, Metabolites, Class_new2, beta_ci, adj_p) %>% arrange(Phthalates, Class_new2) 
  

########################################################################################################################################################################
#Hierchical clustering


#Create correlation matrix

mets_forCor <- phthalates_met2[, 4:326]
cor_met <- rcorr(as.matrix(mets_forCor), type="spearman")
cor_met_sp <- cor_met$r
cor_met_sp2 = 1-(cor_met_sp)


#Correlation with data that excludes one participant

mets_forCor1 <- phthalates_met2_short[, 4:326]
cor_met1 <- rcorr(as.matrix(mets_forCor1), type="spearman")
cor_met_sp1 <- cor_met1$r
cor_met_sp2_1 = 1-(cor_met_sp1)




#Hierarchial clustering

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 50)
# heatmap.2(cor_met_sp, dendrogram="col",  density.info="none", trace="none",
#           main="Metabolite Correlations", margins =c(12,9), col=my_palette)


#Function to create a cluster

getclus <- function(x, data){
  hc.rows <- hclust((data))
  ct <- cutree(hc.rows, h=x)  # cut the dendrogram into clusters
  plot(hc.rows,hang = -2, cex = 0.4)
  dend <- rect.hclust(hc.rows, h = x) # visualize the groups
  sum_ct <- table(ct)
  return(list(dend, ct, sum_ct))
}

#Mean correlation per cluster 

get_mean <- function(clus){
  mean_all <- NULL
  for(i in 1: length(clus[[1]])){
    # gp <- names(clus[[1]][[i]])
    gp <- names(which(clus[[2]]==i))
    gp_cor <-cor_met_sp[which(rownames(cor_met_sp) %in% gp ==TRUE),which(colnames(cor_met_sp) %in% gp ==TRUE)]
    # gp_cor <- abs(cor_met_sp[which(rownames(cor_met_sp) %in% gp ==TRUE),which(colnames(cor_met_sp) %in% gp ==TRUE)])
    cor_mean <- mean(gp_cor[lower.tri(gp_cor)],na.rm=TRUE)
    row=data.frame(group=i, num=length(gp),cor_mean=round(cor_mean,3))
    mean_all <- rbind(mean_all, row)
  }
  return(mean_all)
}

#See phthalates_metabolites_62023 for details. However, we will use h = 1.1

##############################################################################################################################################################
#where to cut the tree


seqval <- seq(0.5,2, 0.05)
matval <- matrix(nrow = length(seqval), ncol = 4)
# corchk_temp <- vector(nro)


for (i in seqval){
  optclus_temp <- getclus(x=i, data = as.dist(cor_met_sp2))
  optclus_mean_temp <- get_mean(optclus_temp)
  pos <- match(i, seqval)
  matval[pos, 1] <- i
  matval[pos, 2] <- length(which(optclus_mean_temp$cor_mean < .3)) #number of clusters with avg r mean lt 0.15
  matval[pos, 3] <- length(which(optclus_mean_temp$num < 5)) #number of clusters with <5 metabolic markers
  matval[pos, 4] <- length(optclus_temp[[3]]) #number of clusters
  
}



matval <- as.data.frame(matval)
names(matval) <- c("len", "n_rlt2", "n_nlth5", "n_clus")

matval_n <- matval %>% select(len, n_nlth5, n_clus) %>% mutate(Criteria = "ClusSize >= 5")
names(matval_n) <- c("length", "n_clus_bad", "n_cluster", "Criteria") #n_clus_bad = clusters not meeting the requirement
matval_r <- matval %>% select(len, n_rlt2, n_clus) %>% mutate(Criteria = "Pairwise_r >= 0.3")
names(matval_r) <- c("length", "n_clus_bad", "n_cluster", "Criteria") 

matval_long <- rbind.data.frame(matval_n, matval_r)


ggplot(matval_long, aes(x = length, y = n_clus_bad, col = Criteria)) + geom_point()
# ggplot(matval, aes(x = len, y = n_nlth5)) + geom_point()

#Also use elbow method 


plot(matval$len, matval$n_nlth5)
abline(v = 1.1)

##############################################################################################################################################################


#Get clusters with at least 5 metabolites and r < 0.3
set.seed(1991)
cluslen <- 1.1
optclus <- getclus(cluslen, data = as.dist(cor_met_sp2))
optclus_mean <- get_mean(optclus)


##############################################################################################################################################################

#Better way to visualize cluster 

hc.rows <- hclust((as.dist(cor_met_sp2)))
plot(hc.rows, hang = -2, main = "", cex = 0.3, sub = "", xlab = "") 
dend <- rect.hclust(hc.rows, h = 1.1) # visualize the groups
title(main = "Hierarchical clustering for maternal plasma metabolic features", cex.main = 1)

clust_prenat <- cutree(hc.rows, h = 1.1)
clust_prenat
table(clust_prenat)

#cluster 1, 2, 5, 8, 10, 18, 19, 23 do not meet criteria 


##############################################################################################################################################################

corchk <- which(optclus_mean$cor_mean <0.3)
sizechk <- which(optclus_mean$num <5)

metabolites_t3 <- phthalates_met2[,c(2,4:326)]
t_metabolites_t3 <- data.table::transpose(metabolites_t3)
rownames(t_metabolites_t3) <- colnames(metabolites_t3)
colnames(t_metabolites_t3) <- t_metabolites_t3[1,]
t_metabolites_t4 <- t_metabolites_t3[-1,]
t_metabolites_t4 <- t_metabolites_t4 %>% mutate(metabolite = rownames(t_metabolites_t4), Cluslabel = optclus[[2]]) %>% 
  select(metabolite, Cluslabel, everything())



#After removing one participant
# 
# optclus1 <- getclus(cluslen, data = as.dist(cor_met_sp2_1))
# optclus_mean1 <- get_mean(optclus1)
# 
# corchk1 <- which(optclus_mean1$cor_mean <0.3)
# sizechk1 <- which(optclus_mean1$num <5)
# 
# metabolites_t3_1 <- phthalates_met2_short[,c(2,4:326)]
# t_metabolites_t3_1 <- data.table::transpose(metabolites_t3_1)
# rownames(t_metabolites_t3_1) <- colnames(metabolites_t3_1)
# colnames(t_metabolites_t3_1) <- t_metabolites_t3_1[1,]
# t_metabolites_t4_1 <- t_metabolites_t3_1[-1,]
# t_metabolites_t4_1 <- t_metabolites_t4_1 %>% mutate(metabolite = rownames(t_metabolites_t4_1), Cluslabel = optclus1[[2]]) %>% 
#   select(metabolite, Cluslabel, everything())
# 
# 
# cluschk_all <- t_metabolites_t4 %>% select(metabolite, Cluslabel)
# cluschk_all1 <- t_metabolites_t4_1 %>% select(metabolite, Cluslabel)
# 
# 
# cluschk_chk <- cluschk_all %>% inner_join(cluschk_all1, by = "metabolite")
# 
# cluschk_chkb <- cluschk_chk %>% group_by(Cluslabel.x, Cluslabel.y) %>% tally() %>%
#   filter(Cluslabel.x %in% c(3,4,6,7,9,12,14,15,16,17,21)) #Clusters that fulfill the criteria (see below)

########################################################################################################################################################################


###Multivariate regression using the cluster memberships

# Following Margaret's code

#list of clusterid
clusterid <- t_metabolites_t4$Cluslabel

#number of metabolites in a cluster
clustersum <- aggregate(rep(1,nrow(t_metabolites_t4)) ~ clusterid, FUN = sum)

#Choose the clusters to use 

clustertouse <- optclus_mean %>% filter(num >=5, cor_mean >= 0.3)
sum(clustertouse$num)



clus_chk1 <- t_metabolites_t4 %>% select(metabolite, Cluslabel) %>% rename("group" = "Cluslabel")

clus_chk1 <- clus_chk1 %>% inner_join(optclus_mean, by = c("group")) %>% arrange(as.numeric(group))

#Changing the format to store export 

clus_chk1_toexport <- clus_chk1 %>%
  group_by(group) %>%
  summarise(
    Group_Correlation = first(cor_mean), # Group correlation for each cluster
    NumberofMetabolites = first(num),
    Metabolites = paste(metabolite, collapse = ", ")
  ) %>%
  arrange(group) %>%
  print(n = Inf) # To print all rows

write.csv(clus_chk1_toexport, "prenatal_metabolites_cluster_assignment.csv")


ttmetabolites_t4 <- as.data.frame(t(t_metabolites_t4))
ttmetabolites_t5 <- ttmetabolites_t4[3:nrow(ttmetabolites_t4),]


names(phthalates_met2) <- toupper(names(phthalates_met2))


row_all2 <- NULL

for (j in 327:334){
  z = data.frame(Phth=phthalates_met2[,j] , AGE = phthalates_met2$ISAGE, BMI=phthalates_met2$PREBMI)
  zz = data.matrix(z)

  for(i in clustertouse[,1]){
    ik = i
    x = ttmetabolites_t5[, clusterid == ik]
    xx = x
    xxd = data.frame(xx)
    pp = ncol(xxd)
    names(xxd) = c(paste0("M", 1:pp))

    data1 = data.frame(Phth = zz[,1], AGE = scale(zz[,2]), BMI=scale(zz[,3]), xxd)
    yvec <- as.numeric(as.matrix(cbind(data1[,4:ncol(data1)])))
    yvec <- matrix(yvec, ncol = ncol(data1)-3)
    xvec <- as.matrix(cbind(data1[,1:3]))
    my.model <- lm(yvec ~ xvec) 
    
    xvec2 <- as.matrix(cbind(data1[,2:3]))
    mlm2 <- lm(yvec ~ xvec2)
    pval <- as.numeric(anova(my.model, mlm2)$'Pr(>F)'[2])

    row <- cbind(Metab_Cluster=ik, Phthalate = colnames(phthalates_met2[j]),pval = pval, cluster_size = pp)

    row_all2 <- rbind(row_all2, row)
    # row_all2 <- as.data.frame(row_all2) %>% filter(pval < 0.05)
  }
}



# Now look at the phthalate metabolites that were significant after adjusting for multiple testing.

row_all2_bhcor_chk <- data.frame(row_all2) %>% 
  mutate(pval = as.numeric(as.character(pval))) %>% 
  group_by(Phthalate) %>% 
  mutate(BH = p.adjust(pval, method = "BH"), BH_round = round(BH, 3))

row_all2_bhcor <- row_all2_bhcor_chk %>% 
 filter(BH <= pthres) #This is done to make sure that the cluster with p-value of 0.05 is not excluded. 
  # dplyr::mutate(q = qvalue(pval)$qvalues)  # Not enough tests per phthalate for qvalue

row_all2_bhcor_chk %>% arrange(as.numeric(Metab_Cluster))

#export the results file 
row_all2_bhcor_chk_wide <- row_all2_bhcor_chk %>% select(Metab_Cluster, Phthalate, BH_round) %>% 
  pivot_wider(names_from = Phthalate, values_from = BH_round)

write.csv(row_all2_bhcor_chk_wide, "hierarchical_cluster_results_prenatal.csv")

########################################################################################################################################################################


# Now, add class label to significantly associated metabolites


row_all2_bhcor_sum <- row_all2_bhcor %>% group_by(Metab_Cluster) %>% tally() 

met_clus <- t_metabolites_t4 %>% select(metabolite, Cluslabel) %>% dplyr::rename("Metab_Cluster" = Cluslabel, "Metabolites" = metabolite) %>%
  mutate(Metab_Cluster = as.character(Metab_Cluster))

met_clus <- met_clus %>% mutate(Metabolites = toupper(Metabolites))


met_clus2 <- met_clus %>% inner_join(row_all2_bhcor_sum, by = "Metab_Cluster")

met_clus2 <- met_clus2 %>% mutate(Metabolites2 = case_when(Metabolites == "AMINO FATTY ACIDS (C8H17NO2)" ~ "2-Aminooctanoic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C16H30O4)" ~ "Hexadecanedioic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C12H22O4)" ~ "Dodecanedioic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C18H32O4)" ~ "Octadecenedioic acid",
                                                       Metabolites == "FA(18:1(Ke)).1" ~ "FA(18:1(Ke))"))

met_clus2 <- met_clus2 %>% mutate(Metabolites = ifelse(is.na(Metabolites2), Metabolites, Metabolites2)) %>%
  mutate(Metabolites = toupper(Metabolites))

met_clus2 <- met_clus2 %>% left_join(classes, by = c("Metabolites"))

# write.csv(met_clus2, "significant_cluster_metabolites_list.csv")


# Check the correlation structure for each clusters


# Split the data frame into a list of data frames
list_of_dfs <- split(met_clus2, met_clus2$Metab_Cluster)


# Initialize a list to store the correlation matrices
cor_matrices <- list()
dev.off()

# Loop over each data frame in the list
# for(i in seq_along(list_of_dfs)) {
#   # Get the metabolites for the current cluster
#   metabolites <- list_of_dfs[[i]]$Metabolites
#   
#   # Subset phthalates_met2 to include only the columns corresponding to the metabolites
#   subsetdat <- phthalates_met2[, metabolites]
#   
#   # Create a correlation matrix for the subset and store it in the list
#   cor_matrices[[i]] <- cor(subsetdat, use = "complete.obs", method = "spearman")
#   
#   # Plot the correlation matrix
#   corrplot(cor_matrices[[i]], method = "color")
# }

#########################################################################################################################################################################

# Let us look at what cluster has the highest number of signficant associations 

row_all2_bhcor_sum %>% arrange(desc(n))


#Also look at what class is most represented
met_clus2 %>% group_by(Metab_Cluster, Class_new2) %>% tally() %>% arrange(Metab_Cluster, desc(n))


# phthalates_met2[,which(clusterid == 10)]

met_clus2 <- met_clus2 %>% mutate(Metab_Cluster2 = paste("Cluster", Metab_Cluster))
met_clus2$Metab_Cluster3 <- factor(met_clus2$Metab_Cluster2, levels = c("Cluster 4", "Cluster 6", "Cluster 7", "Cluster 12",
                                                                      "Cluster 15", "Cluster 16", "Cluster 17", "Cluster 20", "Cluster 21"))


ggplot(data = met_clus2, aes(y = Class_new2, fill = Class_new2)) +
  geom_bar() +
  facet_wrap(~as.factor(Metab_Cluster3), nrow = 2) +
  theme_classic() +
  ylab("Classes") +
  theme(
    legend.position = "right", legend.title = element_blank(),  legend.text=element_text(size=6),
    axis.text.y = element_text(size = 7),
    strip.background = element_blank(), # Remove the box around the facet name
    strip.text = element_text(vjust = 0), # Adjust text positioning if needed
    panel.border = element_rect(colour = "black", fill = NA, size = 1) # Add axis line to each facet
  )

  
  




########################################################################################################################################################################


# Now combine results from cluster and single pollutant analysis 

#Get rid of duplicate metabolites

met_clus2 <- met_clus2 %>% filter(substr(Metabolites, nchar(Metabolites), nchar(Metabolites)) != 1)

met_clus2b <- met_clus2 %>% select(Metabolites, Metab_Cluster, n, Class_new2)
all_sig2_class2b <- all_sig2_class2 %>% select(Metabolites, Phthalates, adj_p, beta_sign)
fit_lin34_adjust_longb <- fit_lin34_adjust_long %>% select(samp, Phthalates, adj_p, beta_sign, pval_ind)
fit_lin34_adjust_longb <- fit_lin34_adjust_longb %>% rename("Metabolites" = "samp") %>% 
  mutate(sig_sign = paste(pval_ind, beta_sign, sep = "_")) %>% select(-adj_p, -beta_sign, -pval_ind)
fit_lin34_adjust_longb_wide <- pivot_wider(fit_lin34_adjust_longb, names_from = Phthalates, values_from = sig_sign) %>%
  mutate(Metabolites = toupper(Metabolites))



fit_lin34_adjust_longb_wide <- fit_lin34_adjust_longb_wide %>% mutate(Metabolites2 = case_when(Metabolites == "AMINO FATTY ACIDS (C8H17NO2)" ~ "2-Aminooctanoic acid",
                                                          Metabolites == "DICARBOXYLIC ACIDS (C16H30O4)" ~ "Hexadecanedioic acid",
                                                          Metabolites == "DICARBOXYLIC ACIDS (C12H22O4)" ~ "Dodecanedioic acid",
                                                          Metabolites == "DICARBOXYLIC ACIDS (C18H32O4)" ~ "Octadecenedioic acid",
                                                          Metabolites == "FA(18:1(Ke)).1" ~ "FA(18:1(Ke))"))

fit_lin34_adjust_longb_wide <- fit_lin34_adjust_longb_wide %>% mutate(Metabolites = ifelse(is.na(Metabolites2), Metabolites, Metabolites2)) %>%
  mutate(Metabolites = toupper(Metabolites))

#Get rid of duplicate metabolites

fit_lin34_adjust_longb_wide <- fit_lin34_adjust_longb_wide %>% filter(substr(Metabolites, nchar(Metabolites), nchar(Metabolites)) != 1)


#Create subsets of metabolites for each cluster

for (i in row_all2_bhcor_sum$Metab_Cluster){
  clusname <- paste0("met_clus_ind", i)
  met_clusi <- met_clus2b %>% filter(Metab_Cluster == i) %>% 
    left_join(fit_lin34_adjust_longb_wide, by = "Metabolites") %>% arrange(Class_new2)
  assign(clusname, met_clusi)
  
}

# write.csv(met_clus_ind10, "cluster10.csv")
# write.csv(met_clus_ind15, "cluster15.csv")
# write.csv(met_clus_ind16, "cluster16.csv")
# write.csv(met_clus_ind19, "cluster19.csv")
# write.csv(met_clus_ind5, "cluster5.csv")
# write.csv(met_clus_ind6, "cluster6.csv")


#Combine three clusters
# met_clus_ind_top3 <- rbind.data.frame(met_clus_ind10, met_clus_ind1, met_clus_ind13)


# Now, let us plot metabolites from each cluster for Figure 2

sigmet_clus4 <- c("MCPP")
sigmet_clus6 <- c("MBP", "MCPP")
sigmet_clus7 <- c("MCPP")
sigmet_clus12 <- c("MCPP")
sigmet_clus15 <-c("MBP", "MCPP")
sigmet_clus16 <- c("MCPP")
sigmet_clus17 <- c("MCPP", "MCOP")
sigmet_clus20 <-c("MCPP")
sigmet_clus21 <- c("MCPP")



datprep_tile <- function(data, vars, title){
  
  #Create the data
  dat <- data %>% select(Metabolites, Class_new2, all_of(vars))
  dat_long <- dat %>% pivot_longer(cols = vars, names_to = "Phthalates", values_to = "Beta coefficient",
                                   values_ptypes  = list(`Beta coefficient` = 'character'))
  dat_long2 <- dat_long %>% mutate(`Beta coefficient2` = substr(`Beta coefficient`, 3, 10),
                                   signif = substr(`Beta coefficient`, 1, 1),
                                   signif2 = ifelse(signif == 1, "*", ""))
  
  dat_long2$`Beta coefficient2` = factor(dat_long2$`Beta coefficient2`, levels = c("Negative", "Positive"))
  dat_long2 <- dat_long2 %>% arrange(Class_new2)
  
  #Plot 
  
  plot1 <- ggplot(dat_long2, aes(x = Phthalates, y = Metabolites, fill = `Beta coefficient2`, label = signif2)) + 
    geom_tile() + geom_text(aes(label = signif2), col = "white") + theme_classic() + 
    scale_fill_manual(values = group.colors) +
    ggtitle(title) + theme(axis.title  = element_blank(), legend.position = "none", axis.ticks = element_blank(),
                           axis.text.y = element_text(size = 6))
  
  return(plot1)
  
}

clus4_plot <- datprep_tile(data = met_clus_ind4, vars = sigmet_clus4, title = "Cluster 4")
clus4_plot

clus6_plot <- datprep_tile(data = met_clus_ind6, vars = sigmet_clus6, title = "Cluster 6") 
clus6_plot

clus7_plot <- datprep_tile(data = met_clus_ind7, vars = sigmet_clus7, title = "Cluster 7") 
clus7_plot

clus12_plot <- datprep_tile(data = met_clus_ind12, vars = sigmet_clus12, title = "Cluster 12")
clus12_plot

clus15_plot <- datprep_tile(data = met_clus_ind15, vars = sigmet_clus15, title = "Cluster 15") 
clus15_plot

clus16_plot <- datprep_tile(data = met_clus_ind16, vars = sigmet_clus16, title = "Cluster 16") 
clus16_plot

clus17_plot <- datprep_tile(data = met_clus_ind17, vars = sigmet_clus17, title = "Cluster 17")
clus17_plot

clus20_plot <- datprep_tile(data = met_clus_ind20, vars = sigmet_clus20, title = "Cluster 20")
clus20_plot

clus21_plot <- datprep_tile(data = met_clus_ind21, vars = sigmet_clus21, title = "Cluster 21") 
clus21_plot



# text <- "Table 2. List of metabolic markers in significantly associated clusters. * represents metabolic markers that were significantly associated in MWAS."
# text.p <- ggpubr::ggparagraph(text = text, face = "italic", size = 11, color = "black")
# 
ggpubr::ggarrange(clus4_plot,  clus6_plot, clus7_plot, clus12_plot,
                  align = "v")


ggpubr::ggarrange(clus15_plot, clus16_plot, clus17_plot, clus20_plot, clus21_plot,
                  align = "v")

library(patchwork)

plotarrange <- (clus4_plot + clus6_plot + clus7_plot + clus12_plot + clus15_plot + 
                  clus16_plot + clus17_plot + clus20_plot + clus21_plot)
plotarrange

plotarrange2 <- (clus17_plot + clus20_plot + clus21_plot)
plotarrange2


plotarrange2 <- (clus6_plot + clus7_plot + clus12_plot)
plotarrange2

#######################################################################################################################
#Cluster stability check

#Another way to look at cluster data

clusfit <- hclust(as.dist(cor_met_sp2))
clusfit_groups <- cutree(clusfit, h=1.1) 

#bootstrap

library(fpc)
clusters <- 23
clus.boot <- clusterboot(as.dist(cor_met_sp2), 
                         B=1000, # Number of bootstrap resamples
                         clustermethod=hclustCBI, # for hierarchical clustering 
                         method="ward.D", # use what we used in "hclust"
                         k = clusters, 
                         count=FALSE) # Show progress on screen?
# clus.boot

AvgJaccard <- clus.boot$bootmean
Instability <- clus.boot$bootbrd/1000
Clusters <- c(1:clusters)
Eval <- cbind(Clusters, AvgJaccard, Instability)
Eval

#We need to check clusters with low jaccard index to see what is going on. Perhaps those are the ones that
#were not selected

# 1. Hierarchical clustering
hclust_obj <- hclust(as.dist(cor_met_sp2), method="ward.D")

# 2. Get cluster assignments
clus_assignment <- cutree(hclust_obj, k=clusters)

# 3. Identify clusters with a low Jaccard index
threshold <- 0.5  # Adjust this value as per your definition of 'low'
Eval <- as.data.frame(Eval)
low_jaccard_clusters <- Eval$Clusters[Eval$AvgJaccard < 0.5]

# 4. Extract members of those clusters
low_jaccard_data <- list()

for (i in low_jaccard_clusters) {
  cluster_members <- cor_met_sp2[clus_assignment == i, ]
  low_jaccard_data[[i]] <- cluster_members
}

#Low jaccard clusters: 1  2  3  4  7 10 12 13 16 19 22

#However, only Cluster 4, 7, 12, 16 were selected 
# low_jaccard_data now contains the data points for clusters with a low Jaccard index



#################################################################################################################################
#The next analytical approach is to look at correlation network. I will export data to do that

# phthalates_met2_toexport <- phthalates_met2 %>% rename("2-Aminooctanoic acid" = "AMINO FATTY ACIDS (C8H17NO2)", 
#                                                        "Hexadecanedioic acid" = "DICARBOXYLIC ACIDS (C16H30O4)",
#                                                        "Dodecanedioic acid" = "DICARBOXYLIC ACIDS (C12H22O4)",
#                                                        "Octadecenedioic acid" = "DICARBOXYLIC ACIDS (C18H32O4)")
# write.csv(phthalates_met2_toexport[, 4:326], "prenatal_metfile.csv")

# metlist <- data.frame(Metabolites = toupper(names(phthalates_met2_toexport[, 4:326])))
# 
# classes_toexport <- classes %>% rename("class" = "Class_new2") %>% select(-Metabolites2) %>% mutate(Metabolites = toupper(Metabolites))
# classes_toexport <- classes_toexport %>% inner_join(metlist, by = c("Metabolites"))
write.csv(classes_toexport, "classes_prenatal_metfile.csv")



#Classes of metabolic features


#Distribution of classes

class_sum_prenat <- classes %>% group_by(Class_new2) %>% tally()



