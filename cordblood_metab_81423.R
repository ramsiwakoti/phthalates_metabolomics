#Date: 8/16/2023

#In this file, I will run main analysis for cordblood metabolomics data for manuscript production. Most of the code 
#will come from phthalates_metabolites_cordblood_62273.


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
conflicted::conflicts_prefer(dplyr::first)

########################################################################################################################################################################

#Obtain the analytic file created in phthalates_metabolites_62023. 

#This file has phthalates data, metabolites data, as well as covariates. Phthalates have been log-transformed
#whereas metabolites data have been auto-scaled. 

load("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/cord_blood_data/phthalates_met2_cordblood.RData")
classes <- read.csv("C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/cord_blood_data/list_of_cordblood_metabolites2_csv.csv") %>% distinct()
classes <- classes %>% mutate(class = trimws(class), Class_new2 = class) %>% rename("Metabolites" = "name")

phthalates_all <- read.csv("C:\\Users\\rsiwa\\Dropbox (University of Michigan)\\Ram_shared\\data\\2_v54_006_Phthalates_no_dup_03092021_20210625153132.csv") %>%
  select(StudyID, VISITID, CONC, BLOD, LOD, ANALYTE)
phthalates_all$StudyID<- as.character(phthalates_all$StudyID)

names(phthalates_all) <- c("Study ID", "VISITID", "CONC", "BLOD", "LOD", "ANALYTE")

phthalates_cordblood <- phthalates_met2 %>% select(`Study ID`) %>% inner_join(phthalates_all, by = "Study ID") 


########################################################################################################################################################################

#First check detection rate and missing rate for each phthalate metabolite


# Convert BLOD to a logical column (1 = TRUE, 0 = FALSE)
phthalates_cordblood$BLOD <- phthalates_cordblood$BLOD == 1

# Calculate detection rate
detection_rate <- phthalates_cordblood %>%
  filter(!is.na(CONC)) %>%  # Remove rows where CONC is NA
  group_by(ANALYTE, VISITID) %>%  # Group by analyte and visit ID
  summarise(
    Total = n(),  # Total number of samples
    Detected = sum(!BLOD),  # Number of samples above LOD
    DetectionRate = Detected / Total  # Calculate detection rate
  )

# View the results
print(detection_rate)



# Assuming your original data is in a dataframe named phthalates_cordblood
# Pivot data to wide format
phthalates_cordblood_wide <- phthalates_cordblood %>% select(-BLOD, -LOD) %>%
  pivot_wider(names_from = ANALYTE, values_from = CONC)

# Calculate missing values per analyte per visit for each participant
# missing_values <- phthalates_cordblood_wide %>%
#   group_by(Study_ID, VISITID) %>%
#   summarise(across(everything(), ~sum(is.na(.))))
# 
# # View the results
# print(missing_values)



###########################################################################################################################################################################

#Correlations between different phthalates

phthalates <- phthalates_met2[,338:347]
pattern <- "(log_SG_)|(\\_GM)"
phthalates_names <- toupper(gsub(pattern = pattern, replacement = "", names(phthalates)))
phthalates_names

names(phthalates) <- phthalates_names

corph <- corrplot(corr = cor(phthalates, method = "spearman", use = "complete"), method = "number", type = "lower")




##########################################################################################################################################
#Run MWAS 


metstart <- 4
metend <- 276
# phstart <- 333
# phend <- 344
phstart <- 338
phend <- 347


fit_lin3 <- matrix(0, nrow = 273, ncol = 10)
fit_lin4 <- matrix(0, nrow = 273, ncol = 10)
fit_lin5 <- matrix(0, nrow = 273, ncol = 10)

phthalates_met2[, 338:347] <- scale(phthalates_met2[,338:347]) #scaling does not change the final results - changes betas but not p-values



for (i in metstart:metend)
  for (j in phstart:phend){
    
    fit_lin3[i-metstart+1,j-phstart+1] <- summary(lm(phthalates_met2[,i] ~ phthalates_met2[, j] + ISAGE + as.factor(PPSEX) +  prebmi, data = phthalates_met2))$coefficients[17] #p-values
    
    fit_lin4[i-metstart+1,j-phstart+1] <- summary(lm(phthalates_met2[,i] ~ phthalates_met2[, j] + ISAGE + as.factor(PPSEX) +  prebmi, data = phthalates_met2))$coefficients[2] #betas
    
    fit_lin5[i-metstart+1,j-phstart+1] <- summary(lm(phthalates_met2[,i] ~ phthalates_met2[, j] + ISAGE + as.factor(PPSEX) +  prebmi, data = phthalates_met2))$coefficients[7] #se
    

  }


fit_lin3 <- as.data.frame(fit_lin3)
fit_lin4 <- as.data.frame(fit_lin4)
fit_lin5<- as.data.frame(fit_lin5)

# fit_lin3 <- as.data.frame(fit_lin3)

colnames(fit_lin3) <- names(phthalates_met2[,phstart:phend])
colnames(fit_lin4) <- names(phthalates_met2[,phstart:phend])
colnames(fit_lin5) <- names(phthalates_met2[,phstart:phend])

fit_lin3 <- cbind.data.frame(fit_lin3, Metabolite = names(phthalates_met2[,metstart:metend]))
fit_lin4 <- cbind.data.frame(fit_lin4, Metabolite = names(phthalates_met2[,metstart:metend]))
fit_lin5 <- cbind.data.frame(fit_lin5, Metabolite = names(phthalates_met2[,metstart:metend]))



# Now, I will derive FDR adjusted p-values again. Using Q value approach would result in same conclusions regarding the significant metabolites.

# ```{r checking for significant metabolites after adding covariates}

#Define FDR threshold value

pthres = 0.050
# pthres = 0.100

fit_lin3_adjust <- cbind.data.frame(Metabolite = fit_lin3[,11], 
                                    apply(fit_lin3[,1:10], 2, function(x) p.adjust(x, method = "BH")))

origpvalue_cord <- cbind.data.frame(Metabolite = fit_lin3[,11], fit_lin3[,1:10])

fit_lin3_adjust_long <- fit_lin3_adjust %>% pivot_longer(!Metabolite, names_to = "Phthalates", values_to = "adj_p") %>% 
  mutate(pval_ind = ifelse(adj_p > pthres, 0, 1))

fit_lin4_adjust_long <- fit_lin4 %>% pivot_longer(!Metabolite, names_to = "Phthalates", values_to = "betas") 
fit_lin5_adjust_long <- fit_lin5 %>% pivot_longer(!Metabolite, names_to = "Phthalates", values_to = "se") 


fit_lin34_adjust_long <- cbind.data.frame(fit_lin3_adjust_long, betas = fit_lin4_adjust_long$betas, se = fit_lin5_adjust_long$se)
fit_lin34_adjust_long <- fit_lin34_adjust_long %>% 
  mutate(beta_sign = ifelse(betas < 0, "Negative", "Positive")) %>%
  mutate(Phthalate =  stringi::stri_replace_all_regex(Phthalates, pattern = c("log_", "_GM", "SG_"), replacement = c("", ""), vectorize_all = FALSE)) %>% select(-Phthalates)



########################################################################################################################################################################

# Now draw a Manhattan plot of p-values


group.colors <- c(Negative = "skyblue4", Positive = "tomato1")


geom_text_data2 <- fit_lin34_adjust_long %>% filter(pval_ind == 1)

fit_lin34_adjust_long <- fit_lin34_adjust_long %>% mutate(`Beta coefficient` = beta_sign, 
                                                          Label = ifelse(adj_p <= pthres, Metabolite, ""), 
                                                          Metabolite = toupper(Metabolite), Phthalate = toupper(Phthalate))




fit_lin34_adjust_long <- fit_lin34_adjust_long %>% 
  mutate(adj_p = round(adj_p, 3)) %>%
  mutate(`Beta coefficient` = beta_sign, Label = ifelse(adj_p <= pthres, "1", "")) 

fit_lin34_adjust_long <- fit_lin34_adjust_long %>% mutate(Phthalate = ifelse(Phthalate == "MBZP", "MBzP", Phthalate))




# fit_lin34_adjust_long <- fit_lin34_adjust_long %>% 
#   mutate(Phthalate2 = substr(Phthalate, 4, nchar(Phthalate)))

man_plot2 <- ggplot(fit_lin34_adjust_long, aes(x = Metabolite, y = -log10(adj_p), color = `Beta coefficient`)) + 
  geom_point(shape = 15) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(fill = 'white'), strip.placement = "outside", axis.line = element_line(),
        legend.position = "bottom") +
  scale_y_continuous(expression(paste(-log10('p-value')))) +
  scale_color_discrete("Beta coefficient") +
  scale_color_manual(values = c("Positive" = "skyblue4", "Negative" = "red")) +
  facet_wrap(~Phthalate, strip.position = "bottom", nrow = 1) 


man_plot2 <- man_plot2  + theme(panel.spacing.x = unit(0, "cm")) +
  geom_hline(yintercept = -log10(0.1), col = "blue", linetype = 2) +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 2) 
  # ggtitle("P-values from MLR models for each phthalate metabolite and \nmetabolic marker") 


man_plot2

########################################################################################################################################################################

# Now, I will tabulate and plot only the significant metabolites

fit_lin34_adjust_long <- fit_lin34_adjust_long %>% 
  rename("Phthalates" = "Phthalate", "Metabolites" = "Metabolite")

all_sig2 <- fit_lin34_adjust_long %>% mutate(adj_p = round(adj_p, 3)) %>%
  filter(adj_p <= pthres) %>% arrange(Phthalates, adj_p) %>% 
  select(Phthalates, Metabolites, adj_p, betas, se, beta_sign) #2 other metabolic markers with adj p-values = 0.06

# all_sig2 %>% group_by(Phthalates) %>% dplyr::summarize(Nsig = n())

ggplot(all_sig2, aes(x = Phthalates, group = Phthalates, fill = Phthalates)) + geom_bar() + 
  theme(legend.position = "none", axis.title.x  = element_blank()) + 
  ggtitle("Number of significant metabolic markers per phthalate metabolite")  + 
  facet_wrap(~beta_sign)


########################################################################################################################################################################
# Plot significant metabolites using a heat plot. Here, I am only plotting metabolites with 
# adjusted p-values \< 0.1 to make it easier to visuaulize the result.


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 50)

all_sig2_plot <- all_sig2  %>% mutate(`Beta coefficient` = betas) 

# Create the heatmap
# all_sig2_plot <- all_sig2_plot %>%  
#   mutate(Metabolites2 = case_when(Metabolites == "AMINO FATTY ACIDS (C8H17NO2)" ~ "2-Aminooctanoic acid",
#                                   Metabolites == "DICARBOXYLIC ACIDS (C16H30O4)" ~ "Hexadecanedioic acid",
#                                   Metabolites == "DICARBOXYLIC ACIDS (C12H22O4)" ~ "Dodecanedioic acid",
#                                   Metabolites == "DICARBOXYLIC ACIDS (C18H32O4)" ~ "Octadecenedioic acid",
#                                   Metabolites == "FA(18:1(Ke)).1" ~ "FA(18:1(Ke))"))
# 
# all_sig2_plot <- all_sig2_plot %>% mutate(Metabolites = ifelse(is.na(Metabolites2), Metabolites, Metabolites2)) %>%
#   mutate(Metabolites = toupper(Metabolites))
# 


#Check # of significantly associated phthalates per metabolite


sum_all_sig2 <- all_sig2 %>% group_by(Metabolites) %>% summarise(NumPhthalates = n()) %>% 
  arrange(desc(NumPhthalates))
sum_all_sig2

########################################################################################################################################################################

#Add classes of metabolites


###Add classes of significant metabolites

all_sig2_plot$Metabolites <- trimws(all_sig2_plot$Metabolites, which = c("both"))
all_sig2_plot$Metabolites <- toupper(all_sig2_plot$Metabolites)
# all_sig2_plot <- all_sig2_plot %>% select(-Metabolites2)


##Change names of some metabolites in Classes file too

classes <- classes %>% mutate(Metabolites2 = case_when(Metabolites == "AMINO FATTY ACIDS (C8H17NO2)" ~ "2-Aminooctanoic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C16H30O4)" ~ "Hexadecanedioic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C12H22O4)" ~ "Dodecanedioic acid",
                                                       Metabolites == "DICARBOXYLIC ACIDS (C18H32O4)" ~ "Octadecenedioic acid",
                                                       Metabolites == "FA(18:1(Ke)).1" ~ "FA(18:1(Ke))"))
classes <- classes %>% mutate(Metabolites = ifelse(is.na(Metabolites2), Metabolites, Metabolites2)) %>%
  mutate(Metabolites = toupper(Metabolites))

all_sig2_class <- all_sig2_plot %>% left_join(classes, by = c("Metabolites"))

all_sig2_class2 <- all_sig2_class 

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
  labs(x = "Phthalates", y = "Cord blood metabolic features")



#Also export significant file 

sigfile_to_export <- all_sig2_class2 %>%
  mutate(
    lower_ci = round(betas - 1.96 * se, 2),
    upper_ci = round(betas + 1.96 * se, 2),
    beta_ci = sprintf("%.2f (%.2f, %.2f)", betas, lower_ci, upper_ci), # Using sprintf to format the string
    adj_p = round(adj_p, 3)) %>%  # Round adjusted p-values to 3 decimal points
  select(Phthalates, Metabolites, Class_new2, beta_ci, adj_p) %>% arrange(Phthalates, Class_new2) 


write.csv(sigfile_to_export, file = "C:/Users/rsiwa/OneDrive/Michigan_course/research_rotation2/manuscript_writing/significant_metabolites_by_phthalates_cordblood.csv")

#Look at original p-value 

origpvalue_cord2 <- origpvalue_cord %>% 
  rename("Metabolites" = "Metabolite") %>% 
  pivot_longer(cols = log_SG_MBzP_GM:log_SG_DEHTP_GM, names_to = "Phthalates", values_to = "Unadj_p") %>%
  mutate(Metabolites = toupper(Metabolites), Phthalates = toupper(Phthalates),  
  Phthalates = stringi::stri_replace_all_regex(Phthalates, pattern = c("LOG_", "_GM", "SG_"), replacement = c("", ""), vectorize_all = FALSE)) %>%
  inner_join(sigfile_to_export, by = c("Metabolites", "Phthalates"))


#Add class to a file with all metabolites
fit_lin34_adjust_long_chk <- fit_lin34_adjust_long %>% mutate(Metabolites = toupper(Metabolites)) %>%
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

mets_forCor <- phthalates_met2[, 4:276]
cor_met <- rcorr(as.matrix(mets_forCor), type="spearman")
cor_met_sp <- cor_met$r
cor_met_sp2 = 1-cor_met_sp 

#Hierarchial clustering

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 50)
# heatmap.2(cor_met_sp, dendrogram="col",  density.info="none", trace="none",
#           main="Metabolite Correlations", margins =c(12,9), col=my_palette)


#Function to create a cluster

getclus <- function(x, data){
  hc.rows <- hclust((data))
  ct <- cutree(hc.rows, h=x)  # cut the dendrogram into clusters
  plot(hc.rows,hang = -1, cex = 0.4) #plot dendrogram
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
    cor_mean <- mean(gp_cor[lower.tri(gp_cor)],na.rm=TRUE)
    row=data.frame(group=i, num=length(gp),cor_mean=round(cor_mean,3))
    mean_all <- rbind(mean_all, row)
  }
  return(mean_all)
}

########################################################################################################################################################################

#See phthalates_metabolites_62023 for details. However, we will use h = 1.1

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


plot(matval$len, matval$n_clus)
abline(v = 1.1)

#Get clusters with at least 5 metabolites and r < 0.3
set.seed(1991)
cluslen <- 1.1
optclus <- getclus(cluslen, data = as.dist(cor_met_sp2))
optclus_mean <- get_mean(optclus)

########################################################################################################################################################################


#Better way to visualize cluster 

hc.rows <- hclust((as.dist(cor_met_sp2)))
plot(hc.rows, hang = -1, main = "", cex = 0.3, sub = "", xlab = "") 
dend <- rect.hclust(hc.rows, h = 1.1) # visualize the groups
title(main = "Hierarchical clustering for cord blood metabolic features", cex.main = 1)

clust_cord <- cutree(hc.rows, h = 1.1)
clust_cord
table(clust_cord)

#cluster 1, 10, 11, 15, 16 do not meet criteria 

########################################################################################################################################################################


corchk <- which(optclus_mean$cor_mean <0.3)
sizechk <- which(optclus_mean$num <5)

metabolites_t3 <- phthalates_met2[,c(2,4:276)]
t_metabolites_t3 <- data.table::transpose(metabolites_t3)
rownames(t_metabolites_t3) <- colnames(metabolites_t3)
colnames(t_metabolites_t3) <- t_metabolites_t3[1,]
t_metabolites_t4 <- t_metabolites_t3[-1,]
t_metabolites_t4 <- t_metabolites_t4 %>% mutate(metabolite = rownames(t_metabolites_t4), Cluslabel = optclus[[2]]) %>% 
  select(metabolite, Cluslabel, everything())


########################################################################################################################################################################


###Multivariate regression using the cluster memberships

# Following Margaret's code

#list of clusterid
clusterid <- t_metabolites_t4$Cluslabel

#number of metabolites in a cluster
clustersum <- aggregate(rep(1,nrow(t_metabolites_t4)) ~ clusterid, FUN = sum)

#Choose the clusters to use 

clustertouse <- optclus_mean %>% filter(num >=5, cor_mean >= 0.3) 


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

write.csv(clus_chk1_toexport, "cordblood_metabolites_cluster_assignment.csv")


ttmetabolites_t4 <- as.data.frame(t(t_metabolites_t4))
ttmetabolites_t5 <- ttmetabolites_t4[3:nrow(ttmetabolites_t4),]


names(phthalates_met2) <- toupper(names(phthalates_met2))


row_all2 <- NULL

for (j in phstart:phend){
  z = data.frame(Phth=phthalates_met2[,j] , AGE = phthalates_met2$ISAGE, BMI=phthalates_met2$PREBMI,
                 SEX = phthalates_met2$PPSEX)
  zz = data.matrix(z)
  
  for(i in clustertouse[,1]){
    ik = i
    x = ttmetabolites_t5[, clusterid == ik]
    xx = x
    xxd = data.frame(xx)
    pp = ncol(xxd)
    names(xxd) = c(paste0("M", 1:pp))
    
    data1 = data.frame(Phth = zz[,1], AGE = scale(zz[,2]), BMI=zz[,3], SEX = zz[,4], xxd)
    yvec <- as.numeric(as.matrix(cbind(data1[,5:ncol(data1)])))
    yvec <- matrix(yvec, ncol = ncol(data1)-4)
    xvec <- as.matrix(cbind(data1[,1:4]))
    my.model <- lm(yvec ~ xvec) 
    
    xvec2 <- as.matrix(cbind(data1[,2:4]))
    mlm2 <- lm(yvec ~ xvec2)
    pval <- as.numeric(anova(my.model, mlm2)$'Pr(>F)'[2])
    
    row <- cbind(Metab_Cluster=ik, Phthalate = colnames(phthalates_met2[j]),pval = pval, cluster_size = pp)
    
    row_all2 <- rbind(row_all2, row)
    # row_all2 <- as.data.frame(row_all2) %>% filter(pval < 0.05)
  }
}

pthres2 = 0.050
# Now look at the phthalate metabolites that were significant after adjusting for multiple testing.

row_all2_bhcor_chk <- data.frame(row_all2) %>% 
  mutate(pval = as.numeric(as.character(pval))) %>% 
  group_by(Phthalate) %>% 
  mutate(BH = p.adjust(pval, method = "BH"), BH_round = round(BH, 3))

row_all2_bhcor <- row_all2_bhcor_chk %>% 
 filter(BH <= pthres2) #This is done to make sure that the cluster with p-value of 0.05 is not excluded. 
  # dplyr::mutate(q = qvalue(pval)$qvalues)  # Not enough tests per phthalate for qvalue

row_all2_bhcor_chk %>% arrange(as.numeric(Metab_Cluster))

#export the results file 
row_all2_bhcor_chk_wide <- row_all2_bhcor_chk %>% select(Metab_Cluster, Phthalate, BH_round) %>% 
  pivot_wider(names_from = Phthalate, values_from = BH_round)

write.csv(row_all2_bhcor_chk_wide, "hierarchical_cluster_results_cordblood.csv")

########################################################################################################################################################################


# Now, add class label to significantly associated metabolites


row_all2_bhcor_sum <- row_all2_bhcor %>% group_by(Metab_Cluster) %>% tally() 

met_clus <- t_metabolites_t4 %>% select(metabolite, Cluslabel) %>% dplyr::rename("Metab_Cluster" = Cluslabel, "Metabolites" = metabolite) %>%
  mutate(Metab_Cluster = as.character(Metab_Cluster))

met_clus <- met_clus %>% mutate(Metabolites = toupper(Metabolites))


met_clus2 <- met_clus %>% inner_join(row_all2_bhcor_sum, by = "Metab_Cluster")

# met_clus2 <- met_clus2 %>% mutate(Metabolites2 = case_when(Metabolites == "AMINO FATTY ACIDS (C8H17NO2)" ~ "2-Aminooctanoic acid",
#                                                        Metabolites == "DICARBOXYLIC ACIDS (C16H30O4)" ~ "Hexadecanedioic acid",
#                                                        Metabolites == "DICARBOXYLIC ACIDS (C12H22O4)" ~ "Dodecanedioic acid",
#                                                        Metabolites == "DICARBOXYLIC ACIDS (C18H32O4)" ~ "Octadecenedioic acid",
#                                                        Metabolites == "FA(18:1(Ke)).1" ~ "FA(18:1(Ke))"))
# 
# met_clus2 <- met_clus2 %>% mutate(Metabolites = ifelse(is.na(Metabolites2), Metabolites, Metabolites2)) %>%
#   mutate(Metabolites = toupper(Metabolites))

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
met_clus2$Metab_Cluster3 <- factor(met_clus2$Metab_Cluster2, levels = c("Cluster 3", "Cluster 4"))


ggplot(data = met_clus2, aes(y = Class_new2, fill = Class_new2)) +
  geom_bar() +
  facet_wrap(~as.factor(Metab_Cluster3), nrow = 1) +
  theme_classic() +
  ylab("Metabolite classes") +
  theme(
    legend.position = "right", legend.title = element_blank(),
    axis.text.y = element_text(size = 7),
    strip.background = element_blank(), # Remove the box around the facet name
    strip.text = element_text(vjust = 0), # Adjust text positioning if needed
    panel.border = element_rect(colour = "black", fill = NA, size = 1) # Add axis line to each facet
  )

  
  #Need to change the names of classes. For example, purines and derivatives should be purines and purine derivatives to be consistent with plasma metabolic data




########################################################################################################################################################################


# Now combine results from cluster and single pollutant analysis 

#Get rid of duplicate metabolites

met_clus2 <- met_clus2 %>% filter(substr(Metabolites, nchar(Metabolites), nchar(Metabolites)) != 1)

met_clus2b <- met_clus2 %>% select(Metabolites, Metab_Cluster, n, Class_new2)
all_sig2_class2b <- all_sig2_class2 %>% select(Metabolites, Phthalates, adj_p, beta_sign)
fit_lin34_adjust_longb <- fit_lin34_adjust_long %>% select(Metabolites, Phthalates, adj_p, beta_sign, pval_ind)
fit_lin34_adjust_longb <- fit_lin34_adjust_longb  %>% 
  mutate(sig_sign = paste(pval_ind, beta_sign, sep = "_")) %>% select(-adj_p, -beta_sign, -pval_ind)
fit_lin34_adjust_longb_wide <- pivot_wider(fit_lin34_adjust_longb, names_from = Phthalates, values_from = sig_sign) %>%
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

write.csv(met_clus_ind3, "cluster3.csv")
write.csv(met_clus_ind4, "cluster4.csv")
# write.csv(met_clus_ind16, "cluster16.csv")
# write.csv(met_clus_ind19, "cluster19.csv")
# write.csv(met_clus_ind5, "cluster5.csv")
# write.csv(met_clus_ind6, "cluster6.csv")


#Combine three clusters
# met_clus_ind_top3 <- rbind.data.frame(met_clus_ind10, met_clus_ind1, met_clus_ind13)


# Now, let us plot metabolites from each cluster for Figure 2


sigmet_clus3 <- c("DEHTP")
sigmet_clus4 <- c("DEHTP")
# sigmet_clus16 <- c("MCPP")
# sigmet_clus19 <- c("MCPP")
# sigmet_clus5 <- c("MCPP")
# sigmet_clus6 <- c("DEHP")



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
    ggtitle(title) + theme(axis.title  = element_blank(), legend.position = "none", axis.ticks = element_blank())
  
  return(plot1)
  
}


clus3_plot <- datprep_tile(data = met_clus_ind3, vars = sigmet_clus3, title = "Cluster 3")
clus3_plot

clus4_plot <- datprep_tile(data = met_clus_ind4, vars = sigmet_clus4, title = "Cluster 4") 
clus4_plot
# 
# clus16_plot <- datprep_tile(data = met_clus_ind16, vars = sigmet_clus16, title = "Cluster 16") 
# clus16_plot
# 
# clus19_plot <- datprep_tile(data = met_clus_ind19, vars = sigmet_clus19, title = "Cluster 19")
# clus19_plot
# 
# clus5_plot <- datprep_tile(data = met_clus_ind5, vars = sigmet_clus5, title = "Cluster 5") 
# clus5_plot
# 
# clus6_plot <- datprep_tile(data = met_clus_ind6, vars = sigmet_clus6, title = "Cluster 6") 
# clus6_plot

# text <- "Table 2. List of metabolic markers in significantly associated clusters. * represents metabolic markers that were significantly associated in MWAS."
# text.p <- ggpubr::ggparagraph(text = text, face = "italic", size = 11, color = "black")
# 
ggpubr::ggarrange(clus3_plot, clus4_plot, align = "v")

#################################################################################################################################
#The next analytical approach is to look at correlation network. I will export data to do that


metlist <- data.frame(Metabolites = toupper(names(phthalates_met2[, 4:276])))

classes_toexport <- classes %>%  select(-Metabolites2) %>% mutate(Metabolites = toupper(Metabolites))
classes_toexport <- classes_toexport %>% inner_join(metlist, by = c("Metabolites"))
# write.csv(classes_toexport, "classes_cordblood_metfile.csv")


#################################################################################################################################

#To get %> LOD for cord blood data 

cord_orig_data <- read.csv("C:\\Users\\rsiwa\\Dropbox (University of Michigan)\\Ram_shared\\data\\2_v54_006_Phthalates_no_dup_03092021_20210625153132.csv") 
cord_orig_data$StudyID <- as.character(cord_orig_data$StudyID)

cord_orig_data2 <- phthalates_met2 %>% select(`Study ID`) %>% rename(StudyID = "Study ID") %>% 
  inner_join(cord_orig_data, by = "StudyID")



#Detection rate 

table(cord_orig_data2$VISITID, cord_orig_data2$ANALYTE)

blod_summary1 <- cord_orig_data2 %>% group_by(ANALYTE, VISITID) %>% dplyr::summarize(detrate = 100-sum(BLOD)/81*100) %>% dplyr::filter(VISITID == 1)
blod_summary2 <- cord_orig_data2 %>% group_by(ANALYTE, VISITID) %>% dplyr::summarize(detrate = 100-sum(BLOD)/79*100) %>% dplyr::filter(VISITID == 2)
blod_summary3 <- cord_orig_data2 %>% group_by(ANALYTE, VISITID) %>% dplyr::summarize(detrate = 100-sum(BLOD)/71*100) %>% dplyr::filter(VISITID == 3)

blod_summary <- rbind.data.frame(blod_summary1, blod_summary2, blod_summary3) 

blod_summary_chk <- cord_orig_data2 %>% group_by(ANALYTE) %>% dplyr::summarize(detrate = 100-sum(BLOD)/231*100, nasum = sum(is.na(BLOD))/231*100) 



ggplot(blod_summary, aes(y = ANALYTE, x = detrate/3, fill = ANALYTE)) + 
  geom_bar(stat = "identity") + theme(legend.position = "none") + xlab("det rate") + ylab("Phthalates") +
  ggtitle("Number of detects for each Phthalate metoabolite")


#Putting data in wide format

cord_orig_data_wide <- cord_orig_data2 %>% select(StudyID, VISITID, ANALYTE, BLOD, CONC) %>% group_by(StudyID, ANALYTE) %>%
  dplyr::summarize(ConcGM = exp(mean(log(CONC), na.rm = TRUE))) %>% 
  pivot_wider(names_from = c("ANALYTE"), # The columns to use as column names
                         values_from = c("ConcGM")) # The columns to use as cell values

#Distribution of classes

class_sum_cord <- classes %>% group_by(Class_new2) %>% tally()






