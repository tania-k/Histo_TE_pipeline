#Large repetative and messy script which will help you collect and assess your TE results to a stacked barplot visual.
#First calculate TEs, collect all varieties per isolate, combine into large table
#I worked the final table on my computer, had the script read that in to do final TE analysis.
#Trees are built here too, but had two CDS and AA trees, went with the AA tree.
#Created trees using Aspergillus fumigatus, Emmonsia and Blastomyces as the root species. Final tree uses Blastomyces as root.
#Also added Phylogenetic signal to each rooted species. Make sure to run the phylogenetic signal analysis for each TE type.
#I had 7 groups done, RNA TEs, DNA TEs, Simple Repeats, Low Complexity, MITEs, Satellites, and Unknowns. I've seen many disregard Low complexity, simple repeats and unknowns. Up to your discretion. 


library(tidyverse)
library(readr)
library(plyr)
library(dplyr)
library(ggsci)
library(gridExtra)
library(ggplot2)
library(hrbrthemes)

setwd("/nas/longleaf/home/taniak/taniak/TE_Analysis")

Histoplasma_capsulatum_WU24 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_WU24.RM/Histoplasma_capsulatum_WU24.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_WU24 ) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_WU24  <- Histoplasma_capsulatum_WU24  %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_WU24,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_WU24.tsv")

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
Histoplasma_capsulatum_WU24class <- Histoplasma_capsulatum_WU24 %>% group_by(class) %>% 
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/32531515)
write_tsv(Histoplasma_capsulatum_WU24class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Histoplasma_capsulatum_WU24class.tsv")
Histoplasma_capsulatum_WU24class <- Histoplasma_capsulatum_WU24class %>% mutate('Histoplasma capsulatum WU24'= percentage)

#Histoplasma_capsulatum_WU24class_N <- Histoplasma_capsulatum_WU24class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_WU24class1 <- Histoplasma_capsulatum_WU24class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_G217B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_G217B.RM/Histoplasma_ohiense_G217B.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_G217B ) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_G217B  <- Histoplasma_ohiense_G217B  %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_G217B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_ohiense_G217B.tsv")

Histoplasma_ohiense_G217Bclass <- Histoplasma_ohiense_G217B %>% group_by(class) %>% 
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/39447273)
write_tsv(Histoplasma_ohiense_G217Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Histoplasma_ohiense_G217Bclass.tsv")
Histoplasma_ohiense_G217Bclass <- Histoplasma_ohiense_G217Bclass %>% mutate('Histoplasma ohiense G217B'= percentage)
#Histoplasma_capsulatum_WU24class_N <- Histoplasma_capsulatum_WU24class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_G217Bclass1 <- Histoplasma_ohiense_G217Bclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_var.duboisii_H88 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_var.duboisii_H88.RM/Histoplasma_capsulatum_var.duboisii_H88.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_var.duboisii_H88 ) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_var.duboisii_H88  <- Histoplasma_capsulatum_var.duboisii_H88  %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_var.duboisii_H88,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_var.duboisii_H88.tsv")

Histoplasma_capsulatum_var.duboisii_H88class <- Histoplasma_capsulatum_var.duboisii_H88 %>% group_by(class) %>% 
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/37996987)
write_tsv(Histoplasma_capsulatum_var.duboisii_H88class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Histoplasma_capsulatum_var.duboisii_H88class.tsv")
Histoplasma_capsulatum_var.duboisii_H88class <- Histoplasma_capsulatum_var.duboisii_H88class %>% mutate('Histoplasma capsulatum var.duboisii H88'= percentage)
Histoplasma_capsulatum_var.duboisii_H88class1 <- Histoplasma_capsulatum_var.duboisii_H88class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_G186AR <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_G186AR.RM/Histoplasma_capsulatum_G186AR.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_G186AR ) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_G186AR  <- Histoplasma_capsulatum_G186AR  %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_G186AR,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_G186AR.tsv")

Histoplasma_capsulatum_G186ARclass <- Histoplasma_capsulatum_G186AR %>% group_by(class) %>% 
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/31111494)
write_tsv(Histoplasma_capsulatum_G186ARclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Histoplasma_capsulatum_G186ARclass.tsv")
Histoplasma_capsulatum_G186ARclass <- Histoplasma_capsulatum_G186ARclass %>% mutate('Histoplasma capsulatum G186AR'= percentage)
Histoplasma_capsulatum_G186ARclass1 <- Histoplasma_capsulatum_G186ARclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_G184AR <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_G184AR.RM/Histoplasma_capsulatum_G184AR.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_G184AR) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_G184AR  <- Histoplasma_capsulatum_G184AR  %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_G184AR,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_G184AR.tsv")

Histoplasma_capsulatum_G184ARclass <- Histoplasma_capsulatum_G184AR %>% group_by(class) %>% 
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/30991520)
write_tsv(Histoplasma_capsulatum_G184ARclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Histoplasma_capsulatum_G184ARclass.tsv")
Histoplasma_capsulatum_G184ARclass <- Histoplasma_capsulatum_G184ARclass %>% mutate('Histoplasma capsulatum G184AR'= percentage)
Histoplasma_capsulatum_G184ARclass1 <- Histoplasma_capsulatum_G184ARclass %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_HCCI_17 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_HCCI_17.RM/Histoplasma_ohiense_HCCI_17.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_HCCI_17) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_HCCI_17 <- Histoplasma_ohiense_HCCI_17 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_HCCI_171 <- subset(Histoplasma_ohiense_HCCI_17,class != "Low_complexity")
#Histoplasma_ohiense_HCCI_172 <- subset(Histoplasma_ohiense_HCCI_171,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_HCCI_17,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_HCCI_17_version2.tsv")

Histoplasma_ohiense_HCCI_17class<- Histoplasma_ohiense_HCCI_17 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/34006939)
write_tsv(Histoplasma_ohiense_HCCI_17class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_HCCI_17class_version2.tsv")
Histoplasma_ohiense_HCCI_17class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCCI_17class_version2.tsv")
Histoplasma_ohiense_HCCI_17class <- Histoplasma_ohiense_HCCI_17class %>% mutate('Histoplasma ohiense HCCI_17'= percentage)
Histoplasma_ohiense_HCCI_17class_N <- Histoplasma_ohiense_HCCI_17class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_HCCI_17class1 <- Histoplasma_ohiense_HCCI_17class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_HCCI_6 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_HCCI_6.RM/Histoplasma_ohiense_HCCI_6.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_HCCI_6) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_HCCI_6 <- Histoplasma_ohiense_HCCI_6 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_HCCI_61 <- subset(Histoplasma_ohiense_HCCI_6,class != "Low_complexity")
#Histoplasma_ohiense_HCCI_62 <- subset(Histoplasma_ohiense_HCCI_61,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_HCCI_6,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_HCCI_6.tsv")

Histoplasma_ohiense_HCCI_6class<- Histoplasma_ohiense_HCCI_62 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28445673)
write_tsv(Histoplasma_ohiense_HCCI_6class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_HCCI_6class.tsv")
Histoplasma_ohiense_HCCI_6class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCCI_6class.tsv")
Histoplasma_ohiense_HCCI_6class <- Histoplasma_ohiense_HCCI_6class %>% mutate('Histoplasma ohiense HCCI_6'= percentage)
Histoplasma_ohiense_HCCI_6class1 <- Histoplasma_ohiense_HCCI_6class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_HCG217B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_HCG217B.RM/Histoplasma_ohiense_HCG217B.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_HCG217B) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_HCG217B <- Histoplasma_ohiense_HCG217B %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_HCG217B1 <- subset(Histoplasma_ohiense_HCG217B,class != "Low_complexity")
#Histoplasma_ohiense_HCG217B2 <- subset(Histoplasma_ohiense_HCG217B1,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_HCG217B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_HCG217B.tsv")

Histoplasma_ohiense_HCG217Bclass<- Histoplasma_ohiense_HCG217B2 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/33872679)
write_tsv(Histoplasma_ohiense_HCG217Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCG217Bclass.tsv")
Histoplasma_ohiense_HCG217Bclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCG217Bclass.tsv")
Histoplasma_ohiense_HCG217Bclass <- Histoplasma_ohiense_HCG217Bclass %>% mutate('Histoplasma ohiense HCG217B'= percentage)
Histoplasma_ohiense_HCG217Bclass1 <- Histoplasma_ohiense_HCG217Bclass %>% select(-c(n,total,percentage))
Histoplasma_ohiense_HCG217Bclass_N <- Histoplasma_ohiense_HCG217Bclass %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_82.Nam1_CI_4 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_82-Nam1_CI_4.RM/Histoplasma_ohiense_SECH_82-Nam1_CI_4.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_82.Nam1_CI_4) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_82.Nam1_CI_4 <- Histoplasma_ohiense_SECH_82.Nam1_CI_4 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_82.Nam1_CI_41 <- subset(Histoplasma_ohiense_SECH_82.Nam1_CI_4,class != "Low_complexity")
#Histoplasma_ohiense_SECH_82.Nam1_CI_42 <- subset(Histoplasma_ohiense_SECH_82.Nam1_CI_41,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_82.Nam1_CI_4,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_82-Nam1_CI_4.tsv")

Histoplasma_ohiense_SECH_82.Nam1_CI_4class<- Histoplasma_ohiense_SECH_82.Nam1_CI_42 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/30880349)
write_tsv(Histoplasma_ohiense_SECH_82.Nam1_CI_4class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_82-Nam1_CI_4class.tsv")
Histoplasma_ohiense_SECH_82.Nam1_CI_4class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_82-Nam1_CI_4class.tsv")
Histoplasma_ohiense_SECH_82.Nam1_CI_4class <- Histoplasma_ohiense_SECH_82.Nam1_CI_4class %>% mutate('Histoplasma ohiense SECH_82-Nam1_CI_4'= percentage)
Histoplasma_ohiense_SECH_82.Nam1_CI_4class1 <- Histoplasma_ohiense_SECH_82.Nam1_CI_4class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_82.Nam1_CI_4class_N <- Histoplasma_ohiense_SECH_82.Nam1_CI_4class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_91.Nam2_G217B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_91-Nam2_G217B.RM/Histoplasma_ohiense_SECH_91-Nam2_G217B.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_91.Nam2_G217B) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_91.Nam2_G217B <- Histoplasma_ohiense_SECH_91.Nam2_G217B %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_91.Nam2_G217B1 <- subset(Histoplasma_ohiense_SECH_91.Nam2_G217B,class != "Low_complexity")
#Histoplasma_ohiense_SECH_91.Nam2_G217B2 <- subset(Histoplasma_ohiense_SECH_91.Nam2_G217B1,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_91.Nam2_G217B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_91-Nam2_G217B.tsv")

Histoplasma_ohiense_SECH_91.Nam2_G217Bclass<- Histoplasma_ohiense_SECH_91.Nam2_G217B2 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29816969)
write_tsv(Histoplasma_ohiense_SECH_91.Nam2_G217Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_91-Nam2_G217Bclass.tsv")
Histoplasma_ohiense_SECH_91.Nam2_G217Bclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_91-Nam2_G217Bclass.tsv")
Histoplasma_ohiense_SECH_91.Nam2_G217Bclass <- Histoplasma_ohiense_SECH_91.Nam2_G217Bclass %>% mutate('Histoplasma ohiense SECH_91-Nam2_G217B'= percentage)
Histoplasma_ohiense_SECH_91.Nam2_G217Bclass1 <- Histoplasma_ohiense_SECH_91.Nam2_G217Bclass %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_91.Nam2_G217Bclass_N <- Histoplasma_ohiense_SECH_91.Nam2_G217Bclass %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_92.Nam2_G222B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_92-Nam2_G222B.RM/Histoplasma_ohiense_SECH_92-Nam2_G222B.scaffolds.fa.out",skip=3,col_names = F)

colnames(Histoplasma_ohiense_SECH_92.Nam2_G222B) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_92.Nam2_G222B <- Histoplasma_ohiense_SECH_92.Nam2_G222B %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_92.Nam2_G222B1 <- subset(Histoplasma_ohiense_SECH_92.Nam2_G222B,class != "Low_complexity")
#Histoplasma_ohiense_SECH_92.Nam2_G222B2 <- subset(Histoplasma_ohiense_SECH_92.Nam2_G222B1,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_92.Nam2_G222B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_92-Nam2_G222B.tsv")

Histoplasma_ohiense_SECH_92.Nam2_G222Bclass<- Histoplasma_ohiense_SECH_92.Nam2_G222B2 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29639711)
write_tsv(Histoplasma_ohiense_SECH_92.Nam2_G222Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_92-Nam2_G222Bclass.tsv")
Histoplasma_ohiense_SECH_92.Nam2_G222Bclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_92-Nam2_G222Bclass.tsv")
Histoplasma_ohiense_SECH_92.Nam2_G222Bclass <- Histoplasma_ohiense_SECH_92.Nam2_G222Bclass %>% mutate('Histoplasma ohiense SECH_92-Nam2_G222B'= percentage)
Histoplasma_ohiense_SECH_92.Nam2_G222Bclass1 <- Histoplasma_ohiense_SECH_92.Nam2_G222Bclass %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_92.Nam2_G222Bclass_N <- Histoplasma_ohiense_SECH_92.Nam2_G222Bclass %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_93.Nam2_CI_6 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_93-Nam2_CI_6.RM/Histoplasma_ohiense_SECH_93-Nam2_CI_6.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_93.Nam2_CI_6) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_93.Nam2_CI_6 <- Histoplasma_ohiense_SECH_93.Nam2_CI_6 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_93.Nam2_CI_61 <- subset(Histoplasma_ohiense_SECH_93.Nam2_CI_6,class != "Low_complexity")
#Histoplasma_ohiense_SECH_93.Nam2_CI_62 <- subset(Histoplasma_ohiense_SECH_93.Nam2_CI_61,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_93.Nam2_CI_6,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_93-Nam2_CI_6.tsv")

Histoplasma_ohiense_SECH_93.Nam2_CI_6class<- Histoplasma_ohiense_SECH_93.Nam2_CI_62 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23646057)
write_tsv(Histoplasma_ohiense_SECH_93.Nam2_CI_6class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_93-Nam2_CI_6class.tsv")
Histoplasma_ohiense_SECH_93.Nam2_CI_6class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_93-Nam2_CI_6class.tsv")
Histoplasma_ohiense_SECH_93.Nam2_CI_6class <- Histoplasma_ohiense_SECH_93.Nam2_CI_6class %>% mutate('Histoplasma ohiense SECH_93.Nam2_CI_6'= percentage)
Histoplasma_ohiense_SECH_93.Nam2_CI_6class1 <- Histoplasma_ohiense_SECH_93.Nam2_CI_6class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_93.Nam2_CI_6class_N <- Histoplasma_ohiense_SECH_93.Nam2_CI_6class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_94.Nam2_CI_9 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_94-Nam2_CI_9.RM/Histoplasma_ohiense_SECH_94-Nam2_CI_9.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_94.Nam2_CI_9) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_94.Nam2_CI_9 <- Histoplasma_ohiense_SECH_94.Nam2_CI_9 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_94.Nam2_CI_91 <- subset(Histoplasma_ohiense_SECH_94.Nam2_CI_9,class != "Low_complexity")
#Histoplasma_ohiense_SECH_94.Nam2_CI_92 <- subset(Histoplasma_ohiense_SECH_94.Nam2_CI_91,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_94.Nam2_CI_9,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_94.Nam2_CI_9.tsv")

Histoplasma_ohiense_SECH_94.Nam2_CI_9class<- Histoplasma_ohiense_SECH_94.Nam2_CI_92 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/20233006)
write_tsv(Histoplasma_ohiense_SECH_94.Nam2_CI_9class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_94.Nam2_CI_9class.tsv")
Histoplasma_ohiense_SECH_94.Nam2_CI_9class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_94.Nam2_CI_9class.tsv")
Histoplasma_ohiense_SECH_94.Nam2_CI_9class <- Histoplasma_ohiense_SECH_94.Nam2_CI_9class %>% mutate('Histoplasma ohiense SECH_94-Nam2_CI_9'= percentage)
Histoplasma_ohiense_SECH_94.Nam2_CI_9class1 <- Histoplasma_ohiense_SECH_94.Nam2_CI_9class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_94.Nam2_CI_9class_N <- Histoplasma_ohiense_SECH_94.Nam2_CI_9class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_95.Nam2_CI_10 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_95-Nam2_CI_10.RM/Histoplasma_ohiense_SECH_95-Nam2_CI_10.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_95.Nam2_CI_10) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_95.Nam2_CI_10 <- Histoplasma_ohiense_SECH_95.Nam2_CI_10 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_95.Nam2_CI_101 <- subset(Histoplasma_ohiense_SECH_95.Nam2_CI_10,class != "Low_complexity")
#Histoplasma_ohiense_SECH_95.Nam2_CI_102 <- subset(Histoplasma_ohiense_SECH_95.Nam2_CI_101,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_95.Nam2_CI_10,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_95.Nam2_CI_10.tsv")

Histoplasma_ohiense_SECH_95.Nam2_CI_10class<- Histoplasma_ohiense_SECH_95.Nam2_CI_102 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25510614)
write_tsv(Histoplasma_ohiense_SECH_95.Nam2_CI_10class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_95.Nam2_CI_10class.tsv")
Histoplasma_ohiense_SECH_95.Nam2_CI_10class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_95.Nam2_CI_10class.tsv")
Histoplasma_ohiense_SECH_95.Nam2_CI_10class <- Histoplasma_ohiense_SECH_95.Nam2_CI_10class %>% mutate('Histoplasma ohiense SECH_95.Nam2_CI_10'= percentage)
Histoplasma_ohiense_SECH_95.Nam2_CI_10class1 <- Histoplasma_ohiense_SECH_95.Nam2_CI_10class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_95.Nam2_CI_10class_N <- Histoplasma_ohiense_SECH_95.Nam2_CI_10class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_96.Nam2_CI_17 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_96-Nam2_CI_17.RM/Histoplasma_ohiense_SECH_96-Nam2_CI_17.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_96.Nam2_CI_17) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_96.Nam2_CI_17 <- Histoplasma_ohiense_SECH_96.Nam2_CI_17 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_96.Nam2_CI_171 <- subset(Histoplasma_ohiense_SECH_96.Nam2_CI_17,class != "Low_complexity")
#Histoplasma_ohiense_SECH_96.Nam2_CI_172 <- subset(Histoplasma_ohiense_SECH_96.Nam2_CI_171,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_96.Nam2_CI_17,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_96.Nam2_CI_17.tsv")

Histoplasma_ohiense_SECH_96.Nam2_CI_17class<- Histoplasma_ohiense_SECH_96.Nam2_CI_172 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25893578)
write_tsv(Histoplasma_ohiense_SECH_96.Nam2_CI_17class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_96.Nam2_CI_17class.tsv")
Histoplasma_ohiense_SECH_96.Nam2_CI_17class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_96.Nam2_CI_17class.tsv")
Histoplasma_ohiense_SECH_96.Nam2_CI_17class <- Histoplasma_ohiense_SECH_96.Nam2_CI_17class %>% mutate('Histoplasma ohiense SECH_96.Nam2_CI_17'= percentage)
Histoplasma_ohiense_SECH_96.Nam2_CI_17class1 <- Histoplasma_ohiense_SECH_96.Nam2_CI_17class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_96.Nam2_CI_17class_N <- Histoplasma_ohiense_SECH_96.Nam2_CI_17class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_97.Nam2_CI_18 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_97-Nam2_CI_18.RM/Histoplasma_ohiense_SECH_97-Nam2_CI_18.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_97.Nam2_CI_18) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_97.Nam2_CI_18 <- Histoplasma_ohiense_SECH_97.Nam2_CI_18 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_97.Nam2_CI_181 <- subset(Histoplasma_ohiense_SECH_97.Nam2_CI_18,class != "Low_complexity")
#Histoplasma_ohiense_SECH_97.Nam2_CI_182 <- subset(Histoplasma_ohiense_SECH_97.Nam2_CI_181,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_97.Nam2_CI_18,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_97-Nam2_CI_18.tsv")

Histoplasma_ohiense_SECH_97.Nam2_CI_18class<- Histoplasma_ohiense_SECH_97.Nam2_CI_182 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28259788)
write_tsv(Histoplasma_ohiense_SECH_97.Nam2_CI_18class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_97-Nam2_CI_18class.tsv")
Histoplasma_ohiense_SECH_97.Nam2_CI_18class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_97-Nam2_CI_18class.tsv")
Histoplasma_ohiense_SECH_97.Nam2_CI_18class <- Histoplasma_ohiense_SECH_97.Nam2_CI_18class %>% mutate('Histoplasma ohiense SECH_97-Nam2_CI_18'= percentage)
Histoplasma_ohiense_SECH_97.Nam2_CI_18class1 <- Histoplasma_ohiense_SECH_97.Nam2_CI_18class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_97.Nam2_CI_18class_N <- Histoplasma_ohiense_SECH_97.Nam2_CI_18class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_98.Nam2_CI_30 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_98-Nam2_CI_30.RM/Histoplasma_ohiense_SECH_98-Nam2_CI_30.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_98.Nam2_CI_30) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_98.Nam2_CI_30 <- Histoplasma_ohiense_SECH_98.Nam2_CI_30 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_98.Nam2_CI_301 <- subset(Histoplasma_ohiense_SECH_98.Nam2_CI_30,class != "Low_complexity")
#Histoplasma_ohiense_SECH_98.Nam2_CI_302 <- subset(Histoplasma_ohiense_SECH_98.Nam2_CI_301,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_98.Nam2_CI_30,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_98.Nam2_CI_30.tsv")

Histoplasma_ohiense_SECH_98.Nam2_CI_30class<- Histoplasma_ohiense_SECH_98.Nam2_CI_302 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24259982)
write_tsv(Histoplasma_ohiense_SECH_98.Nam2_CI_30class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_98.Nam2_CI_30class.tsv")
Histoplasma_ohiense_SECH_98.Nam2_CI_30class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_98.Nam2_CI_30class.tsv")
Histoplasma_ohiense_SECH_98.Nam2_CI_30class <- Histoplasma_ohiense_SECH_98.Nam2_CI_30class %>% mutate('Histoplasma ohiense SECH_98.Nam2_CI_30'= percentage)
Histoplasma_ohiense_SECH_98.Nam2_CI_30class1 <- Histoplasma_ohiense_SECH_98.Nam2_CI_30class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_98.Nam2_CI_30class_N <- Histoplasma_ohiense_SECH_98.Nam2_CI_30class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_99.Nam2_CI_35 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_99-Nam2_CI_35.RM/Histoplasma_ohiense_SECH_99-Nam2_CI_35.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_99.Nam2_CI_35) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_99.Nam2_CI_35 <- Histoplasma_ohiense_SECH_99.Nam2_CI_35 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiense_SECH_99.Nam2_CI_351 <- subset(Histoplasma_ohiense_SECH_99.Nam2_CI_35,class != "Low_complexity")
#Histoplasma_ohiense_SECH_99.Nam2_CI_352 <- subset(Histoplasma_ohiense_SECH_99.Nam2_CI_351,class != "Simple_repeat")
write_tsv(Histoplasma_ohiense_SECH_99.Nam2_CI_35,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_99.Nam2_CI_35.tsv")

Histoplasma_ohiense_SECH_99.Nam2_CI_35class<- Histoplasma_ohiense_SECH_99.Nam2_CI_352 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29591009)
write_tsv(Histoplasma_ohiense_SECH_99.Nam2_CI_35class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiense_SECH_99.Nam2_CI_35class.tsv")
Histoplasma_ohiense_SECH_99.Nam2_CI_35class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_99.Nam2_CI_35class.tsv")
Histoplasma_ohiense_SECH_99.Nam2_CI_35class <- Histoplasma_ohiense_SECH_99.Nam2_CI_35class %>% mutate('Histoplasma ohiense SECH_99.Nam2_CI_35'= percentage)
Histoplasma_ohiense_SECH_99.Nam2_CI_35class1 <- Histoplasma_ohiense_SECH_99.Nam2_CI_35class %>% select(-c(n,total,percentage))
Histoplasma_ohiense_SECH_99.Nam2_CI_35class_N <- Histoplasma_ohiense_SECH_99.Nam2_CI_35class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiensis_Hc1986 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiensis_Hc1986.RM/Histoplasma_ohiensis_Hc1986.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiensis_Hc1986) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiensis_Hc1986 <- Histoplasma_ohiensis_Hc1986 %>% mutate(length=endq-beginq +1) 
#Histoplasma_ohiensis_Hc19861 <- subset(Histoplasma_ohiensis_Hc1986,class != "Low_complexity")
#Histoplasma_ohiensis_Hc19862 <- subset(Histoplasma_ohiensis_Hc19861,class != "Simple_repeat")
write_tsv(Histoplasma_ohiensis_Hc1986,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiensis_Hc1986.tsv")

Histoplasma_ohiensis_Hc1986class<- Histoplasma_ohiensis_Hc19862 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/13148448)
write_tsv(Histoplasma_ohiensis_Hc1986class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table/Histoplasma_ohiensis_Hc1986class.tsv")
Histoplasma_ohiensis_Hc1986class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiensis_Hc1986class.tsv")
Histoplasma_ohiensis_Hc1986class <- Histoplasma_ohiensis_Hc1986class %>% mutate('Histoplasma ohiensis Hc1986'= percentage)
Histoplasma_ohiensis_Hc1986class1 <- Histoplasma_ohiensis_Hc1986class %>% select(-c(n,total,percentage))
Histoplasma_ohiensis_Hc1986class_N <- Histoplasma_ohiensis_Hc1986class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_19VMG15 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_19VMG-15.RM/Histoplasma_capsulatum_19VMG-15.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_19VMG15) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_19VMG15 <- Histoplasma_capsulatum_19VMG15 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_19VMG15,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_19VMG15.tsv")

Histoplasma_capsulatum_19VMG15class<- Histoplasma_capsulatum_19VMG15 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26953388)
write_tsv(Histoplasma_capsulatum_19VMG15class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_19VMG15class.tsv")
Histoplasma_capsulatum_19VMG15class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_19VMG15class.tsv")
Histoplasma_capsulatum_19VMG15class <- Histoplasma_capsulatum_19VMG15class %>% mutate('Histoplasma capsulatum 19VMG15'= percentage)
Histoplasma_capsulatum_19VMG15class1 <- Histoplasma_capsulatum_19VMG15class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_19VMG15class_N <- Histoplasma_capsulatum_19VMG15class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BB <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-CM7256-xx-CL-SEN-xxxx-036-BB.RM/Histoplasma_capsulatum_HISSP-CM7256-xx-CL-SEN-xxxx-036-BB.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BB) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BB <- Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BB %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BB,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM7256-xx-CL-SEN-xxxx-036-BB.tsv")

Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass<- Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BB %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/59797550)
write_tsv(Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM7256-xx-CL-SEN-xxxx-036-BBclass.tsv")
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM7256-xx-CL-SEN-xxxx-036-BBclass.tsv")
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass <- Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass %>% mutate('Histoplasma capsulatum HISSP-CM7256-xx-CL-SEN-xxxx-036-BB'= percentage)
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass <- Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass %>% select (-c('Histoplasma capsulatum HISSP-CM7256-xx-CL-SEN-xxxx-036-BB'))
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass1 <- Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass_N <- Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_Histo485P20 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_Histo-485P20.RM/Histoplasma_capsulatum_Histo-485P20.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_Histo485P20) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_Histo485P20 <- Histoplasma_capsulatum_Histo485P20 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_Histo485P20,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Histo-485P20.tsv")

Histoplasma_capsulatum_Histo485P20class<- Histoplasma_capsulatum_Histo485P20 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28091669)
write_tsv(Histoplasma_capsulatum_Histo485P20class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Histo-485P20class.tsv")
Histoplasma_capsulatum_Histo485P20class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Histo-485P20class.tsv")
Histoplasma_capsulatum_Histo485P20class <- Histoplasma_capsulatum_Histo485P20class %>% mutate('Histoplasma capsulatum Histo-485P20'= percentage)
Histoplasma_capsulatum_Histo485P20class1 <- Histoplasma_capsulatum_Histo485P20class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_Histo485P20class_N <- Histoplasma_capsulatum_Histo485P20class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_JB_01752Hc_01752 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_JB_01752-Hc_01752.RM/Histoplasma_capsulatum_JB_01752-Hc_01752.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_JB_01752Hc_01752) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_JB_01752Hc_01752 <- Histoplasma_capsulatum_JB_01752Hc_01752 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_JB_01752Hc_01752,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_01752-Hc_01752.tsv")

Histoplasma_capsulatum_JB_01752Hc_01752class<- Histoplasma_capsulatum_JB_01752Hc_01752 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23508609)
write_tsv(Histoplasma_capsulatum_JB_01752Hc_01752class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_01752-Hc_01752class.tsv")
Histoplasma_capsulatum_JB_01752Hc_01752class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_01752-Hc_01752class.tsv")
Histoplasma_capsulatum_JB_01752Hc_01752class <- Histoplasma_capsulatum_JB_01752Hc_01752class %>% mutate('Histoplasma capsulatum JB_01752-Hc_01752'= percentage)
Histoplasma_capsulatum_JB_01752Hc_01752class <- Histoplasma_capsulatum_JB_01752Hc_01752class %>% select (-c('Histoplasma capsulatum JB_01752-Hc_01752'))
Histoplasma_capsulatum_JB_01752Hc_01752class1 <- Histoplasma_capsulatum_JB_01752Hc_01752class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_JB_01752Hc_01752class_N <- Histoplasma_capsulatum_JB_01752Hc_01752class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_JB_021091Hc_021091 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_JB_021091-Hc_021091.RM/Histoplasma_capsulatum_JB_021091-Hc_021091.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_JB_021091Hc_021091) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_JB_021091Hc_021091 <- Histoplasma_capsulatum_JB_021091Hc_021091 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_JB_021091Hc_021091,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_021091-Hc_021091.tsv")

Histoplasma_capsulatum_JB_021091Hc_021091class<- Histoplasma_capsulatum_JB_021091Hc_021091 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22583691)
write_tsv(Histoplasma_capsulatum_JB_021091Hc_021091class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_021091-Hc_021091class.tsv")
Histoplasma_capsulatum_JB_021091Hc_021091class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_021091-Hc_021091class.tsv")
Histoplasma_capsulatum_JB_021091Hc_021091class <- Histoplasma_capsulatum_JB_021091Hc_021091class %>% mutate('Histoplasma capsulatum JB_021091-Hc_021091'= percentage)
Histoplasma_capsulatum_JB_021091Hc_021091class1 <- Histoplasma_capsulatum_JB_021091Hc_021091class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_JB_021091Hc_021091class_N <- Histoplasma_capsulatum_JB_021091Hc_021091class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_JB_031837Hc_031837 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_JB_031837-Hc_031837.RM/Histoplasma_capsulatum_JB_031837-Hc_031837.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_JB_031837Hc_031837) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_JB_031837Hc_031837 <- Histoplasma_capsulatum_JB_031837Hc_031837 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_JB_031837Hc_031837,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_031837-Hc_031837.tsv")

Histoplasma_capsulatum_JB_031837Hc_031837class<- Histoplasma_capsulatum_JB_031837Hc_031837 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/27508086)
write_tsv(Histoplasma_capsulatum_JB_031837Hc_031837class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_031837-Hc_031837class.tsv")
Histoplasma_capsulatum_JB_031837Hc_031837class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_031837-Hc_031837class.tsv")
Histoplasma_capsulatum_JB_031837Hc_031837class <- Histoplasma_capsulatum_JB_031837Hc_031837class %>% mutate('Histoplasma capsulatum JB_031837-Hc_031837'= percentage)
Histoplasma_capsulatum_JB_031837Hc_031837class1 <- Histoplasma_capsulatum_JB_031837Hc_031837class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_JB_031837Hc_031837class_N <- Histoplasma_capsulatum_JB_031837Hc_031837class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_JB_042430Hc_042430 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_JB_042430-Hc_042430.RM/Histoplasma_capsulatum_JB_042430-Hc_042430.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_JB_042430Hc_042430) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_JB_042430Hc_042430 <- Histoplasma_capsulatum_JB_042430Hc_042430 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_JB_042430Hc_042430,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_042430-Hc_042430.tsv")

Histoplasma_capsulatum_JB_042430Hc_042430class<- Histoplasma_capsulatum_JB_042430Hc_042430 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26691293)
write_tsv(Histoplasma_capsulatum_JB_042430Hc_042430class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_042430-Hc_042430class.tsv")
Histoplasma_capsulatum_JB_042430Hc_042430class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_042430-Hc_042430class.tsv")
Histoplasma_capsulatum_JB_042430Hc_042430class <- Histoplasma_capsulatum_JB_042430Hc_042430class %>% mutate('Histoplasma capsulatum JB_042430-Hc_042430'= percentage)
Histoplasma_capsulatum_JB_042430Hc_042430class1 <- Histoplasma_capsulatum_JB_042430Hc_042430class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_JB_042430Hc_042430class_N <- Histoplasma_capsulatum_JB_042430Hc_042430class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_JB_062632Hc_062632 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_JB_062632-Hc_062632.RM/Histoplasma_capsulatum_JB_062632-Hc_062632.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_JB_062632Hc_062632) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_JB_062632Hc_062632 <- Histoplasma_capsulatum_JB_062632Hc_062632 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_JB_062632Hc_062632,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_06263-2Hc_062632.tsv")

Histoplasma_capsulatum_JB_062632Hc_062632class<- Histoplasma_capsulatum_JB_062632Hc_062632 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25946185)
write_tsv(Histoplasma_capsulatum_JB_062632Hc_062632class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_062632-Hc_062632class.tsv")
Histoplasma_capsulatum_JB_062632Hc_062632class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_062632-Hc_062632class.tsv")
Histoplasma_capsulatum_JB_062632Hc_062632class <- Histoplasma_capsulatum_JB_062632Hc_062632class %>% mutate('Histoplasma capsulatum JB_062632-Hc_062632'= percentage)
Histoplasma_capsulatum_JB_062632Hc_062632class1 <- Histoplasma_capsulatum_JB_062632Hc_062632class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_JB_062632Hc_062632class_N <- Histoplasma_capsulatum_JB_062632Hc_062632class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_JB_062775Hc_062775 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_JB_062775-Hc_062775.RM/Histoplasma_capsulatum_JB_062775-Hc_062775.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_JB_062775Hc_062775) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_JB_062775Hc_062775 <- Histoplasma_capsulatum_JB_062775Hc_062775 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_JB_062775Hc_062775,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_062775-Hc_062775.tsv")

Histoplasma_capsulatum_JB_062775Hc_062775class<- Histoplasma_capsulatum_JB_062775Hc_062775 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25015646)
write_tsv(Histoplasma_capsulatum_JB_062775Hc_062775class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_062775-Hc_062775class.tsv")
Histoplasma_capsulatum_JB_062775Hc_062775class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_062775-Hc_062775class.tsv")
Histoplasma_capsulatum_JB_062775Hc_062775class <- Histoplasma_capsulatum_JB_062775Hc_062775class %>% mutate('Histoplasma capsulatum JB_062775-Hc_062775'= percentage)
#Histoplasma_capsulatum_JB_062775Hc_062775class <- Histoplasma_capsulatum_JB_062775Hc_062775class %>% select (-c('Histoplasma capsulatum JB_062775-Hc_062775'))
Histoplasma_capsulatum_JB_062775Hc_062775class1 <- Histoplasma_capsulatum_JB_062775Hc_062775class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_JB_062775Hc_062775class_N <- Histoplasma_capsulatum_JB_062775Hc_062775class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_JB_073129Hc_073129 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_JB_073129-Hc_073129.RM/Histoplasma_capsulatum_JB_073129-Hc_073129.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_JB_073129Hc_073129) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_JB_073129Hc_073129 <- Histoplasma_capsulatum_JB_073129Hc_073129 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_JB_073129Hc_073129,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_073129-Hc_073129.tsv")

Histoplasma_capsulatum_JB_073129Hc_073129class<- Histoplasma_capsulatum_JB_073129Hc_073129 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23705142)
write_tsv(Histoplasma_capsulatum_JB_073129Hc_073129class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_073129-Hc_073129class.tsv")
Histoplasma_capsulatum_JB_073129Hc_073129class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_073129-Hc_073129class.tsv")
Histoplasma_capsulatum_JB_073129Hc_073129class <- Histoplasma_capsulatum_JB_073129Hc_073129class %>% mutate('Histoplasma capsulatum JB_073129-Hc_073129'= n)
Histoplasma_capsulatum_JB_073129Hc_073129class1 <- Histoplasma_capsulatum_JB_073129Hc_073129class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_JB_073129Hc_073129class_N <- Histoplasma_capsulatum_JB_073129Hc_073129class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_JB_083285_2Hc_083285_2 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_JB_083285_2-Hc_083285_2.RM/Histoplasma_capsulatum_JB_083285_2-Hc_083285_2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_JB_083285_2Hc_083285_2) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_JB_083285_2Hc_083285_2 <- Histoplasma_capsulatum_JB_083285_2Hc_083285_2 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_JB_083285_2Hc_083285_2,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_083285_2-Hc_083285_2.tsv")

Histoplasma_capsulatum_JB_083285_2Hc_083285_2class<- Histoplasma_capsulatum_JB_083285_2Hc_083285_2 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/27324225)
write_tsv(Histoplasma_capsulatum_JB_083285_2Hc_083285_2class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_083285_2-Hc_083285_2class.tsv")
Histoplasma_capsulatum_JB_083285_2Hc_083285_2class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_083285_2-Hc_083285_2class.tsv")
Histoplasma_capsulatum_JB_083285_2Hc_083285_2class <- Histoplasma_capsulatum_JB_083285_2Hc_083285_2class %>% mutate('Histoplasma capsulatum JB_083285_2-Hc_083285_2'= percentage)
Histoplasma_capsulatum_JB_083285_2Hc_083285_2class1 <- Histoplasma_capsulatum_JB_083285_2Hc_083285_2class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_JB_083285_2Hc_083285_2class_N <- Histoplasma_capsulatum_JB_083285_2Hc_083285_2class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_SECH_101Nam2_G184A <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SECH_101-Nam2_G184A.RM/Histoplasma_capsulatum_SECH_101-Nam2_G184A.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_SECH_101Nam2_G184A) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_SECH_101Nam2_G184A <- Histoplasma_capsulatum_SECH_101Nam2_G184A %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_SECH_101Nam2_G184A,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_101-Nam2_G184A.tsv")

Histoplasma_capsulatum_SECH_101Nam2_G184Aclass<- Histoplasma_capsulatum_SECH_101Nam2_G184A %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23224387)
write_tsv(Histoplasma_capsulatum_SECH_101Nam2_G184Aclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_101-Nam2_G184Aclass.tsv")
Histoplasma_capsulatum_SECH_101Nam2_G184Aclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_101-Nam2_G184Aclass.tsv")
Histoplasma_capsulatum_SECH_101Nam2_G184Aclass <- Histoplasma_capsulatum_SECH_101Nam2_G184Aclass %>% mutate('Histoplasma capsulatum SECH_101-Nam2_G184A'= percentage)
Histoplasma_capsulatum_SECH_101Nam2_G184Aclass1 <- Histoplasma_capsulatum_SECH_101Nam2_G184Aclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_SECH_101Nam2_G184Aclass_N <- Histoplasma_capsulatum_SECH_101Nam2_G184Aclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiB <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SECH_107-mis_Hc_duboisii-B.RM/Histoplasma_capsulatum_SECH_107-mis_Hc_duboisii-B.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiB) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiB <- Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiB %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiB,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_107-mis_Hc_duboisii-B.tsv")

Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass<- Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiB %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22891084)
write_tsv(Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_107-mis_Hc_duboisii-Bclass.tsv")
Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_107-mis_Hc_duboisii-Bclass.tsv")
Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass <- Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass %>% mutate('Histoplasma capsulatum SECH_107-mis_Hc_duboisii-B'= percentage)
Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass1 <- Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass_N <- Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_SECH_109 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SECH_109.RM/Histoplasma_capsulatum_SECH_109.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_SECH_109) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_SECH_109 <- Histoplasma_capsulatum_SECH_109 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_SECH_109,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_109.tsv")

Histoplasma_capsulatum_SECH_109class<- Histoplasma_capsulatum_SECH_109 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23941220)
write_tsv(Histoplasma_capsulatum_SECH_109class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_109class.tsv")
Histoplasma_capsulatum_SECH_109class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_109class.tsv")
Histoplasma_capsulatum_SECH_109class <- Histoplasma_capsulatum_SECH_109class %>% mutate('Histoplasma capsulatum SECH_109'= percentage)
Histoplasma_capsulatum_SECH_109class <- Histoplasma_capsulatum_SECH_109class %>% select (-c('Histoplasma capsulatum SECH_109'))
Histoplasma_capsulatum_SECH_109class1 <- Histoplasma_capsulatum_SECH_109class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_SECH_109class_N <- Histoplasma_capsulatum_SECH_109class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_SECH_110 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SECH_110.RM/Histoplasma_capsulatum_SECH_110.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_SECH_110) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_SECH_110 <- Histoplasma_capsulatum_SECH_110 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_SECH_110,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_110.tsv")

Histoplasma_capsulatum_SECH_110class<- Histoplasma_capsulatum_SECH_110 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22122765)
write_tsv(Histoplasma_capsulatum_SECH_110class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_110class.tsv")
Histoplasma_capsulatum_SECH_110class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_110class.tsv")
Histoplasma_capsulatum_SECH_110class <- Histoplasma_capsulatum_SECH_110class %>% mutate('Histoplasma capsulatum SECH_110'= percentage)
Histoplasma_capsulatum_SECH_110class <- Histoplasma_capsulatum_SECH_110class %>% select (-c('Histoplasma capsulatum SECH_110'))
Histoplasma_capsulatum_SECH_110class1 <- Histoplasma_capsulatum_SECH_110class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_SECH_110class_N <- Histoplasma_capsulatum_SECH_110class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HCH143 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HCH143.RM/Histoplasma_capsulatum_HCH143.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HCH143) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HCH143 <- Histoplasma_capsulatum_HCH143 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HCH143,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HCH143.tsv")

Histoplasma_capsulatum_HCH143class<- Histoplasma_capsulatum_HCH143 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29042875)
write_tsv(Histoplasma_capsulatum_HCH143class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HCH143class.tsv")
Histoplasma_capsulatum_HCH143class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HCH143class.tsv")
Histoplasma_capsulatum_HCH143class <- Histoplasma_capsulatum_HCH143class %>% mutate('Histoplasma capsulatum HCH143'= percentage)
Histoplasma_capsulatum_HCH143class1 <- Histoplasma_capsulatum_HCH143class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HCH143class_N <- Histoplasma_capsulatum_HCH143class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HCG186A <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HCG186A.RM/Histoplasma_capsulatum_HCG186A.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HCG186A) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HCG186A <- Histoplasma_capsulatum_HCG186A %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HCG186A,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HCG186A.tsv")

Histoplasma_capsulatum_HCG186Aclass<- Histoplasma_capsulatum_HCG186A %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/18558703)
write_tsv(Histoplasma_capsulatum_HCG186Aclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HCG186Aclass.tsv")
Histoplasma_capsulatum_HCG186Aclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HCG186Aclass.tsv")
Histoplasma_capsulatum_HCG186Aclass <- Histoplasma_capsulatum_HCG186Aclass %>% mutate('Histoplasma capsulatum HCG186A'= percentage)
Histoplasma_capsulatum_HCG186Aclass1 <- Histoplasma_capsulatum_HCG186Aclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HCG186Aclass_N <- Histoplasma_capsulatum_HCG186Aclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2.RM/Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2 <- Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2.tsv")

Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class<- Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22694737)
write_tsv(Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class.tsv")
Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class.tsv")
Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class <- Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class %>% mutate('Histoplasma capsulatum NACVFR_Histo_HC1070058_2'= percentage)
Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class1 <- Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class_N <- Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGTRO0285 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGTRO0285.RM/Histoplasma_capsulatum_HISSP-FGTRO0285.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGTRO0285) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGTRO0285 <- Histoplasma_capsulatum_HISSPFGTRO0285 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGTRO0285,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGTRO0285.tsv")

Histoplasma_capsulatum_HISSPFGTRO0285class<- Histoplasma_capsulatum_HISSPFGTRO0285 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22088586)
write_tsv(Histoplasma_capsulatum_HISSPFGTRO0285class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGTRO0285class.tsv")
Histoplasma_capsulatum_HISSPFGTRO0285class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGTRO0285class.tsv")
Histoplasma_capsulatum_HISSPFGTRO0285class <- Histoplasma_capsulatum_HISSPFGTRO0285class %>% mutate('Histoplasma capsulatum HISSP-FGTRO0285'= n)
Histoplasma_capsulatum_HISSPFGTRO0285class1 <- Histoplasma_capsulatum_HISSPFGTRO0285class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGTRO0285class_N <- Histoplasma_capsulatum_HISSPFGTRO0285class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGPSO2043 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGPSO2043.RM/Histoplasma_capsulatum_HISSP-FGPSO2043.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGPSO2043) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGPSO2043 <- Histoplasma_capsulatum_HISSPFGPSO2043 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGPSO2043,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPSO2043.tsv")

Histoplasma_capsulatum_HISSPFGPSO2043class<- Histoplasma_capsulatum_HISSPFGPSO2043 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/19677294)
write_tsv(Histoplasma_capsulatum_HISSPFGPSO2043class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPSO2043class.tsv")
Histoplasma_capsulatum_HISSPFGPSO2043class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPSO2043class.tsv")
Histoplasma_capsulatum_HISSPFGPSO2043class <- Histoplasma_capsulatum_HISSPFGPSO2043class %>% mutate('Histoplasma capsulatum HISSP-FGPSO2043'= percentage)
Histoplasma_capsulatum_HISSPFGPSO2043class1 <- Histoplasma_capsulatum_HISSPFGPSO2043class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGPSO2043class_N <- Histoplasma_capsulatum_HISSPFGPSO2043class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGPIE2055 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGPIE2055.RM/Histoplasma_capsulatum_HISSP-FGPIE2055.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGPIE2055) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGPIE2055 <- Histoplasma_capsulatum_HISSPFGPIE2055 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGPIE2055,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPIE2055.tsv")

Histoplasma_capsulatum_HISSPFGPIE2055class<- Histoplasma_capsulatum_HISSPFGPIE2055 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28918902)
write_tsv(Histoplasma_capsulatum_HISSPFGPIE2055class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPIE2055class.tsv")
Histoplasma_capsulatum_HISSPFGPIE2055class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPIE2055class.tsv")
Histoplasma_capsulatum_HISSPFGPIE2055class <- Histoplasma_capsulatum_HISSPFGPIE2055class %>% mutate('Histoplasma capsulatum HISSP-FGPIE2055'= percentage)
Histoplasma_capsulatum_HISSPFGPIE2055class1 <- Histoplasma_capsulatum_HISSPFGPIE2055class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGPIE2055class_N <- Histoplasma_capsulatum_HISSPFGPIE2055class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGPIA2052 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGPIA2052.RM/Histoplasma_capsulatum_HISSP-FGPIA2052.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGPIA2052) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGPIA2052 <- Histoplasma_capsulatum_HISSPFGPIA2052 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGPIA2052,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPIA2052.tsv")

Histoplasma_capsulatum_HISSPFGPIA2052class<- Histoplasma_capsulatum_HISSPFGPIA2052 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/20859097)
write_tsv(Histoplasma_capsulatum_HISSPFGPIA2052class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPIA2052class.tsv")
Histoplasma_capsulatum_HISSPFGPIA2052class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPIA2052class.tsv")
Histoplasma_capsulatum_HISSPFGPIA2052class <- Histoplasma_capsulatum_HISSPFGPIA2052class %>% mutate('Histoplasma capsulatum HISSP-FGPIA2052'= percentage)
Histoplasma_capsulatum_HISSPFGPIA2052class1 <- Histoplasma_capsulatum_HISSPFGPIA2052class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGPIA2052class_N <- Histoplasma_capsulatum_HISSPFGPIA2052class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGPERS2034 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGPERS2034.RM/Histoplasma_capsulatum_HISSP-FGPERS2034.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGPERS2034) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGPERS2034 <- Histoplasma_capsulatum_HISSPFGPERS2034 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGPERS2034,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPERS2034.tsv")

Histoplasma_capsulatum_HISSPFGPERS2034class<- Histoplasma_capsulatum_HISSPFGPERS2034 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24259670)
write_tsv(Histoplasma_capsulatum_HISSPFGPERS2034class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPERS2034class.tsv")
Histoplasma_capsulatum_HISSPFGPERS2034class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPERS2034class.tsv")
Histoplasma_capsulatum_HISSPFGPERS2034class <- Histoplasma_capsulatum_HISSPFGPERS2034class %>% mutate('Histoplasma capsulatum HISSP-FGPERS2034'= percentage)
Histoplasma_capsulatum_HISSPFGPERS2034class1 <- Histoplasma_capsulatum_HISSPFGPERS2034class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGPERS2034class_N <- Histoplasma_capsulatum_HISSPFGPERS2034class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGMAR2044 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGMAR2044.RM/Histoplasma_capsulatum_HISSP-FGMAR2044.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGMAR2044) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGMAR2044 <- Histoplasma_capsulatum_HISSPFGMAR2044 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGMAR2044,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGMAR2044.tsv")

Histoplasma_capsulatum_HISSPFGMAR2044class<- Histoplasma_capsulatum_HISSPFGMAR2044 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24065669)
write_tsv(Histoplasma_capsulatum_HISSPFGMAR2044class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGMAR2044class.tsv")
Histoplasma_capsulatum_HISSPFGMAR2044class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGMAR2044class.tsv")
Histoplasma_capsulatum_HISSPFGMAR2044class <- Histoplasma_capsulatum_HISSPFGMAR2044class %>% mutate('Histoplasma capsulatum HISSP-FGMAR2044'= percentage)
Histoplasma_capsulatum_HISSPFGMAR2044class1 <- Histoplasma_capsulatum_HISSPFGMAR2044class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGMAR2044class_N <- Histoplasma_capsulatum_HISSPFGMAR2044class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGLIN2055 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGLIN2055.RM/Histoplasma_capsulatum_HISSP-FGLIN2055.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGLIN2055) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGLIN2055 <- Histoplasma_capsulatum_HISSPFGLIN2055 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGLIN2055,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGLIN2055.tsv")

Histoplasma_capsulatum_HISSPFGLIN2055class<- Histoplasma_capsulatum_HISSPFGLIN2055 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/47935846)
write_tsv(Histoplasma_capsulatum_HISSPFGLIN2055class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGLIN2055class.tsv")
Histoplasma_capsulatum_HISSPFGLIN2055class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGLIN2055class.tsv")
Histoplasma_capsulatum_HISSPFGLIN2055class <- Histoplasma_capsulatum_HISSPFGLIN2055class %>% mutate('Histoplasma capsulatum HISSP-FGLIN2055'= percentage)
Histoplasma_capsulatum_HISSPFGLIN2055class1 <- Histoplasma_capsulatum_HISSPFGLIN2055class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGLIN2055class_N <- Histoplasma_capsulatum_HISSPFGLIN2055class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGJOS2044 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGJOS2044.RM/Histoplasma_capsulatum_HISSP-FGJOS2044.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGJOS2044) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGJOS2044 <- Histoplasma_capsulatum_HISSPFGJOS2044 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGJOS2044,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGJOS2044.tsv")

Histoplasma_capsulatum_HISSPFGJOS2044class<- Histoplasma_capsulatum_HISSPFGJOS2044 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22245502)
write_tsv(Histoplasma_capsulatum_HISSPFGJOS2044class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGJOS2044class.tsv")
Histoplasma_capsulatum_HISSPFGJOS2044class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGJOS2044class.tsv")
Histoplasma_capsulatum_HISSPFGJOS2044class <- Histoplasma_capsulatum_HISSPFGJOS2044class %>% mutate('Histoplasma capsulatum HISSP-FGJOS2044'= percentage)
Histoplasma_capsulatum_HISSPFGJOS2044class1 <- Histoplasma_capsulatum_HISSPFGJOS2044class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGJOS2044class_N <- Histoplasma_capsulatum_HISSPFGJOS2044class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGGRE2022 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGGRE2022.RM/Histoplasma_capsulatum_HISSP-FGGRE2022.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGGRE2022) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGGRE2022 <- Histoplasma_capsulatum_HISSPFGGRE2022 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGGRE2022,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGGRE2022.tsv")

Histoplasma_capsulatum_HISSPFGGRE2022class<- Histoplasma_capsulatum_HISSPFGGRE2022 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/20040255)
write_tsv(Histoplasma_capsulatum_HISSPFGGRE2022class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGGRE2022class.tsv")
Histoplasma_capsulatum_HISSPFGGRE2022class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGGRE2022class.tsv")
Histoplasma_capsulatum_HISSPFGGRE2022class <- Histoplasma_capsulatum_HISSPFGGRE2022class %>% mutate('Histoplasma capsulatum HISSP-FGGRE2022'= percentage)
Histoplasma_capsulatum_HISSPFGGRE2022class1 <- Histoplasma_capsulatum_HISSPFGGRE2022class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGGRE2022class_N <- Histoplasma_capsulatum_HISSPFGGRE2022class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGFIN2028 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGFIN2028.RM/Histoplasma_capsulatum_HISSP-FGFIN2028.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGFIN2028) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGFIN2028 <- Histoplasma_capsulatum_HISSPFGFIN2028 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGFIN2028,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFIN2028.tsv")

Histoplasma_capsulatum_HISSPFGFIN2028class<- Histoplasma_capsulatum_HISSPFGFIN2028 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/21138080)
write_tsv(Histoplasma_capsulatum_HISSPFGFIN2028class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFIN2028class.tsv")
Histoplasma_capsulatum_HISSPFGFIN2028class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFIN2028class.tsv")
Histoplasma_capsulatum_HISSPFGFIN2028class <- Histoplasma_capsulatum_HISSPFGFIN2028class %>% mutate('Histoplasma capsulatum HISSP-FGFIN2028'= percentage)
Histoplasma_capsulatum_HISSPFGFIN2028class1 <- Histoplasma_capsulatum_HISSPFGFIN2028class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGFIN2028class_N <- Histoplasma_capsulatum_HISSPFGFIN2028class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGFER2036 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGFER2036.RM/Histoplasma_capsulatum_HISSP-FGFER2036.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGFER2036) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGFER2036 <- Histoplasma_capsulatum_HISSPFGFER2036 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGFER2036,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFER2036.tsv")

Histoplasma_capsulatum_HISSPFGFER2036class<- Histoplasma_capsulatum_HISSPFGFER2036 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22918657)
write_tsv(Histoplasma_capsulatum_HISSPFGFER2036class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFER2036class.tsv")
Histoplasma_capsulatum_HISSPFGFER2036class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFER2036class.tsv")
Histoplasma_capsulatum_HISSPFGFER2036class <- Histoplasma_capsulatum_HISSPFGFER2036class %>% mutate('Histoplasma capsulatum HHISSP-FGFER2036'= percentage)
Histoplasma_capsulatum_HISSPFGFER2036class1 <- Histoplasma_capsulatum_HISSPFGFER2036class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGFER2036class_N <- Histoplasma_capsulatum_HISSPFGFER2036class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGFAR0189 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGFAR0189.RM/Histoplasma_capsulatum_HISSP-FGFAR0189.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGFAR0189) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGFAR0189 <- Histoplasma_capsulatum_HISSPFGFAR0189 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGFAR0189,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFAR0189.tsv")

Histoplasma_capsulatum_HISSPFGFAR0189class<- Histoplasma_capsulatum_HISSPFGFAR0189 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22815282)
write_tsv(Histoplasma_capsulatum_HISSPFGFAR0189class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFAR0189class.tsv")
Histoplasma_capsulatum_HISSPFGFAR0189class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFAR0189class.tsv")
Histoplasma_capsulatum_HISSPFGFAR0189class <- Histoplasma_capsulatum_HISSPFGFAR0189class %>% mutate('Histoplasma capsulatum HISSP-FGFAR0189'= percentage)
Histoplasma_capsulatum_HISSPFGFAR0189class1 <- Histoplasma_capsulatum_HISSPFGFAR0189class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGFAR0189class_N <- Histoplasma_capsulatum_HISSPFGFAR0189class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGFAN2059 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGFAN2059.RM/Histoplasma_capsulatum_HISSP-FGFAN2059.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGFAN2059) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGFAN2059 <- Histoplasma_capsulatum_HISSPFGFAN2059 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGFAN2059,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFAN2059.tsv")

Histoplasma_capsulatum_HISSPFGFAN2059class<- Histoplasma_capsulatum_HISSPFGFAN2059 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23710091)
write_tsv(Histoplasma_capsulatum_HISSPFGFAN2059class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFAN2059class.tsv")
Histoplasma_capsulatum_HISSPFGFAN2059class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFAN2059class.tsv")
Histoplasma_capsulatum_HISSPFGFAN2059class <- Histoplasma_capsulatum_HISSPFGFAN2059class %>% mutate('Histoplasma capsulatum HISSP-FGFAN2059'= percentage)
Histoplasma_capsulatum_HISSPFGFAN2059class1 <- Histoplasma_capsulatum_HISSPFGFAN2059class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGFAN2059class_N <- Histoplasma_capsulatum_HISSPFGFAN2059class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGBON2001 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGBON2001.RM/Histoplasma_capsulatum_HISSP-FGBON2001.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGBON2001) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGBON2001 <- Histoplasma_capsulatum_HISSPFGBON2001 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGBON2001,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGBON2001.tsv")

Histoplasma_capsulatum_HISSPFGBON2001class<- Histoplasma_capsulatum_HISSPFGBON2001 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/30753906)
write_tsv(Histoplasma_capsulatum_HISSPFGBON2001class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGBON2001class.tsv")
Histoplasma_capsulatum_HISSPFGBON2001class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGBON2001class.tsv")
Histoplasma_capsulatum_HISSPFGBON2001class <- Histoplasma_capsulatum_HISSPFGBON2001class %>% mutate('Histoplasma capsulatum HISSP-FGBON2001'= percentage)
Histoplasma_capsulatum_HISSPFGBON2001class1 <- Histoplasma_capsulatum_HISSPFGBON2001class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGBON2001class_N <- Histoplasma_capsulatum_HISSPFGBON2001class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGBIK2051 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGBIK2051.RM/Histoplasma_capsulatum_HISSP-FGBIK2051.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGBIK2051) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGBIK2051 <- Histoplasma_capsulatum_HISSPFGBIK2051 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGBIK2051,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGBIK2051.tsv")

Histoplasma_capsulatum_HISSPFGBIK2051class<- Histoplasma_capsulatum_HISSPFGBIK2051 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23534554)
write_tsv(Histoplasma_capsulatum_HISSPFGBIK2051class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGBIK2051class.tsv")
Histoplasma_capsulatum_HISSPFGBIK2051class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGBIK2051class.tsv")
Histoplasma_capsulatum_HISSPFGBIK2051class <- Histoplasma_capsulatum_HISSPFGBIK2051class %>% mutate('Histoplasma capsulatum HISSP-FGBIK2051'= percentage)
Histoplasma_capsulatum_HISSPFGBIK2051class1 <- Histoplasma_capsulatum_HISSPFGBIK2051class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGBIK2051class_N <- Histoplasma_capsulatum_HISSPFGBIK2051class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPFGAMA2041 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-FGAMA2041.RM/Histoplasma_capsulatum_HISSP-FGAMA2041.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPFGAMA2041) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPFGAMA2041 <- Histoplasma_capsulatum_HISSPFGAMA2041 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPFGAMA2041,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGAMA2041.tsv")

Histoplasma_capsulatum_HISSPFGAMA2041class<- Histoplasma_capsulatum_HISSPFGAMA2041 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/20307679)
write_tsv(Histoplasma_capsulatum_HISSPFGAMA2041class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGAMA2041class.tsv")
Histoplasma_capsulatum_HISSPFGAMA2041class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGAMA2041class.tsv")
Histoplasma_capsulatum_HISSPFGAMA2041class <- Histoplasma_capsulatum_HISSPFGAMA2041class %>% mutate('Histoplasma capsulatum HISSP-FGAMA2041'= percentage)
Histoplasma_capsulatum_HISSPFGAMA2041class <- Histoplasma_capsulatum_HISSPFGAMA2041class %>% select (-c('Histoplasma capsulatum HISSP-FGAMA2041'))
Histoplasma_capsulatum_HISSPFGAMA2041class1 <- Histoplasma_capsulatum_HISSPFGAMA2041class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPFGAMA2041class_N <- Histoplasma_capsulatum_HISSPFGAMA2041class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPCM6408 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-CM6408.RM/Histoplasma_capsulatum_HISSP-CM6408.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPCM6408) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPCM6408 <- Histoplasma_capsulatum_HISSPCM6408 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPCM6408,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM6408.tsv")

Histoplasma_capsulatum_HISSPCM6408class<- Histoplasma_capsulatum_HISSPCM6408 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26008567)
write_tsv(Histoplasma_capsulatum_HISSPCM6408class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM6408class.tsv")
Histoplasma_capsulatum_HISSPCM6408class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM6408class.tsv")
Histoplasma_capsulatum_HISSPCM6408class <- Histoplasma_capsulatum_HISSPCM6408class %>% mutate('Histoplasma capsulatum HISSP-CM6408'= percentage)
Histoplasma_capsulatum_HISSPCM6408class <- Histoplasma_capsulatum_HISSPCM6408class %>% select(-c('Histoplasma capsulatum HISSP-CM6408'))
Histoplasma_capsulatum_HISSPCM6408class1 <- Histoplasma_capsulatum_HISSPCM6408class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPCM6408class_N <- Histoplasma_capsulatum_HISSPCM6408class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPCM6015 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-CM6015.RM/Histoplasma_capsulatum_HISSP-CM6015.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPCM6015) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPCM6015 <- Histoplasma_capsulatum_HISSPCM6015 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPCM6015,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM6015.tsv")

Histoplasma_capsulatum_HISSPCM6015class<- Histoplasma_capsulatum_HISSPCM6015 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/19585677)
write_tsv(Histoplasma_capsulatum_HISSPCM6015class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM6015class.tsv")
Histoplasma_capsulatum_HISSPCM6015class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM6015class.tsv")
Histoplasma_capsulatum_HISSPCM6015class <- Histoplasma_capsulatum_HISSPCM6015class %>% mutate('Histoplasma capsulatum HISSP-CM6015'= percentage)
Histoplasma_capsulatum_HISSPCM6015class1 <- Histoplasma_capsulatum_HISSPCM6015class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPCM6015class_N <- Histoplasma_capsulatum_HISSPCM6015class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSPB05821 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-B05821.RM/Histoplasma_capsulatum_HISSP-B05821.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSPB05821) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSPB05821 <- Histoplasma_capsulatum_HISSPB05821 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSPB05821,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-B05821.tsv")

Histoplasma_capsulatum_HISSPB05821class<- Histoplasma_capsulatum_HISSPB05821 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/36858421)
write_tsv(Histoplasma_capsulatum_HISSPB05821class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-B05821class.tsv")
Histoplasma_capsulatum_HISSPB05821class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-B05821class.tsv")
Histoplasma_capsulatum_HISSPB05821class <- Histoplasma_capsulatum_HISSPB05821class %>% mutate('Histoplasma capsulatum HISSP-B05821'= percentage)
Histoplasma_capsulatum_HISSPB05821class <- Histoplasma_capsulatum_HISSPB05821class %>% select(-c('Histoplasma capsulatum HISSP-B05821'))
Histoplasma_capsulatum_HISSPB05821class1 <- Histoplasma_capsulatum_HISSPB05821class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSPB05821class_N <- Histoplasma_capsulatum_HISSPB05821class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSP11571Belem1 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-11571-Belem1.RM/Histoplasma_capsulatum_HISSP-11571-Belem1.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSP11571Belem1) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSP11571Belem1 <- Histoplasma_capsulatum_HISSP11571Belem1 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSP11571Belem1,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-11571-Belem1.tsv")

Histoplasma_capsulatum_HISSP11571Belem1class<- Histoplasma_capsulatum_HISSP11571Belem1 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/21937402)
write_tsv(Histoplasma_capsulatum_HISSP11571Belem1class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-11571-Belem1class.tsv")
Histoplasma_capsulatum_HISSP11571Belem1class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-11571-Belem1class.tsv")
Histoplasma_capsulatum_HISSP11571Belem1class <- Histoplasma_capsulatum_HISSP11571Belem1class %>% mutate('Histoplasma capsulatum HISSP-11571-Belem1'= percentage)
Histoplasma_capsulatum_HISSP11571Belem1class <- Histoplasma_capsulatum_HISSP11571Belem1class %>% select(-c('Histoplasma capsulatum HISSP-11571-Belem1'))
Histoplasma_capsulatum_HISSP11571Belem1class1 <- Histoplasma_capsulatum_HISSP11571Belem1class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSP11571Belem1class_N <- Histoplasma_capsulatum_HISSP11571Belem1class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_HISSP1014Belem3 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HISSP-1014-Belem3.RM/Histoplasma_capsulatum_HISSP-1014-Belem3.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_HISSP1014Belem3) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_HISSP1014Belem3 <- Histoplasma_capsulatum_HISSP1014Belem3 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_HISSP1014Belem3,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-1014-Belem3.tsv")

Histoplasma_capsulatum_HISSP1014Belem3class<- Histoplasma_capsulatum_HISSP1014Belem3 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/31744825)
write_tsv(Histoplasma_capsulatum_HISSP1014Belem3class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-1014-Belem3class.tsv")
Histoplasma_capsulatum_HISSP1014Belem3class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-1014-Belem3class.tsv")
Histoplasma_capsulatum_HISSP1014Belem3class <- Histoplasma_capsulatum_HISSP1014Belem3class %>% mutate('Histoplasma capsulatum HISSP-1014-Belem3'= percentage)
#Histoplasma_capsulatum_HISSP1014Belem3class <- Histoplasma_capsulatum_HISSP1014Belem3class %>% select(-c('Histoplasma capsulatum HISSP-1014-Belem3'))
Histoplasma_capsulatum_HISSP1014Belem3class1 <- Histoplasma_capsulatum_HISSP1014Belem3class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_HISSP1014Belem3class_N <- Histoplasma_capsulatum_HISSP1014Belem3class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_07_12RJ <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_07_12-RJ.RM/Histoplasma_capsulatum_07_12-RJ.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_07_12RJ) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_07_12RJ <- Histoplasma_capsulatum_07_12RJ %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_07_12RJ,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_07_12-RJ.tsv")

Histoplasma_capsulatum_07_12RJclass<- Histoplasma_capsulatum_07_12RJ %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/17244466)
write_tsv(Histoplasma_capsulatum_07_12RJclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_07_12-RJclass.tsv")
Histoplasma_capsulatum_07_12RJclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_07_12-RJclass.tsv")
Histoplasma_capsulatum_07_12RJclass <- Histoplasma_capsulatum_07_12RJclass %>% mutate('Histoplasma capsulatum 07_12-RJ'= percentage)
Histoplasma_capsulatum_07_12RJclass1 <- Histoplasma_capsulatum_07_12RJclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_07_12RJclass_N <- Histoplasma_capsulatum_07_12RJclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_104_p_06 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_104_p_06.RM/Histoplasma_capsulatum_104_p_06.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_104_p_06) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_104_p_06 <- Histoplasma_capsulatum_104_p_06 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_104_p_06,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_104_p_06.tsv")

Histoplasma_capsulatum_104_p_06class<- Histoplasma_capsulatum_104_p_06 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28442945)
write_tsv(Histoplasma_capsulatum_104_p_06class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_104_p_06class.tsv")
Histoplasma_capsulatum_104_p_06class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_104_p_06class.tsv")
Histoplasma_capsulatum_104_p_06class <- Histoplasma_capsulatum_104_p_06class %>% mutate('Histoplasma capsulatum 104_p_06'= percentage)
Histoplasma_capsulatum_104_p_06class1 <- Histoplasma_capsulatum_104_p_06class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_104_p_06class_N <- Histoplasma_capsulatum_104_p_06class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_104_P_19 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_104_P_19.RM/Histoplasma_capsulatum_104_P_19.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_104_P_19) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_104_P_19 <- Histoplasma_capsulatum_104_P_19 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_104_P_19,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_104_P_19.tsv")

Histoplasma_capsulatum_104_P_19class<- Histoplasma_capsulatum_104_P_19 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28325613)
write_tsv(Histoplasma_capsulatum_104_P_19class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_104_P_19class.tsv")
Histoplasma_capsulatum_104_P_19class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_104_P_19class.tsv")
Histoplasma_capsulatum_104_P_19class <- Histoplasma_capsulatum_104_P_19class %>% mutate('Histoplasma capsulatum 104_P_19'= percentage)
Histoplasma_capsulatum_104_P_19class1 <- Histoplasma_capsulatum_104_P_19class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_104_P_19class_N <- Histoplasma_capsulatum_104_P_19class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_117_p_12 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_117_p_12.RM/Histoplasma_capsulatum_117_p_12.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_117_p_12) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_117_p_12 <- Histoplasma_capsulatum_117_p_12 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_117_p_12,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_117_p_12.tsv")

Histoplasma_capsulatum_117_p_12class<- Histoplasma_capsulatum_117_p_12 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28587451)
write_tsv(Histoplasma_capsulatum_117_p_12class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_117_p_12class.tsv")
Histoplasma_capsulatum_117_p_12class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_117_p_12class.tsv")
Histoplasma_capsulatum_117_p_12class <- Histoplasma_capsulatum_117_p_12class %>% mutate('Histoplasma capsulatum 117_p_12'= percentage)
Histoplasma_capsulatum_117_p_12class1 <- Histoplasma_capsulatum_117_p_12class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_117_p_12class_N <- Histoplasma_capsulatum_117_p_12class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_122_p_10_B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_122_p_10_B.RM/Histoplasma_capsulatum_122_p_10_B.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_122_p_10_B) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_122_p_10_B <- Histoplasma_capsulatum_122_p_10_B %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_122_p_10_B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_122_p_10_B.tsv")

Histoplasma_capsulatum_122_p_10_Bclass<- Histoplasma_capsulatum_122_p_10_B %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28585056)
write_tsv(Histoplasma_capsulatum_122_p_10_Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_122_p_10_Bclass.tsv")
Histoplasma_capsulatum_122_p_10_Bclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_122_p_10_Bclass.tsv")
Histoplasma_capsulatum_122_p_10_Bclass <- Histoplasma_capsulatum_122_p_10_Bclass %>% mutate('Histoplasma capsulatum 122_p_10_B'= percentage)
Histoplasma_capsulatum_122_p_10_Bclass1 <- Histoplasma_capsulatum_122_p_10_Bclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_122_p_10_Bclass_N <- Histoplasma_capsulatum_122_p_10_Bclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_136_P_07 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_136_P_07.RM/Histoplasma_capsulatum_136_P_07.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_136_P_07) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_136_P_07 <- Histoplasma_capsulatum_136_P_07 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_136_P_07,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_136_P_07.tsv")

Histoplasma_capsulatum_136_P_07class<- Histoplasma_capsulatum_136_P_07 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28536980)
write_tsv(Histoplasma_capsulatum_136_P_07class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_136_P_07class.tsv")
Histoplasma_capsulatum_136_P_07class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_136_P_07class.tsv")
Histoplasma_capsulatum_136_P_07class <- Histoplasma_capsulatum_136_P_07class %>% mutate('Histoplasma capsulatum 136_P_07'= percentage)
Histoplasma_capsulatum_136_P_07class1 <- Histoplasma_capsulatum_136_P_07class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_136_P_07class_N <- Histoplasma_capsulatum_136_P_07class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_144_p_08 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_144_p_08.RM/Histoplasma_capsulatum_144_p_08.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_144_p_08) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_144_p_08 <- Histoplasma_capsulatum_144_p_08 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_144_p_08,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_144_p_08.tsv")

Histoplasma_capsulatum_144_p_08class<- Histoplasma_capsulatum_144_p_08 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28329803)
write_tsv(Histoplasma_capsulatum_144_p_08class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_144_p_08class.tsv")
Histoplasma_capsulatum_144_p_08class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_144_p_08class.tsv")
Histoplasma_capsulatum_144_p_08class <- Histoplasma_capsulatum_144_p_08class %>% mutate('Histoplasma capsulatum 144_p_08'= percentage)
Histoplasma_capsulatum_144_p_08class1 <- Histoplasma_capsulatum_144_p_08class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_144_p_08class_N <- Histoplasma_capsulatum_144_p_08class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_1517_p_17 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_1517_p_17.RM/Histoplasma_capsulatum_1517_p_17.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_1517_p_17) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_1517_p_17 <- Histoplasma_capsulatum_1517_p_17 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_1517_p_17,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_1517_p_17.tsv")

Histoplasma_capsulatum_1517_p_17class<- Histoplasma_capsulatum_1517_p_17 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28646603)
write_tsv(Histoplasma_capsulatum_1517_p_17class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_1517_p_17class.tsv")
Histoplasma_capsulatum_1517_p_17class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_1517_p_17class.tsv")
Histoplasma_capsulatum_1517_p_17class <- Histoplasma_capsulatum_1517_p_17class %>% mutate('Histoplasma capsulatum 1517_p_17'= percentage)
#Histoplasma_capsulatum_1517_p_17class <- Histoplasma_capsulatum_1517_p_17class %>% select(-c('Histoplasma capsulatum 1517_p_17'))
Histoplasma_capsulatum_1517_p_17class1 <- Histoplasma_capsulatum_1517_p_17class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_1517_p_17class_N <- Histoplasma_capsulatum_1517_p_17class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_256_P_18 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_256_P_18.RM/Histoplasma_capsulatum_256_P_18.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_256_P_18) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_256_P_18 <- Histoplasma_capsulatum_256_P_18 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_256_P_18,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_256_P_18.tsv")

Histoplasma_capsulatum_256_P_18class<- Histoplasma_capsulatum_256_P_18 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28534158)
write_tsv(Histoplasma_capsulatum_256_P_18class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_256_P_18class.tsv")
Histoplasma_capsulatum_256_P_18class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_256_P_18class.tsv")
Histoplasma_capsulatum_256_P_18class <- Histoplasma_capsulatum_256_P_18class %>% mutate('Histoplasma capsulatum 256_P_18'= n)
Histoplasma_capsulatum_256_P_18class1 <- Histoplasma_capsulatum_256_P_18class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_256_P_18class_N <- Histoplasma_capsulatum_256_P_18class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_316_p_10 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_316_p_10.RM/Histoplasma_capsulatum_316_p_10.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_316_p_10) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_316_p_10 <- Histoplasma_capsulatum_316_p_10 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_316_p_10,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_316_p_10.tsv")

Histoplasma_capsulatum_316_p_10class<- Histoplasma_capsulatum_316_p_10 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28625580)
write_tsv(Histoplasma_capsulatum_316_p_10class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_316_p_10class.tsv")
Histoplasma_capsulatum_316_p_10class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_316_p_10class.tsv")
Histoplasma_capsulatum_316_p_10class <- Histoplasma_capsulatum_316_p_10class %>% mutate('Histoplasma capsulatum 316_p_10'= percentage)
#Histoplasma_capsulatum_316_p_10class <- Histoplasma_capsulatum_316_p_10class %>% select(-c('Histoplasma capsulatum 316_p_10'))
Histoplasma_capsulatum_316_p_10class1 <- Histoplasma_capsulatum_316_p_10class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_316_p_10class_N <- Histoplasma_capsulatum_316_p_10class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_327_P_12 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_327_P_12.RM/Histoplasma_capsulatum_327_P_12.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_327_P_12) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_327_P_12 <- Histoplasma_capsulatum_327_P_12 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_327_P_12,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_327_P_12.tsv")

Histoplasma_capsulatum_327_P_12class<- Histoplasma_capsulatum_327_P_12 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28431708)
write_tsv(Histoplasma_capsulatum_327_P_12class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_327_P_12class.tsv")
Histoplasma_capsulatum_327_P_12class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_327_P_12class.tsv")
Histoplasma_capsulatum_327_P_12class <- Histoplasma_capsulatum_327_P_12class %>% mutate('Histoplasma capsulatum 327_P_12'= percentage)
#Histoplasma_capsulatum_327_P_12class <- Histoplasma_capsulatum_327_P_12class %>% select(-c('Histoplasma capsulatum 327_P_12'))
Histoplasma_capsulatum_327_P_12class1 <- Histoplasma_capsulatum_327_P_12class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_327_P_12class_N <- Histoplasma_capsulatum_327_P_12class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_343_p_18 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_343_p_18.RM/Histoplasma_capsulatum_343_p_18.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_343_p_18) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_343_p_18 <- Histoplasma_capsulatum_343_p_18 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_343_p_18,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_343_p_18.tsv")

Histoplasma_capsulatum_343_p_18class<- Histoplasma_capsulatum_343_p_18 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28732569)
write_tsv(Histoplasma_capsulatum_343_p_18class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_343_p_18class.tsv")
Histoplasma_capsulatum_343_p_18class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_343_p_18class.tsv")
Histoplasma_capsulatum_343_p_18class <- Histoplasma_capsulatum_343_p_18class %>% mutate('Histoplasma capsulatum 343_p_18'= n)
Histoplasma_capsulatum_343_p_18class1 <- Histoplasma_capsulatum_343_p_18class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_343_p_18class_N <- Histoplasma_capsulatum_343_p_18class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_388_p_11 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_388_p_11.RM/Histoplasma_capsulatum_388_p_11.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_388_p_11) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_388_p_11 <- Histoplasma_capsulatum_388_p_11 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_388_p_11,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_388_p_11.tsv")

Histoplasma_capsulatum_388_p_11class<- Histoplasma_capsulatum_388_p_11 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28579830)
write_tsv(Histoplasma_capsulatum_388_p_11class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_388_p_11class.tsv")
Histoplasma_capsulatum_388_p_11class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_388_p_11class.tsv")
Histoplasma_capsulatum_388_p_11class <- Histoplasma_capsulatum_388_p_11class %>% mutate('Histoplasma capsulatum 388_p_11'= n)
Histoplasma_capsulatum_388_p_11class1 <- Histoplasma_capsulatum_388_p_11class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_388_p_11class_N <- Histoplasma_capsulatum_388_p_11class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_ES2_83Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_ES2_83Z.RM/Histoplasma_capsulatum_ES2_83Z.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_ES2_83Z) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_ES2_83Z <- Histoplasma_capsulatum_ES2_83Z %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_ES2_83Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_83Z.tsv")

Histoplasma_capsulatum_ES2_83Zclass<- Histoplasma_capsulatum_ES2_83Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/33292477)
write_tsv(Histoplasma_capsulatum_ES2_83Zclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_83Zclass.tsv")
Histoplasma_capsulatum_ES2_83Zclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_83Zclass.tsv")
Histoplasma_capsulatum_ES2_83Zclass <- Histoplasma_capsulatum_ES2_83Zclass %>% mutate('Histoplasma capsulatum ES2_83Z'= percentage)
#Histoplasma_capsulatum_ES2_83Zclass <- Histoplasma_capsulatum_ES2_83Zclass %>% select(-c('Histoplasma capsulatum ES2_83Z'))
Histoplasma_capsulatum_ES2_83Zclass1 <- Histoplasma_capsulatum_ES2_83Zclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_ES2_83Zclass_N <- Histoplasma_capsulatum_ES2_83Zclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_ES2_85Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_ES2_85Z.RM/Histoplasma_capsulatum_ES2_85Z.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_ES2_85Z) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_ES2_85Z <- Histoplasma_capsulatum_ES2_85Z %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_ES2_85Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_85Z.tsv")

Histoplasma_capsulatum_ES2_85Zclass<- Histoplasma_capsulatum_ES2_85Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28184349)
write_tsv(Histoplasma_capsulatum_ES2_85Zclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_85Zclass.tsv")
Histoplasma_capsulatum_ES2_85Zclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_85Zclass.tsv")
Histoplasma_capsulatum_ES2_85Zclass <- Histoplasma_capsulatum_ES2_85Zclass %>% mutate('Histoplasma capsulatum ES2_85Z'= n)
Histoplasma_capsulatum_ES2_85Zclass1 <- Histoplasma_capsulatum_ES2_85Zclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_ES2_85Zclass_N <- Histoplasma_capsulatum_ES2_85Zclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_ES2_86Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_ES2_86Z.RM/Histoplasma_capsulatum_ES2_86Z.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_ES2_86Z) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_ES2_86Z <- Histoplasma_capsulatum_ES2_86Z %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_ES2_86Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_86Z.tsv")

Histoplasma_capsulatum_ES2_86Zclass<- Histoplasma_capsulatum_ES2_86Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/40716201)
write_tsv(Histoplasma_capsulatum_ES2_86Zclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_86Zclass.tsv")
Histoplasma_capsulatum_ES2_86Zclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_86Zclass.tsv")
Histoplasma_capsulatum_ES2_86Zclass <- Histoplasma_capsulatum_ES2_86Zclass %>% mutate('Histoplasma capsulatum ES2_86Z'= n)
Histoplasma_capsulatum_ES2_86Zclass1 <- Histoplasma_capsulatum_ES2_86Zclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_ES2_86Zclass_N <- Histoplasma_capsulatum_ES2_86Zclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_ES2_87 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_ES2_87.RM/Histoplasma_capsulatum_ES2_87.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_ES2_87) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_ES2_87 <- Histoplasma_capsulatum_ES2_87 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_ES2_87,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_87.tsv")

Histoplasma_capsulatum_ES2_87class<- Histoplasma_capsulatum_ES2_87 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/61384325)
write_tsv(Histoplasma_capsulatum_ES2_87class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_87class.tsv")
Histoplasma_capsulatum_ES2_87class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_87class.tsv")
Histoplasma_capsulatum_ES2_87class <- Histoplasma_capsulatum_ES2_87class %>% mutate('Histoplasma capsulatum ES2_87'= n)
Histoplasma_capsulatum_ES2_87class1 <- Histoplasma_capsulatum_ES2_87class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_ES2_87class_N <- Histoplasma_capsulatum_ES2_87class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_ES2_88Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_ES2_88Z.RM/Histoplasma_capsulatum_ES2_88Z.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_ES2_88Z) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_ES2_88Z <- Histoplasma_capsulatum_ES2_88Z %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_ES2_88Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_88Z.tsv")

Histoplasma_capsulatum_ES2_88Zclass<- Histoplasma_capsulatum_ES2_88Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25109007)
write_tsv(Histoplasma_capsulatum_ES2_88Zclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_88Zclass.tsv")
Histoplasma_capsulatum_ES2_88Zclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_88Zclass.tsv")
Histoplasma_capsulatum_ES2_88Zclass <- Histoplasma_capsulatum_ES2_88Zclass %>% mutate('Histoplasma capsulatum ES2_88Z'= n)
Histoplasma_capsulatum_ES2_88Zclass1 <- Histoplasma_capsulatum_ES2_88Zclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_ES2_88Zclass_N <- Histoplasma_capsulatum_ES2_88Zclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_ES2_89Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_ES2_89Z.RM/Histoplasma_capsulatum_ES2_89Z.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_ES2_89Z) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_ES2_89Z <- Histoplasma_capsulatum_ES2_89Z %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_ES2_89Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_89Z.tsv")

Histoplasma_capsulatum_ES2_89Zclass<- Histoplasma_capsulatum_ES2_89Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/34058413)
write_tsv(Histoplasma_capsulatum_ES2_89Zclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_89Zclass.tsv")
Histoplasma_capsulatum_ES2_89Zclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_89Zclass.tsv")
Histoplasma_capsulatum_ES2_89Zclass <- Histoplasma_capsulatum_ES2_89Zclass %>% mutate('Histoplasma capsulatum ES2_89Z'= n)
Histoplasma_capsulatum_ES2_89Zclass1 <- Histoplasma_capsulatum_ES2_89Zclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_ES2_89Zclass_N <- Histoplasma_capsulatum_ES2_89Zclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_ES2_90Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_ES2_90Z.RM/Histoplasma_capsulatum_ES2_90Z.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_ES2_90Z) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_ES2_90Z <- Histoplasma_capsulatum_ES2_90Z %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_ES2_90Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_90Z.tsv")

Histoplasma_capsulatum_ES2_90Zclass<- Histoplasma_capsulatum_ES2_90Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/36554438)
write_tsv(Histoplasma_capsulatum_ES2_90Zclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_90Zclass.tsv")
Histoplasma_capsulatum_ES2_90Zclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_90Zclass.tsv")
Histoplasma_capsulatum_ES2_90Zclass <- Histoplasma_capsulatum_ES2_90Zclass %>% mutate('Histoplasma capsulatum ES2_90Z'= n)
Histoplasma_capsulatum_ES2_90Zclass1 <- Histoplasma_capsulatum_ES2_90Zclass %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_ES2_90Zclass_N <- Histoplasma_capsulatum_ES2_90Zclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_SA15 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SA15.RM/Histoplasma_capsulatum_SA15.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_SA15) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_SA15 <- Histoplasma_capsulatum_ES2_90Z %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_SA15,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SA15.tsv")

Histoplasma_capsulatum_SA15class<- Histoplasma_capsulatum_SA15 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28298332)
write_tsv(Histoplasma_capsulatum_SA15class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SA15class.tsv")
Histoplasma_capsulatum_SA15class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SA15class.tsv")
Histoplasma_capsulatum_SA15class <- Histoplasma_capsulatum_SA15class %>% mutate('Histoplasma capsulatum SA15'= n)
Histoplasma_capsulatum_SA15class1 <- Histoplasma_capsulatum_SA15class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_SA15class_N <- Histoplasma_capsulatum_SA15class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB053 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB053.RM/Histoplasma_capsulatum_CB053.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB053) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB053 <- Histoplasma_capsulatum_CB053%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB053,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB053.tsv")

Histoplasma_capsulatum_CB053class<- Histoplasma_capsulatum_CB053 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26667782)
write_tsv(Histoplasma_capsulatum_CB053class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB053class.tsv")
Histoplasma_capsulatum_CB053class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB053class.tsv")
Histoplasma_capsulatum_CB053class <- Histoplasma_capsulatum_CB053class %>% mutate('Histoplasma capsulatum CB053'= n)
Histoplasma_capsulatum_CB053class1 <- Histoplasma_capsulatum_CB053class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB053class_N <- Histoplasma_capsulatum_CB053class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB062 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB062.RM/Histoplasma_capsulatum_CB062.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB062) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB062 <- Histoplasma_capsulatum_CB062%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB062,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB062.tsv")

Histoplasma_capsulatum_CB062class<- Histoplasma_capsulatum_CB062 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28817514)
write_tsv(Histoplasma_capsulatum_CB062class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB062class.tsv")
Histoplasma_capsulatum_CB062class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB062class.tsv")
Histoplasma_capsulatum_CB062class <- Histoplasma_capsulatum_CB062class %>% mutate('Histoplasma capsulatum CB062'= n)
Histoplasma_capsulatum_CB062class1 <- Histoplasma_capsulatum_CB062class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB062class_N <- Histoplasma_capsulatum_CB062class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB063<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB063.RM/Histoplasma_capsulatum_CB063.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB063) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB063 <- Histoplasma_capsulatum_CB063%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB063,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB063.tsv")

Histoplasma_capsulatum_CB063class<- Histoplasma_capsulatum_CB063 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28864442)
write_tsv(Histoplasma_capsulatum_CB063class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB063class.tsv")
Histoplasma_capsulatum_CB063class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB063class.tsv")
Histoplasma_capsulatum_CB063class <- Histoplasma_capsulatum_CB063class %>% mutate('Histoplasma capsulatum CB063'= percentage)
#Histoplasma_capsulatum_CB063class <- Histoplasma_capsulatum_CB063class %>% select (-c('Histoplasma capsulatum CB063'))
Histoplasma_capsulatum_CB063class1 <- Histoplasma_capsulatum_CB063class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB063class_N <- Histoplasma_capsulatum_CB063class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB066<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB066.RM/Histoplasma_capsulatum_CB066.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB066) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB066 <- Histoplasma_capsulatum_CB066%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB066,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB066.tsv")

Histoplasma_capsulatum_CB066class<- Histoplasma_capsulatum_CB066 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29117842)
write_tsv(Histoplasma_capsulatum_CB066class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB066class.tsv")
Histoplasma_capsulatum_CB066class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB066class.tsv")
Histoplasma_capsulatum_CB066class <- Histoplasma_capsulatum_CB066class %>% mutate('Histoplasma capsulatum CB066'= n)
Histoplasma_capsulatum_CB066class1 <- Histoplasma_capsulatum_CB066class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB066class_N <- Histoplasma_capsulatum_CB066class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB180<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB180.RM/Histoplasma_capsulatum_CB180.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB180) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB180 <- Histoplasma_capsulatum_CB180%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB180,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB180.tsv")

Histoplasma_capsulatum_CB180class<- Histoplasma_capsulatum_CB180 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29316270)
write_tsv(Histoplasma_capsulatum_CB180class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB180class.tsv")
Histoplasma_capsulatum_CB180class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB180class.tsv")
Histoplasma_capsulatum_CB180class <- Histoplasma_capsulatum_CB180class %>% mutate('Histoplasma capsulatum CB180'= n)
Histoplasma_capsulatum_CB180class1 <- Histoplasma_capsulatum_CB180class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB180class_N <- Histoplasma_capsulatum_CB180class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB0522<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB052-2.RM/Histoplasma_capsulatum_CB052-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB0522) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB0522 <- Histoplasma_capsulatum_CB0522%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB0522,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB052-2.tsv")

Histoplasma_capsulatum_CB0522class<- Histoplasma_capsulatum_CB0522 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/17119958)
write_tsv(Histoplasma_capsulatum_CB0522class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB052-2class.tsv")
Histoplasma_capsulatum_CB0522class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB052-2class.tsv")
Histoplasma_capsulatum_CB0522class <- Histoplasma_capsulatum_CB0522class %>% mutate('Histoplasma capsulatum CB052-2'= n)
Histoplasma_capsulatum_CB0522class1 <- Histoplasma_capsulatum_CB0522class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB0522class_N <- Histoplasma_capsulatum_CB0522class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB0532<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB053-2.RM/Histoplasma_capsulatum_CB053-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB0532) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB0532 <- Histoplasma_capsulatum_CB0532%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB0532,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB053-2.tsv")

Histoplasma_capsulatum_CB0532class<- Histoplasma_capsulatum_CB0532 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23649534)
write_tsv(Histoplasma_capsulatum_CB0532class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB053-2class.tsv")
Histoplasma_capsulatum_CB0532class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB053-2class.tsv")
Histoplasma_capsulatum_CB0532class <- Histoplasma_capsulatum_CB0532class %>% mutate('Histoplasma capsulatum CB053-2'= n)
Histoplasma_capsulatum_CB0532class1 <- Histoplasma_capsulatum_CB0532class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB0532class_N <- Histoplasma_capsulatum_CB0532class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB0552<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB055-2.RM/Histoplasma_capsulatum_CB055-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB0552) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB0552 <- Histoplasma_capsulatum_CB0552%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB0552,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB055-2.tsv")

Histoplasma_capsulatum_CB0552class<- Histoplasma_capsulatum_CB0552 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23847241)
write_tsv(Histoplasma_capsulatum_CB0552class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB055-2class.tsv")
Histoplasma_capsulatum_CB0552class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB055-2class.tsv")
Histoplasma_capsulatum_CB0552class <- Histoplasma_capsulatum_CB0552class %>% mutate('Histoplasma capsulatum CB055-2'= percentage)
Histoplasma_capsulatum_CB0552class <- Histoplasma_capsulatum_CB0552class %>% select (-c('Histoplasma capsulatum CB055-2'))
Histoplasma_capsulatum_CB0552class1 <- Histoplasma_capsulatum_CB0552class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB0552class_N <- Histoplasma_capsulatum_CB0552class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB0622<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB062-2.RM/Histoplasma_capsulatum_CB062-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB0622) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB0622 <- Histoplasma_capsulatum_CB0622%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB0622,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB062-2.tsv")

Histoplasma_capsulatum_CB0622class<- Histoplasma_capsulatum_CB0622 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28739478)
write_tsv(Histoplasma_capsulatum_CB0622class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB062-2class.tsv")
Histoplasma_capsulatum_CB0622class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB062-2class.tsv")
Histoplasma_capsulatum_CB0622class <- Histoplasma_capsulatum_CB0622class %>% mutate('Histoplasma capsulatum CB062-2'=percentage)
Histoplasma_capsulatum_CB0622class1 <- Histoplasma_capsulatum_CB0622class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB0622class_N <- Histoplasma_capsulatum_CB0622class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB0632<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB063-2.RM/Histoplasma_capsulatum_CB063-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB0632) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB0632 <- Histoplasma_capsulatum_CB0632%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB0632,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB063-2.tsv")

Histoplasma_capsulatum_CB0632class<- Histoplasma_capsulatum_CB0632 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28655085)
write_tsv(Histoplasma_capsulatum_CB0632class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB063-2class.tsv")
Histoplasma_capsulatum_CB0632class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB063-2class.tsv")
Histoplasma_capsulatum_CB0632class <- Histoplasma_capsulatum_CB0632class %>% mutate('Histoplasma capsulatum CB063-2'= n)
Histoplasma_capsulatum_CB0632class1 <- Histoplasma_capsulatum_CB0632class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB0632class_N <- Histoplasma_capsulatum_CB0632class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB0642<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB064-2.RM/Histoplasma_capsulatum_CB064-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB0642) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB0642 <- Histoplasma_capsulatum_CB0642%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB0642,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB064-2.tsv")

Histoplasma_capsulatum_CB0642class<- Histoplasma_capsulatum_CB0642 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28635935)
write_tsv(Histoplasma_capsulatum_CB0642class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB064-2class.tsv")
Histoplasma_capsulatum_CB0642class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB064-2class.tsv")
Histoplasma_capsulatum_CB0642class <- Histoplasma_capsulatum_CB0642class %>% mutate('Histoplasma capsulatum CB064-2'= n)
Histoplasma_capsulatum_CB0642class1 <- Histoplasma_capsulatum_CB0642class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB0642class_N <- Histoplasma_capsulatum_CB0642class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB0652<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB065-2.RM/Histoplasma_capsulatum_CB065-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB0652) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB0652 <- Histoplasma_capsulatum_CB0652%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB0652,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB065-2.tsv")

Histoplasma_capsulatum_CB0652class<- Histoplasma_capsulatum_CB0652 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29229412)
write_tsv(Histoplasma_capsulatum_CB0652class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB065-2class.tsv")
Histoplasma_capsulatum_CB0652class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB065-2class.tsv")
Histoplasma_capsulatum_CB0652class <- Histoplasma_capsulatum_CB0652class %>% mutate('Histoplasma capsulatum CB065-2'= n)
Histoplasma_capsulatum_CB0652class1 <- Histoplasma_capsulatum_CB0652class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB0652class_N <- Histoplasma_capsulatum_CB0652class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB0662<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB066-2.RM/Histoplasma_capsulatum_CB066-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB0662) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB0662 <- Histoplasma_capsulatum_CB0662%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB0662,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB066-2.tsv")

Histoplasma_capsulatum_CB0662class<- Histoplasma_capsulatum_CB0662 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28630619)
write_tsv(Histoplasma_capsulatum_CB0662class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB066-2class.tsv")
Histoplasma_capsulatum_CB0662class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB066-2class.tsv")
Histoplasma_capsulatum_CB0662class <- Histoplasma_capsulatum_CB0662class %>% mutate('Histoplasma capsulatum CB066-2'= n)
Histoplasma_capsulatum_CB0662class1 <- Histoplasma_capsulatum_CB0662class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB0662class_N <- Histoplasma_capsulatum_CB0662class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB1682<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB168-2.RM/Histoplasma_capsulatum_CB168-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB1682) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB1682 <- Histoplasma_capsulatum_CB1682%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB1682,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB168-2.tsv")

Histoplasma_capsulatum_CB1682class<- Histoplasma_capsulatum_CB1682 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28670342)
write_tsv(Histoplasma_capsulatum_CB1682class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB168-2class.tsv")
Histoplasma_capsulatum_CB1682class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB168-2class.tsv")
Histoplasma_capsulatum_CB1682class <- Histoplasma_capsulatum_CB1682class %>% mutate('Histoplasma capsulatum CB168-2'= n)
Histoplasma_capsulatum_CB1682class1 <- Histoplasma_capsulatum_CB1682class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB1682class_N <- Histoplasma_capsulatum_CB1682class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB1742<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB174-2.RM/Histoplasma_capsulatum_CB174-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB1742) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB1742 <- Histoplasma_capsulatum_CB1742%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB1742,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB174-2.tsv")

Histoplasma_capsulatum_CB1742class<- Histoplasma_capsulatum_CB1742 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28472442)
write_tsv(Histoplasma_capsulatum_CB1742class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB174-2class.tsv")
Histoplasma_capsulatum_CB1742class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB174-2class.tsv")
Histoplasma_capsulatum_CB1742class <- Histoplasma_capsulatum_CB1742class %>% mutate('Histoplasma capsulatum CB174-2'= n)
Histoplasma_capsulatum_CB1742class1 <- Histoplasma_capsulatum_CB1742class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB1742class_N <- Histoplasma_capsulatum_CB1742class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB1802<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB180-2.RM/Histoplasma_capsulatum_CB180-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB1802) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB1802 <- Histoplasma_capsulatum_CB1802%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB1802,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB180-2.tsv")

Histoplasma_capsulatum_CB1802class<- Histoplasma_capsulatum_CB1802 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28608116)
write_tsv(Histoplasma_capsulatum_CB1802class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB180-2class.tsv")
Histoplasma_capsulatum_CB1802class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB180-2class.tsv")
Histoplasma_capsulatum_CB1802class <- Histoplasma_capsulatum_CB1802class %>% mutate('Histoplasma capsulatum CB180-2'= n)
Histoplasma_capsulatum_CB1802class1 <- Histoplasma_capsulatum_CB1802class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB1802class_N <- Histoplasma_capsulatum_CB1802class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB1862<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB186-2.RM/Histoplasma_capsulatum_CB186-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB1862) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB1862 <- Histoplasma_capsulatum_CB1862%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB1862,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB186-2.tsv")

Histoplasma_capsulatum_CB1862class<- Histoplasma_capsulatum_CB1862 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28220179)
write_tsv(Histoplasma_capsulatum_CB1862class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB186-2class.tsv")
Histoplasma_capsulatum_CB1862class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB186-2class.tsv")
Histoplasma_capsulatum_CB1862class <- Histoplasma_capsulatum_CB1862class %>% mutate('Histoplasma capsulatum CB186-2'= n)
Histoplasma_capsulatum_CB1862class1 <- Histoplasma_capsulatum_CB1862class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB1862class_N <- Histoplasma_capsulatum_CB1862class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_CB1922<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_CB192-2.RM/Histoplasma_capsulatum_CB192-2.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_CB1922) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_CB1922 <- Histoplasma_capsulatum_CB1922%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_CB1922,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB192-2.tsv")

Histoplasma_capsulatum_CB1922class<- Histoplasma_capsulatum_CB1922 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28781169)
write_tsv(Histoplasma_capsulatum_CB1922class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB192-2class.tsv")
Histoplasma_capsulatum_CB1922class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB192-2class.tsv")
Histoplasma_capsulatum_CB1922class <- Histoplasma_capsulatum_CB1922class %>% mutate('Histoplasma capsulatum CB192-2'= n)
Histoplasma_capsulatum_CB1922class1 <- Histoplasma_capsulatum_CB1922class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_CB1922class_N <- Histoplasma_capsulatum_CB1922class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11.RM/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11 <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11.tsv")

Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class<- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28690043)
write_tsv(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class.tsv")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class.tsv")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class %>% mutate('Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11'= percentage)
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class1 <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class_N <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14.RM/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14 <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14.tsv")

Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class<- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28326430)
write_tsv(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class.tsv")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class.tsv")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class %>% mutate('Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14'= n)
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class1 <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class_N <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16.RM/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16 <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16.tsv")

Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class<- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28332259)
write_tsv(Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class.tsv")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class.tsv")
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class %>% mutate('Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16'= n)
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class1 <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class %>% select(-c(n,total,percentage))
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class_N <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_102Nam2_505<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_102-Nam2_505.RM/Histoplasma_mississippiense_II_SECH_102-Nam2_505.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_102Nam2_505) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_102Nam2_505 <- Histoplasma_mississippienseII_SECH_102Nam2_505%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_102Nam2_505,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_102Nam2_505.tsv")

Histoplasma_mississippienseII_SECH_102Nam2_505class<- Histoplasma_mississippienseII_SECH_102Nam2_505 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/27544444)
write_tsv(Histoplasma_mississippienseII_SECH_102Nam2_505class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_102Nam2_505class.tsv")
Histoplasma_mississippienseII_SECH_102Nam2_505class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_102Nam2_505class.tsv")
Histoplasma_mississippienseII_SECH_102Nam2_505class <- Histoplasma_mississippienseII_SECH_102Nam2_505class %>% mutate('Histoplasma mississippiense II SECH_102Nam2_505'= n)
Histoplasma_mississippienseII_SECH_102Nam2_505class1 <- Histoplasma_mississippienseII_SECH_102Nam2_505class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_102Nam2_505class_N <- Histoplasma_mississippienseII_SECH_102Nam2_505class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_81Nam1_WU24<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_81-Nam1_WU24.RM/Histoplasma_mississippiense_II_SECH_81-Nam1_WU24.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_81Nam1_WU24) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_81Nam1_WU24 <- Histoplasma_mississippienseII_SECH_81Nam1_WU24%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_81Nam1_WU24,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_81-Nam1_WU24.tsv")

Histoplasma_mississippienseII_SECH_81Nam1_WU24class<- Histoplasma_mississippienseII_SECH_81Nam1_WU24 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25518040)
write_tsv(Histoplasma_mississippienseII_SECH_81Nam1_WU24class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_81-Nam1_WU24class.tsv")
Histoplasma_mississippienseII_SECH_81Nam1_WU24class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_81-Nam1_WU24class.tsv")
Histoplasma_mississippienseII_SECH_81Nam1_WU24class <- Histoplasma_mississippienseII_SECH_81Nam1_WU24class %>% mutate('Histoplasma mississippiense II SECH_81-Nam1_WU24'= n)
Histoplasma_mississippienseII_SECH_81Nam1_WU24class1 <- Histoplasma_mississippienseII_SECH_81Nam1_WU24class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_81Nam1_WU24class_N <- Histoplasma_mississippienseII_SECH_81Nam1_WU24class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_83Nam1_CI_7<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_83-Nam1_CI_7.RM/Histoplasma_mississippiense_II_SECH_83-Nam1_CI_7.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_83Nam1_CI_7) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_83Nam1_CI_7 <- Histoplasma_mississippienseII_SECH_83Nam1_CI_7%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_83Nam1_CI_7,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_83-Nam1_CI_7.tsv")

Histoplasma_mississippienseII_SECH_83Nam1_CI_7class<- Histoplasma_mississippienseII_SECH_83Nam1_CI_7 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25776464)
write_tsv(Histoplasma_mississippienseII_SECH_83Nam1_CI_7class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_83-Nam1_CI_7class.tsv")
Histoplasma_mississippienseII_SECH_83Nam1_CI_7class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_83-Nam1_CI_7class.tsv")
Histoplasma_mississippienseII_SECH_83Nam1_CI_7class <- Histoplasma_mississippienseII_SECH_83Nam1_CI_7class %>% mutate('Histoplasma mississippiense II SECH_83-Nam1_CI_7'= n)
Histoplasma_mississippienseII_SECH_83Nam1_CI_7class1 <- Histoplasma_mississippienseII_SECH_83Nam1_CI_7class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_83Nam1_CI_7class_N <- Histoplasma_mississippienseII_SECH_83Nam1_CI_7class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_84Nam1_CI_19<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_84-Nam1_CI_19.RM/Histoplasma_mississippiense_II_SECH_84-Nam1_CI_19.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_84Nam1_CI_19) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_84Nam1_CI_19 <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_84Nam1_CI_19,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_84-Nam1_CI_19.tsv")

Histoplasma_mississippienseII_SECH_84Nam1_CI_19class<- Histoplasma_mississippienseII_SECH_84Nam1_CI_19 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26524728)
write_tsv(Histoplasma_mississippienseII_SECH_84Nam1_CI_19class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_84-Nam1_CI_19class.tsv")
Histoplasma_mississippienseII_SECH_84Nam1_CI_19class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_84-Nam1_CI_19class.tsv")
Histoplasma_mississippienseII_SECH_84Nam1_CI_19class <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19class %>% mutate('Histoplasma mississippiense II SECH_84-Nam1_CI_19'= n)
Histoplasma_mississippienseII_SECH_84Nam1_CI_19class1 <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_84Nam1_CI_19class_N <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_85Nam1_CI_22<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_85-Nam1_CI_22.RM/Histoplasma_mississippiense_II_SECH_85-Nam1_CI_22.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_85Nam1_CI_22) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_85Nam1_CI_22 <- Histoplasma_mississippienseII_SECH_85Nam1_CI_22%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_85Nam1_CI_22,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_85-Nam1_CI_22.tsv")

Histoplasma_mississippienseII_SECH_85Nam1_CI_22class<- Histoplasma_mississippienseII_SECH_85Nam1_CI_22 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23145789)
write_tsv(Histoplasma_mississippienseII_SECH_85Nam1_CI_22class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_85-Nam1_CI_22class.tsv")
Histoplasma_mississippienseII_SECH_85Nam1_CI_22class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SSECH_85-Nam1_CI_22class.tsv")
Histoplasma_mississippienseII_SECH_85Nam1_CI_22class <- Histoplasma_mississippienseII_SECH_85Nam1_CI_22class %>% mutate('Histoplasma mississippiense II SECH_85-Nam1_CI_22'= n)
Histoplasma_mississippienseII_SECH_85Nam1_CI_22class1 <- Histoplasma_mississippienseII_SECH_85Nam1_CI_22class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_85Nam1_CI_22class_N <- Histoplasma_mississippienseII_SECH_85Nam1_CI_22class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_86Nam1_CI_24<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_86-Nam1_CI_24.RM/Histoplasma_mississippiense_II_SECH_86-Nam1_CI_24.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_86Nam1_CI_24) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_86Nam1_CI_24 <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_86Nam1_CI_24,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_86-Nam1_CI_24.tsv")

Histoplasma_mississippienseII_SECH_86Nam1_CI_24class<- Histoplasma_mississippienseII_SECH_86Nam1_CI_24 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23809594)
write_tsv(Histoplasma_mississippienseII_SECH_86Nam1_CI_24class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_86-Nam1_CI_24class.tsv")
Histoplasma_mississippienseII_SECH_86Nam1_CI_24class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_86-Nam1_CI_24class.tsv")
Histoplasma_mississippienseII_SECH_86Nam1_CI_24class <- Histoplasma_mississippienseII_SECH_86Nam1_CI_24class %>% mutate('Histoplasma mississippiense II SECH_86-Nam1_CI_24'= n)
Histoplasma_mississippienseII_SECH_86Nam1_CI_24class1 <- Histoplasma_mississippienseII_SECH_86Nam1_CI_24class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_86Nam1_CI_24class_N <- Histoplasma_mississippienseII_SECH_86Nam1_CI_24class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_87Nam1_CI_42<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_87-Nam1_CI_42.RM/Histoplasma_mississippiense_II_SECH_87-Nam1_CI_42.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_87Nam1_CI_42) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_87Nam1_CI_42 <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_87Nam1_CI_42,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_87-Nam1_CI_42.tsv")

Histoplasma_mississippienseII_SECH_87Nam1_CI_42class<- Histoplasma_mississippienseII_SECH_87Nam1_CI_42 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26345314)
write_tsv(Histoplasma_mississippienseII_SECH_87Nam1_CI_42class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_87-Nam1_CI_42class.tsv")
Histoplasma_mississippienseII_SECH_87Nam1_CI_42class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_87-Nam1_CI_42class.tsv")
Histoplasma_mississippienseII_SECH_87Nam1_CI_42class <- Histoplasma_mississippienseII_SECH_87Nam1_CI_42class %>% mutate('Histoplasma mississippiense II SECH_87-Nam1_CI_42'= n)
Histoplasma_mississippienseII_SECH_87Nam1_CI_42class1 <- Histoplasma_mississippienseII_SECH_87Nam1_CI_42class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_87Nam1_CI_42class_N <- Histoplasma_mississippienseII_SECH_87Nam1_CI_42class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_88Nam1_CI_43<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_88-Nam1_CI_43.RM/Histoplasma_mississippiense_II_SECH_88-Nam1_CI_43.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_88Nam1_CI_43) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_88Nam1_CI_43 <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_88Nam1_CI_43,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_88-Nam1_CI_43.tsv")

Histoplasma_mississippienseII_SECH_88Nam1_CI_43class<- Histoplasma_mississippienseII_SECH_88Nam1_CI_43 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24925558)
write_tsv(Histoplasma_mississippienseII_SECH_88Nam1_CI_43class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_88-Nam1_CI_43class.tsv")
Histoplasma_mississippienseII_SECH_88Nam1_CI_43class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_88-Nam1_CI_43class.tsv")
Histoplasma_mississippienseII_SECH_88Nam1_CI_43class <- Histoplasma_mississippienseII_SECH_88Nam1_CI_43class %>% mutate('Histoplasma mississippiense II SECH_88-Nam1_CI_43'= n)
Histoplasma_mississippienseII_SECH_88Nam1_CI_43class1 <- Histoplasma_mississippienseII_SECH_88Nam1_CI_43class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_88Nam1_CI_43class_N <- Histoplasma_mississippienseII_SECH_88Nam1_CI_43class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_89Nam1_UCLA531<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_I_SECH_89-Nam1_UCLA-531.RM/Histoplasma_mississippiense_I_SECH_89-Nam1_UCLA-531.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_89Nam1_UCLA531) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_89Nam1_UCLA531 <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_89Nam1_UCLA531,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_89-Nam1_UCLA-531.tsv")

Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class<- Histoplasma_mississippienseII_SECH_89Nam1_UCLA531 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24430709)
write_tsv(Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_89-Nam1_UCLA-531class.tsv")
Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_89-Nam1_UCLA-531class.tsv")
Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class <- Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class %>% mutate('Histoplasma mississippiense II SECH_89-Nam1_UCLA-531'= n)
Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class1 <- Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class_N <- Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_SECH_90Nam1_DOWNS<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_SECH_90-Nam1_DOWNS.RM/Histoplasma_mississippiense_II_SECH_90-Nam1_DOWNS.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_SECH_90Nam1_DOWNS) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_SECH_90Nam1_DOWNS <- Histoplasma_mississippienseII_SECH_90Nam1_DOWNS%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_SECH_90Nam1_DOWNS,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_90-Nam1_DOWNS.tsv")

Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass<- Histoplasma_mississippienseII_SECH_90Nam1_DOWNS %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24054288)
write_tsv(Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_90-Nam1_DOWNSclass.tsv")
Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_90-Nam1_DOWNSclass.tsv")
Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass <- Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass %>% mutate('Histoplasma mississippiense II SECH_90-Nam1_DOWNS'= percentage)
Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass1 <- Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass_N <- Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_HCWU24<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_HCWU24.RM/Histoplasma_mississippiense_II_HCWU24.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_HCWU24) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_HCWU24 <- Histoplasma_mississippienseII_SECH_90Nam1_DOWNS%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_HCWU24,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCWU24.tsv")

Histoplasma_mississippienseII_HCWU24class<- Histoplasma_mississippienseII_SECH_90Nam1_DOWNS %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28071263)
write_tsv(Histoplasma_mississippienseII_HCWU24class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCWU24class.tsv")
Histoplasma_mississippienseII_HCWU24class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCWU24class.tsv")
Histoplasma_mississippienseII_HCWU24class <- Histoplasma_mississippienseII_HCWU24class %>% mutate('Histoplasma mississippiense II HCWU24'= percentage)
Histoplasma_mississippienseII_HCWU24class1 <- Histoplasma_mississippienseII_HCWU24class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_HCWU24class_N <- Histoplasma_mississippienseII_HCWU24class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_HCCI_43<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_HCCI_43.RM/Histoplasma_mississippiense_II_HCCI_43.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_HCCI_43) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_HCCI_43 <- Histoplasma_mississippienseII_HCCI_43%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_HCCI_43,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCCI_43.tsv")

Histoplasma_mississippienseII_HCCI_43class<- Histoplasma_mississippienseII_HCCI_43 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28137748)
write_tsv(Histoplasma_mississippienseII_HCCI_43class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCCI_43class.tsv")
Histoplasma_mississippienseII_HCCI_43class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCCI_43class.tsv")
Histoplasma_mississippienseII_HCCI_43class <- Histoplasma_mississippienseII_HCCI_43class %>% mutate('Histoplasma mississippiense II HCCI_43'= n)
Histoplasma_mississippienseII_HCCI_43class1 <- Histoplasma_mississippienseII_HCCI_43class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_HCCI_43class_N <- Histoplasma_mississippienseII_HCCI_43class %>% select(-c(n,total,percentage))

##
Histoplasma_mississippienseII_HCCI_19<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_mississippiense_II_HCCI_19.RM/Histoplasma_mississippiense_II_HCCI_19.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_mississippienseII_HCCI_19) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_mississippienseII_HCCI_19 <- Histoplasma_mississippienseII_HCCI_19%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_mississippienseII_HCCI_19,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCCI_19.tsv")

Histoplasma_mississippienseII_HCCI_19class<- Histoplasma_mississippienseII_HCCI_19 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/27643440)
write_tsv(Histoplasma_mississippienseII_HCCI_19class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCCI_19class.tsv")
Histoplasma_mississippienseII_HCCI_19class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_HCCI_19class.tsv")
Histoplasma_mississippienseII_HCCI_19class <- Histoplasma_mississippienseII_HCCI_19class %>% mutate('Histoplasma mississippiense II HCCI_19'= n)
Histoplasma_mississippienseII_HCCI_19class1 <- Histoplasma_mississippienseII_HCCI_19class %>% select(-c(n,total,percentage))
Histoplasma_mississippienseII_HCCI_19class_N <- Histoplasma_mississippienseII_HCCI_19class %>% select(-c(n,total,percentage))

##
Histoplasma_suramericanum_SECH_103Nam2_3_11G<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SECH_103-Nam2_3_11G.RM/Histoplasma_capsulatum_SECH_103-Nam2_3_11G.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_suramericanum_SECH_103Nam2_3_11G) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_suramericanum_SECH_103Nam2_3_11G <- Histoplasma_suramericanum_SECH_103Nam2_3_11G%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_suramericanum_SECH_103Nam2_3_11G,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_103-Nam2_3_11G.tsv")

Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass<- Histoplasma_suramericanum_SECH_103Nam2_3_11G %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22453624)
write_tsv(Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_103-Nam2_3_11Gclass.tsv")
Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_103-Nam2_3_11Gclass.tsv")
Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass <- Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass %>% mutate('Histoplasma suramericanum SECH_103-Nam2_3_11G'= percentage)
Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass1 <- Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass %>% select(-c(n,total,percentage))
Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass_N <- Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass %>% select(-c(n,total,percentage))

##
Histoplasma_suramericanum_SECH_104Nam2_27_14<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SECH_104-Nam2_27_14.RM/Histoplasma_capsulatum_SECH_104-Nam2_27_14.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_suramericanum_SECH_104Nam2_27_14) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_suramericanum_SECH_104Nam2_27_14 <- Histoplasma_suramericanum_SECH_104Nam2_27_14%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_suramericanum_SECH_104Nam2_27_14,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_104-Nam2_27_14.tsv")

Histoplasma_suramericanum_SECH_104Nam2_27_14class<- Histoplasma_suramericanum_SECH_104Nam2_27_14 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22213937)
write_tsv(Histoplasma_suramericanum_SECH_104Nam2_27_14class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_104-Nam2_27_14class.tsv")
Histoplasma_suramericanum_SECH_104Nam2_27_14class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_104-Nam2_27_14class.tsv")
Histoplasma_suramericanum_SECH_104Nam2_27_14class <- Histoplasma_suramericanum_SECH_104Nam2_27_14class %>% mutate('Histoplasma suramericanum SECH_104-Nam2_27_14'= n)
Histoplasma_suramericanum_SECH_104Nam2_27_14class1 <- Histoplasma_suramericanum_SECH_104Nam2_27_14class %>% select(-c(n,total,percentage))
Histoplasma_suramericanum_SECH_104Nam2_27_14class_N <- Histoplasma_suramericanum_SECH_104Nam2_27_14class %>% select(-c(n,total,percentage))

##
Histoplasma_suramericanum_SECH_105Nam2_21_14<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SECH_105-Nam2_21_14.RM/Histoplasma_capsulatum_SECH_105-Nam2_21_14.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_suramericanum_SECH_105Nam2_21_14) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_suramericanum_SECH_105Nam2_21_14 <- Histoplasma_suramericanum_SECH_105Nam2_21_14%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_suramericanum_SECH_105Nam2_21_14,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_105-Nam2_21_14.tsv")

Histoplasma_suramericanum_SECH_105Nam2_21_14class<- Histoplasma_suramericanum_SECH_105Nam2_21_14 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28048040)
write_tsv(Histoplasma_suramericanum_SECH_105Nam2_21_14class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_105-Nam2_21_14class.tsv")
Histoplasma_suramericanum_SECH_105Nam2_21_14class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_105-Nam2_21_14class.tsv")
Histoplasma_suramericanum_SECH_105Nam2_21_14class <- Histoplasma_suramericanum_SECH_105Nam2_21_14class %>% mutate('Histoplasma suramericanum SECH_105-Nam2_21_14'= n)
Histoplasma_suramericanum_SECH_105Nam2_21_14class1 <- Histoplasma_suramericanum_SECH_105Nam2_21_14class %>% select(-c(n,total,percentage))
Histoplasma_suramericanum_SECH_105Nam2_21_14class_N <- Histoplasma_suramericanum_SECH_105Nam2_21_14class %>% select(-c(n,total,percentage))

##
Histoplasma_suramericanum_HC27_14<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HC27_14.RM/Histoplasma_capsulatum_HC27_14.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_suramericanum_HC27_14) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_suramericanum_HC27_14 <- Histoplasma_suramericanum_HC27_14%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_suramericanum_HC27_14,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_HC27_14.tsv")

Histoplasma_suramericanum_HC27_14class<- Histoplasma_suramericanum_HC27_14 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24215091)
write_tsv(Histoplasma_suramericanum_HC27_14class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_HC27_14class.tsv")
Histoplasma_suramericanum_HC27_14class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_HC27_14class.tsv")
Histoplasma_suramericanum_HC27_14class <- Histoplasma_suramericanum_HC27_14class %>% mutate('Histoplasma suramericanum HC27_14'= n)
Histoplasma_suramericanum_HC27_14class1 <- Histoplasma_suramericanum_HC27_14class %>% select(-c(n,total,percentage))
Histoplasma_suramericanum_HC27_14class_N <- Histoplasma_suramericanum_HC27_14class %>% select(-c(n,total,percentage))

##
Histoplasma_suramericanum_HC21_14<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_HC21_14.RM/Histoplasma_capsulatum_HC21_14.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_suramericanum_HC21_14) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_suramericanum_HC21_14 <- Histoplasma_suramericanum_HC21_14%>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_suramericanum_HC21_14,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_HC21_14.tsv")

Histoplasma_suramericanum_HC21_14class<- Histoplasma_suramericanum_HC21_14 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28610665)
write_tsv(Histoplasma_suramericanum_HC21_14class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_HC21_14class.tsv")
Histoplasma_suramericanum_HC21_14class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_HC21_14class.tsv")
Histoplasma_suramericanum_HC21_14class <- Histoplasma_suramericanum_HC21_14class %>% mutate('Histoplasma suramericanum HC21_14'= percentage)
Histoplasma_suramericanum_HC21_14class1 <- Histoplasma_suramericanum_HC21_14class %>% select(-c(n,total,percentage))
Histoplasma_suramericanum_HC21_14class_N <- Histoplasma_suramericanum_HC21_14class %>% select(-c(n,total,percentage))

Histoplasma_suramericanum_HC21_14SUM<- Histoplasma_suramericanum_HC21_14 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length))
write_tsv(Histoplasma_suramericanum_HC21_14class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_HC21_14SUM.tsv")

##
Histoplasma_capsulatum_SECH_100Nam2_G186A<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SECH_100-Nam2_G186A.RM/Histoplasma_capsulatum_SECH_100-Nam2_G186A.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_SECH_100Nam2_G186A) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_SECH_100Nam2_G186A <- Histoplasma_capsulatum_SECH_100Nam2_G186A %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_SECH_100Nam2_G186A,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_100-Nam2_G186A.tsv")

Histoplasma_capsulatum_SECH_100Nam2_G186Aclass<- Histoplasma_capsulatum_SECH_100Nam2_G186A %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/15609890)
write_tsv(Histoplasma_capsulatum_SECH_100Nam2_G186Aclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_100-Nam2_G186Aclass.tsv")
Histoplasma_capsulatum_SECH_100Nam2_G186Aclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_100-Nam2_G186Aclass.tsv")
Histoplasma_capsulatum_SECH_100Nam2_G186Aclass <- Histoplasma_capsulatum_SECH_100Nam2_G186Aclass %>% mutate('Histoplasma capsulatum SECH_100-Nam2_G186A'= percentage)
Histoplasma_capsulatum_SECH_100Nam2_G186Aclass1 <- Histoplasma_capsulatum_SECH_100Nam2_G186Aclass %>% select(-c(n,total,percentage))
write_tsv(Histoplasma_capsulatum_SECH_100Nam2_G186Aclass1,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_100-Nam2_G186ASUM.tsv")

########## New Root Organisms
Blastomyces_parvus_UAMH130<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Blastomyces_parvus_UAMH130.RM/Blastomyces_parvus_UAMH130.scaffolds.fa.out",skip=3,col_names = F)
colnames(Blastomyces_parvus_UAMH130) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Blastomyces_parvus_UAMH130 <- Blastomyces_parvus_UAMH130 %>% mutate(length=endq-beginq +1) 
write_tsv(Blastomyces_parvus_UAMH130,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Blastomyces_parvus_UAMH130.tsv")

Blastomyces_parvus_UAMH130class<- Blastomyces_parvus_UAMH130 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/27837035)
write_tsv(Blastomyces_parvus_UAMH130class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Blastomyces_parvus_UAMH130class.tsv")
Blastomyces_parvus_UAMH130class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Blastomyces_parvus_UAMH130class.tsv")
Blastomyces_parvus_UAMH130class <- Blastomyces_parvus_UAMH130class %>% mutate('Blastomyces parvus UAMH130'= percentage)
Blastomyces_parvus_UAMH130class1 <- Blastomyces_parvus_UAMH130class %>% select(-c(n,total,percentage))
write_tsv(Blastomyces_parvus_UAMH130class1,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Blastomyces_parvus_UAMH130SUM.tsv")

##
Emmonsia_crescens_UAMH4076<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Emmonsia_crescens_UAMH4076.RM/Emmonsia_crescens_UAMH4076.scaffolds.fa.out",skip=3,col_names = F)
colnames(Emmonsia_crescens_UAMH4076) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Emmonsia_crescens_UAMH4076 <- Emmonsia_crescens_UAMH4076 %>% mutate(length=endq-beginq +1) 
write_tsv(Emmonsia_crescens_UAMH4076,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Emmonsia_crescens_UAMH4076.tsv")

Emmonsia_crescens_UAMH4076class<- Emmonsia_crescens_UAMH4076 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/33754519)
write_tsv(Emmonsia_crescens_UAMH4076class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Emmonsia_crescens_UAMH4076class.tsv")
Emmonsia_crescens_UAMH4076class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Emmonsia_crescens_UAMH4076class.tsv")
Emmonsia_crescens_UAMH4076class <- Emmonsia_crescens_UAMH4076class %>% mutate('Emmonsia crescens UAMH4076'= percentage)
Emmonsia_crescens_UAMH4076class1 <- Emmonsia_crescens_UAMH4076class %>% select(-c(n,total,percentage))
write_tsv(Emmonsia_crescens_UAMH4076class1,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Emmonsia_crescens_UAMH4076SUM.tsv")

#Combine results into one large table
Total_table<-join_all(list(Histoplasma_ohiense_HCCI_17class_N, Histoplasma_ohiense_HCCI_6class_N, Histoplasma_ohiense_HCG217Bclass_N, Histoplasma_ohiense_SECH_82.Nam1_CI_4class_N, Histoplasma_ohiense_SECH_91.Nam2_G217Bclass_N, Histoplasma_ohiense_SECH_92.Nam2_G222Bclass_N, Histoplasma_ohiense_SECH_93.Nam2_CI_6class_N, Histoplasma_ohiense_SECH_94.Nam2_CI_9class_N, Histoplasma_ohiense_SECH_95.Nam2_CI_10class_N, Histoplasma_ohiense_SECH_96.Nam2_CI_17class_N, Histoplasma_ohiense_SECH_97.Nam2_CI_18class_N, Histoplasma_ohiense_SECH_98.Nam2_CI_30class_N, Histoplasma_ohiense_SECH_99.Nam2_CI_35class_N, Histoplasma_ohiensis_Hc1986class_N, Histoplasma_capsulatum_19VMG15class_N, Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass_N, Histoplasma_capsulatum_Histo485P20class_N, Histoplasma_capsulatum_JB_01752Hc_01752class_N, Histoplasma_capsulatum_JB_021091Hc_021091class_N, Histoplasma_capsulatum_JB_031837Hc_031837class_N, Histoplasma_capsulatum_JB_042430Hc_042430class_N, Histoplasma_capsulatum_JB_062632Hc_062632class_N, Histoplasma_capsulatum_JB_062775Hc_062775class_N, Histoplasma_capsulatum_JB_073129Hc_073129class_N, Histoplasma_capsulatum_JB_083285_2Hc_083285_2class_N, Histoplasma_capsulatum_SECH_101Nam2_G184Aclass_N, Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass_N, Histoplasma_capsulatum_SECH_109class_N, Histoplasma_capsulatum_SECH_110class_N, Histoplasma_capsulatum_HCH143class_N, Histoplasma_capsulatum_HCG186Aclass_N, Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class_N, Histoplasma_capsulatum_HISSPFGTRO0285class_N, Histoplasma_capsulatum_HISSPFGPSO2043class_N, Histoplasma_capsulatum_HISSPFGPIE2055class_N, Histoplasma_capsulatum_HISSPFGPIA2052class_N,Histoplasma_capsulatum_HISSPFGPERS2034class_N, Histoplasma_capsulatum_HISSPFGMAR2044class_N, Histoplasma_capsulatum_HISSPFGLIN2055class_N, Histoplasma_capsulatum_HISSPFGJOS2044class_N, Histoplasma_capsulatum_HISSPFGGRE2022class_N, Histoplasma_capsulatum_HISSPFGFIN2028class_N, Histoplasma_capsulatum_HISSPFGFER2036class_N, Histoplasma_capsulatum_HISSPFGFAR0189class_N, Histoplasma_capsulatum_HISSPFGFAN2059class_N, Histoplasma_capsulatum_HISSPFGBON2001class_N, Histoplasma_capsulatum_HISSPFGBIK2051class_N, Histoplasma_capsulatum_HISSPFGAMA2041class_N, Histoplasma_capsulatum_HISSPCM6408class_N, Histoplasma_capsulatum_HISSPCM6015class_N, Histoplasma_capsulatum_HISSPB05821class_N, Histoplasma_capsulatum_HISSP11571Belem1class_N, Histoplasma_capsulatum_HISSP1014Belem3class_N, Histoplasma_capsulatum_07_12RJclass_N, Histoplasma_capsulatum_104_p_06class_N, Histoplasma_capsulatum_104_P_19class_N, Histoplasma_capsulatum_117_p_12class_N, Histoplasma_capsulatum_122_p_10_Bclass_N, Histoplasma_capsulatum_136_P_07class_N, Histoplasma_capsulatum_144_p_08class_N, Histoplasma_capsulatum_1517_p_17class_N, Histoplasma_capsulatum_256_P_18class_N, Histoplasma_capsulatum_CB063class_N, Histoplasma_capsulatum_CB066class_N, Histoplasma_capsulatum_CB180class_N, Histoplasma_capsulatum_CB0522class_N, Histoplasma_capsulatum_CB0532class_N, Histoplasma_capsulatum_CB0552class_N, Histoplasma_capsulatum_CB0622class_N, Histoplasma_capsulatum_CB0632class_N, Histoplasma_capsulatum_CB0642class_N, Histoplasma_capsulatum_CB0652class_N, Histoplasma_capsulatum_CB0662class_N, Histoplasma_capsulatum_CB1682class_N, Histoplasma_capsulatum_CB1742class_N, Histoplasma_capsulatum_CB1802class_N, Histoplasma_capsulatum_CB1862class_N, Histoplasma_capsulatum_CB1922class_N, Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class_N, Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class_N, Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class_N, Histoplasma_mississippienseII_SECH_102Nam2_505class_N, Histoplasma_mississippienseII_SECH_81Nam1_WU24class_N, Histoplasma_mississippienseII_SECH_83Nam1_CI_7class_N, Histoplasma_mississippienseII_SECH_84Nam1_CI_19class_N, Histoplasma_mississippienseII_SECH_85Nam1_CI_22class_N, Histoplasma_mississippienseII_SECH_86Nam1_CI_24class_N, Histoplasma_mississippienseII_SECH_87Nam1_CI_42class_N, Histoplasma_mississippienseII_SECH_88Nam1_CI_43class_N),by="class", type="full")
Total_table2<-join_all(list(Histoplasma_mississippienseII_SECH_89class_N,Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass_N, Histoplasma_mississippienseII_HCWU24class_N, Histoplasma_mississippienseII_HCCI_43class_N, Histoplasma_mississippienseII_HCCI_19class_N, Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass_N, Histoplasma_suramericanum_SECH_104Nam2_27_14class_N, Histoplasma_suramericanum_SECH_105Nam2_21_14class_N, Histoplasma_suramericanum_HC27_14class_N, Histoplasma_suramericanum_HC21_14class_N, Histoplasma_capsulatum_CB053class_N, Histoplasma_capsulatum_CB062class_N, Histoplasma_capsulatum_316_p_10class_N, Histoplasma_capsulatum_327_P_12class_N, Histoplasma_capsulatum_343_p_18class_N, Histoplasma_capsulatum_388_p_11class_N, Histoplasma_capsulatum_ES2_83Zclass_N,Histoplasma_capsulatum_ES2_85Zclass_N,Histoplasma_capsulatum_ES2_86Zclass_N,Histoplasma_capsulatum_ES2_87class_N,Histoplasma_capsulatum_ES2_88Zclass_N,Histoplasma_capsulatum_ES2_89Zclass_N,Histoplasma_capsulatum_ES2_90Zclass_N, Histoplasma_capsulatum_SA15class_N),by="class", type="full")
Table_ALL<-join_all(list(Total_table, Total_table2), by="class", type="full")
write_tsv(Table_ALL,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/All_Histo.tsv")
d <- as.data.frame(Table_ALL)
#convert NAs to 0
d[is.na(d)] <- 0

#Manipulated the table outputs to seven types of TEs instead of the 34 initially described, then loaded these back in.
Percentage_new5 <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Percentage_REDO.tsv")

##Building the stacked barplot
#Found this description to follow https://zebrabi.com/guide/how-to-customize-stacked-bar-chart-in-r-ggplot2/#:~:text=Creating%20a%20basic%20stacked%20bar%20chart%20using%20ggplot2%20in%20R&text=For%20instance%2C%20we%20can%20change,bars%20using%20the%20width%20parameter

library(dplyr)
library(tidyr)
library(stringr)
library(reshape2) 
library(reshape)
library(RColorBrewer)
library(ggplot2)

Percentage_new5<- as.data.frame(Percentage_new5)
Percentage_new5[is.na(Percentage_new5)] <- 0

#Used this tutorial to exlpain melt https://cran.r-project.org/web/packages/data.table/vignettes/datatable-reshape.html
#list the name of your species as seen in your table, include spaces
PERC_SUM_melt<- melt(Percentage_new5, id.vars="class", 
                 measure.vars=c("Histoplasma ohiense SECH_82-Nam1_CI_4","Histoplasma ohiense SECH_92-Nam2_G222B",
                                "Histoplasma ohiense SECH_93.Nam2_CI_6","Histoplasma ohiense SECH_94-Nam2_CI_9",
                                "Histoplasma ohiense SECH_95.Nam2_CI_10","Histoplasma ohiense SECH_96.Nam2_CI_17",
                                "Histoplasma ohiense SECH_97-Nam2_CI_18","Histoplasma ohiense SECH_98.Nam2_CI_30",
                                "Histoplasma ohiense SECH_99.Nam2_CI_35","Histoplasma capsulatum 19VMG15",
                                "Histoplasma capsulatum HISSP-CM7256","Histoplasma capsulatum Histo-485P20",
                                "Histoplasma capsulatum JB_01752-Hc_01752","Histoplasma capsulatum JB_021091-Hc_021091",
                                "Histoplasma capsulatum JB_031837-Hc_031837","Histoplasma capsulatum JB_042430-Hc_042430",
                                "Histoplasma capsulatum JB_062632-Hc_062632","Histoplasma capsulatum JB_062775-Hc_062775",
                                "Histoplasma capsulatum JB_073129-Hc_073129","Histoplasma capsulatum JB_083285_2-Hc_083285_2",
                                "Histoplasma capsulatum SECH_109","Histoplasma capsulatum SECH_110",
                                "Histoplasma capsulatum HCH143","Histoplasma capsulatum NACVFR_Histo_HC1070058",
                                "Histoplasma capsulatum HISSP-FGTRO0285","Histoplasma capsulatum HISSP-FGPSO2043",
                                "Histoplasma capsulatum HISSP-FGPIE2055","Histoplasma capsulatum HISSP-FGPIA2052",
                                "Histoplasma capsulatum HISSP-FGPERS2034","Histoplasma capsulatum HISSP-FGMAR2044",
                                "Histoplasma capsulatum HISSP-FGLIN2055","Histoplasma capsulatum HISSP-FGJOS2044",
                                "Histoplasma capsulatum HISSP-FGFIN2028","Histoplasma capsulatum HISSP-FGFER2036",
                                "Histoplasma capsulatum HISSP-FGFAR0189","Histoplasma capsulatum HISSP-FGFAN2059",
                                "Histoplasma capsulatum HISSP-FGBON2001","Histoplasma capsulatum HISSP-FGBIK2051",
                                "Histoplasma capsulatum HISSP-FGAMA2041","Histoplasma capsulatum HISSP-CM6408",
                                "Histoplasma capsulatum HISSP-CM6015","Histoplasma capsulatum HISSP-B05821",
                                "Histoplasma capsulatum HISSP-11571-Belem1","Histoplasma capsulatum HISSP-1014-Belem3",
                                "Histoplasma capsulatum 07_12-RJ","Histoplasma capsulatum 104_p_06",
                                "Histoplasma capsulatum 117_p_12","Histoplasma capsulatum 122_p_10_B",
                                "Histoplasma capsulatum 136_P_07","Histoplasma capsulatum 144_p_08",
                                "Histoplasma capsulatum 1517_p_17","Histoplasma capsulatum 256_P_18",
                                "Histoplasma capsulatum CB053","Histoplasma capsulatum CB055","Histoplasma capsulatum CB062",
                                "Histoplasma capsulatum CB063","Histoplasma capsulatum CB064","Histoplasma capsulatum CB065",
                                "Histoplasma capsulatum CB066","Histoplasma capsulatum CB168","Histoplasma capsulatum CB174",
                                "Histoplasma capsulatum CB180","Histoplasma capsulatum CB186","Histoplasma capsulatum CB192",
                                "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11",
                                "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14",
                                "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16",
                                "Histoplasma mississippiense II SECH_102_Nam2_505",
                                "Histoplasma mississippiense II SECH_83-Nam1_CI_7",
                                "Histoplasma mississippiense II SECH_84-Nam1_CI_19",
                                "Histoplasma mississippiense II SECH_85-Nam1_CI_22",
                                "Histoplasma mississippiense II SECH_86-Nam1_CI_24",
                                "Histoplasma mississippiense II SECH_87-Nam1_CI_42",
                                "Histoplasma mississippiense II SECH_88-Nam1_CI_43",
                                "Histoplasma mississippiense II SECH_89-Nam1_UCLA-531",
                                "Histoplasma mississippiense II SECH_90-Nam1_DOWNS",
                                "Histoplasma suramericanum SECH_103-Nam2_3_11G","Histoplasma suramericanum SECH_104-Nam2_27_14",
                                "Histoplasma capsulatum 316_p_10", "Histoplasma suramericanum SECH_105-Nam2_21_14",
                                "Histoplasma capsulatum 327_P_12","Histoplasma capsulatum 343_p_18",
                                "Histoplasma capsulatum 388_p_11","Histoplasma capsulatum ES2_83Z",
                                "Histoplasma capsulatum ES2_85Z","Histoplasma capsulatum ES2_86Z",
                                "Histoplasma capsulatum ES2_88Z","Histoplasma capsulatum ES2_89Z",
                                "Histoplasma capsulatum ES2_90Z","Histoplasma capsulatum SA15","Aspergillus fumigatus Af293",
                                "Histoplasma mississippiense II SECH_81-Nam1_WU24","Histoplasma ohiense SECH_91-Nam2_G217B",
                                "Histoplasma capsulatum SECH_107-mis_Hc_duboisii-B",
                                "Histoplasma capsulatum SECH_101-Nam2_G184A", "Histoplasma capsulatum SECH_100-Nam2_G186A"))


#Saving results to a PDF file
pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Stacked_barplot_REDO.pdf', width=15, height=20)

ggplot(PERC_SUM_melt, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of TE per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired")

dev.off()

###
#Add table to the phylotree
library(ggplot2)
library(ggtree)
library(dplyr)
library(ape)
library(phytools)
library(geiger)
library(phangorn)
library(RRphylo)

#Ran this analysis a couple of times, latest iteration had CDS and AA trees. 
#CDS TREE
tree<-read.tree("/nas/longleaf/home/taniak/taniak/Tree_Construction_PHYling/Trial4.CDS/combined_top50.partition.ULTRA.contree.nw")
tree$tip.label

#Making sure to indicate my root species
trees <- root(tree, outgroup="Aspergillus_fumigatus_Af293", resolve.root=TRUE)
trees

#Tree file needed underscores, below script to adjust to how table is represented.
#Rename the species in the tree
new_seq <- tree$tip.label
dd = data.frame(new_seq)
dd$Species = 0
dd$Species[grep("Af293", dd$new_seq)] = "Aspergillus fumigatus Af293"
dd$Species[grep("07_12-RJ", dd$new_seq)] = "Histoplasma capsulatum 07_12-RJ"
dd$Species[grep("HISSP-11571-Belem1", dd$new_seq)] = "Histoplasma capsulatum HISSP-11571-Belem1"
dd$Species[grep("HISSP-FGFER2036", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGFER2036"
dd$Species[grep("HISSP-CM6015", dd$new_seq)] = "Histoplasma capsulatum HISSP-CM6015"
dd$Species[grep("HISSP-FGAMA2041", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGAMA2041"
dd$Species[grep("HISSP-FGLIN2055", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGLIN2055"
dd$Species[grep("HISSP-FGFAN2059", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGFAN2059"
dd$Species[grep("HISSP-FGPIA2052", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGPIA2052"
dd$Species[grep("HISSP-FGPIE2055", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGPIE2055"
dd$Species[grep("HISSP-FGBON2001", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGBON2001"
dd$Species[grep("HISSP-FGBIK2051", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGBIK2051"
dd$Species[grep("HISSP-FGJOS2044", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGJOS2044"
dd$Species[grep("JB_062632-Hc_062632", dd$new_seq)] = "Histoplasma capsulatum JB_062632-Hc_062632" 
dd$Species[grep("NACVFR_Histo_HC1070058_2", dd$new_seq)] = "Histoplasma capsulatum NACVFR_Histo_HC1070058"
dd$Species[grep("SECH_103-Nam2_3_11G", dd$new_seq)] = "Histoplasma suramericanum SECH_103-Nam2_3_11G"
dd$Species[grep("HISSP-FGPERS2034", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGPERS2034"
dd$Species[grep("HISSP-B05821", dd$new_seq)] = "Histoplasma capsulatum HISSP-B05821"
dd$Species[grep("JB_031837-Hc_031837", dd$new_seq)] = "Histoplasma capsulatum JB_031837-Hc_031837"
dd$Species[grep("JB_042430-Hc_042430", dd$new_seq)] = "Histoplasma capsulatum JB_042430-Hc_042430"
dd$Species[grep("JB_062775-Hc_062775", dd$new_seq)] = "Histoplasma capsulatum JB_062775-Hc_062775"
dd$Species[grep("JB_083285_2-Hc_083285_2", dd$new_seq)] = "Histoplasma capsulatum JB_083285_2-Hc_083285_2"
dd$Species[grep("CB053-2", dd$new_seq)] = "Histoplasma capsulatum CB053" 
dd$Species[grep("CB055-2", dd$new_seq)] = "Histoplasma capsulatum CB055"
dd$Species[grep("G184A", dd$new_seq)] = "Histoplasma capsulatum SECH_101-Nam2_G184A"
dd$Species[grep("HISSP-CM6408", dd$new_seq)] = "Histoplasma capsulatum HISSP-CM6408"
dd$Species[grep("HISSP-1014-Belem3", dd$new_seq)] = "Histoplasma capsulatum HISSP-1014-Belem3"
dd$Species[grep("ES2_83Z", dd$new_seq)] = "Histoplasma capsulatum ES2_83Z"
dd$Species[grep("ES2_90Z", dd$new_seq)] = "Histoplasma capsulatum ES2_90Z"
dd$Species[grep("ES2_86Z", dd$new_seq)] = "Histoplasma capsulatum ES2_86Z"
dd$Species[grep("ES2_88Z", dd$new_seq)] = "Histoplasma capsulatum ES2_88Z"
dd$Species[grep("HCH143", dd$new_seq)] = "Histoplasma capsulatum HCH143"
dd$Species[grep("duboisii", dd$new_seq)] = "Histoplasma capsulatum SECH_107-mis_Hc_duboisii-B" 
dd$Species[grep("ES2_85Z", dd$new_seq)] = "Histoplasma capsulatum ES2_85Z"
dd$Species[grep("HISSP-CM7256", dd$new_seq)] = "Histoplasma capsulatum HISSP-CM7256"
dd$Species[grep("ES2_89Z", dd$new_seq)] = "Histoplasma capsulatum ES2_89Z"
dd$Species[grep("SA15", dd$new_seq)] = "Histoplasma capsulatum SA15"
dd$Species[grep("SECH_98-Nam2_CI_30", dd$new_seq)] = "Histoplasma ohiense SECH_98.Nam2_CI_30"
dd$Species[grep("SECH_110", dd$new_seq)] = "Histoplasma capsulatum SECH_110"
dd$Species[grep("G217B", dd$new_seq)] = "Histoplasma ohiense SECH_91-Nam2_G217B"
dd$Species[grep("SECH_92-Nam2_G222B", dd$new_seq)] = "Histoplasma ohiense SECH_92-Nam2_G222B"
dd$Species[grep("Hc1986", dd$new_seq)] = "Histoplasma ohiense Hc1986"
dd$Species[grep("SECH_82-Nam1_CI_4", dd$new_seq)] = "Histoplasma ohiense SECH_82-Nam1_CI_4"
dd$Species[grep("SECH_95-Nam2_CI_10", dd$new_seq)] = "Histoplasma ohiense SECH_95.Nam2_CI_10"
dd$Species[grep("SECH_94-Nam2_CI_9", dd$new_seq)] = "Histoplasma ohiense SECH_94-Nam2_CI_9"
dd$Species[grep("SECH_96-Nam2_CI_17", dd$new_seq)] = "Histoplasma ohiense SECH_96.Nam2_CI_17"
dd$Species[grep("SECH_93-Nam2_CI_6", dd$new_seq)] = "Histoplasma ohiense SECH_93.Nam2_CI_6"
dd$Species[grep("SECH_97-Nam2_CI_18", dd$new_seq)] = "Histoplasma ohiense SECH_97-Nam2_CI_18"
dd$Species[grep("SECH_99-Nam2_CI_35", dd$new_seq)] = "Histoplasma ohiense SECH_99.Nam2_CI_35"
dd$Species[grep("SECH_104-Nam2_27_14", dd$new_seq)] = "Histoplasma suramericanum SECH_104-Nam2_27_14" 
dd$Species[grep("CB174-2", dd$new_seq)] = "Histoplasma capsulatum CB174"
dd$Species[grep("SECH_109", dd$new_seq)] = "Histoplasma capsulatum SECH_109"
dd$Species[grep("CB062-2", dd$new_seq)] = "Histoplasma capsulatum CB062"
dd$Species[grep("CB065-2", dd$new_seq)] = "Histoplasma capsulatum CB065"
dd$Species[grep("CB192-2", dd$new_seq)] = "Histoplasma capsulatum CB192"
dd$Species[grep("CB064-2", dd$new_seq)] = "Histoplasma capsulatum CB064"
dd$Species[grep("CB180-2", dd$new_seq)] = "Histoplasma capsulatum CB180"
dd$Species[grep("CB063-2", dd$new_seq)] = "Histoplasma capsulatum CB063"
dd$Species[grep("CB186-2", dd$new_seq)] = "Histoplasma capsulatum CB186" 
dd$Species[grep("SECH_81", dd$new_seq)] = "Histoplasma mississippiense II SECH_81-Nam1_WU24"
dd$Species[grep("CB066-2", dd$new_seq)] = "Histoplasma capsulatum CB066"
dd$Species[grep("SECH_88-Nam1_CI_43", dd$new_seq)] = "Histoplasma mississippiense II SECH_88-Nam1_CI_43"
dd$Species[grep("SECH_86-Nam1_CI_24", dd$new_seq)] = "Histoplasma mississippiense II SECH_86-Nam1_CI_24"
dd$Species[grep("CB168-2", dd$new_seq)] = "Histoplasma capsulatum CB168"
dd$Species[grep("SECH_85-Nam1_CI_22", dd$new_seq)] = "Histoplasma mississippiense II SECH_85-Nam1_CI_22"
dd$Species[grep("SECH_83-Nam1_CI_7", dd$new_seq)] = "Histoplasma mississippiense II SECH_83-Nam1_CI_7" 
dd$Species[grep("SECH_87-Nam1_CI_42", dd$new_seq)] = "Histoplasma mississippiense II SECH_87-Nam1_CI_42"
dd$Species[grep("SECH_90-Nam1_DOWNS", dd$new_seq)] = "Histoplasma mississippiense II SECH_90-Nam1_DOWNS"
dd$Species[grep("SECH_102", dd$new_seq)] = "Histoplasma mississippiense II SECH_102_Nam2_505"
dd$Species[grep("SECH_100", dd$new_seq)] = "Histoplasma capsulatum SECH_100-Nam2_G186A"
dd$Species[grep("SECH_84-Nam1_CI_19", dd$new_seq)] = "Histoplasma mississippiense II SECH_84-Nam1_CI_19"
dd$Species[grep("SECH_89-Nam1_UCLA-531", dd$new_seq)] = "Histoplasma mississippiense II SECH_89-Nam1_UCLA-531"
dd$Species[grep("JB_01752-Hc_01752", dd$new_seq)] = "Histoplasma capsulatum JB_01752-Hc_01752"
dd$Species[grep("104_p_06", dd$new_seq)] = "Histoplasma capsulatum 104_p_06"
dd$Species[grep("1517_p_17", dd$new_seq)] = "Histoplasma capsulatum 1517_p_17"
dd$Species[grep("117_p_12", dd$new_seq)] = "Histoplasma capsulatum 117_p_12"
dd$Species[grep("388_p_11", dd$new_seq)] = "Histoplasma capsulatum 388_p_11" 
dd$Species[grep("144_p_08", dd$new_seq)] = "Histoplasma capsulatum 144_p_08"
dd$Species[grep("327_P_12", dd$new_seq)] = "Histoplasma capsulatum 327_P_12"
dd$Species[grep("316_p_10", dd$new_seq)] = "Histoplasma capsulatum 316_p_10"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S14", dd$new_seq)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14"
dd$Species[grep("Histo-485P20", dd$new_seq)] = "Histoplasma capsulatum Histo-485P20"
dd$Species[grep("136_P_07", dd$new_seq)] = "Histoplasma capsulatum 136_P_07"
dd$Species[grep("343_p_18", dd$new_seq)] = "Histoplasma capsulatum 343_p_18"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S11", dd$new_seq)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S16", dd$new_seq)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16"
dd$Species[grep("122_p_10_B", dd$new_seq)] = "Histoplasma capsulatum 122_p_10_B"
dd$Species[grep("256_P_18", dd$new_seq)] = "Histoplasma capsulatum 256_P_18" 
dd$Species[grep("JB_021091-Hc_021091", dd$new_seq)] = "Histoplasma capsulatum JB_021091-Hc_021091"
dd$Species[grep("JB_073129-Hc_073129", dd$new_seq)] = "Histoplasma capsulatum JB_073129-Hc_073129"
dd$Species[grep("HISSP-FGFAR0189", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGFAR0189"
dd$Species[grep("HISSP-FGFIN2028", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGFIN2028"
dd$Species[grep("HISSP-FGPSO2043", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGPSO2043"
dd$Species[grep("HISSP-FGMAR2044", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGMAR2044"
dd$Species[grep("HISSP-FGTRO0285", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGTRO0285"
tree$tip.label <- dd$Species

# Visualize the tree
p <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 3, align=TRUE) +
  theme_tree2() +
  theme(legend.position = "none")

print(p)

# Save the tree as a PDF
ggsave("species_tree_with_node_numbers.pdf", p, width = 20, height = 5, dpi = 150)
ggsave(filename = "species_tree_with_node_numbers.png", plot = p,width = 9, height = 6, units = c("in"), dpi = 300, scale =2, device = "png")

#YOU NEED TO REARRANGE THE COLUMNS SO PROGRAM CAN READ IT PROPERLY
colnames(PERC_SUM_melt)
PERC_SUM_melt2 <- PERC_SUM_melt[,c(2,1,3)]

#Combining Tree with Stacked barplot results
#pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_Stacked_barplot_REDO.pdf', width=15, height=15)

p2 <- p + geom_facet(panel = "TE Abundance", data = PERC_SUM_melt2, geom = geom_col, 
           aes(x = value, color = class, 
               fill = class), orientation = 'y', width = .6) + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.05) +
  theme_tree2(legend.position=c(.95, .25)) + xlim_expand(c(0,0.8), 'TE Abundance') + xlim_expand(c(0, 3), 'Tree')

#dev.off()

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_Stacked_barplot_new5-2.pdf', width=15, height=15)

facet_widths(p2, widths = c(1, 0.40))

dev.off()


############~AMINO ACID TREE~###################################

##Building the stacked barplot
#Found this description to follow https://zebrabi.com/guide/how-to-customize-stacked-bar-chart-in-r-ggplot2/#:~:text=Creating%20a%20basic%20stacked%20bar%20chart%20using%20ggplot2%20in%20R&text=For%20instance%2C%20we%20can%20change,bars%20using%20the%20width%20parameter

#Add table to the phylotree
library(ggplot2)
library(ggtree)
library(dplyr)
library(ape)
library(phytools)
# AA TREE
tree_AA<-read.nexus("/nas/longleaf/home/taniak/taniak/Tree_Construction_PHYling/Trial4.AA/combined_top50.partition.ULTRA.contree.nw")
tree_AA
print(tree_AA)
ggtree(tree_AA)
tree_AA <- root(tree_AA, outgroup="Aspergillus_fumigatus_Af293", resolve.root=TRUE)

#Rename the species in the tree
new_seq2 <- tree_AA$tip.label
dd = data.frame(new_seq2)
dd$Species = 0
dd$Species[grep("Af293", dd$new_seq2)] = "Aspergillus fumigatus Af293"
dd$Species[grep("19VMG-15", dd$new_seq2)] = "Histoplasma capsulatum 19VMG15"
dd$Species[grep("07_12-RJ", dd$new_seq2)] = "Histoplasma capsulatum 07_12-RJ"
dd$Species[grep("HISSP-11571-Belem1", dd$new_seq2)] = "Histoplasma capsulatum HISSP-11571-Belem1"
dd$Species[grep("HISSP-FGFER2036", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGFER2036"
dd$Species[grep("HISSP-CM6015", dd$new_seq2)] = "Histoplasma capsulatum HISSP-CM6015"
dd$Species[grep("HISSP-FGAMA2041", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGAMA2041"
dd$Species[grep("HISSP-FGLIN2055", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGLIN2055"
dd$Species[grep("HISSP-FGFAN2059", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGFAN2059"
dd$Species[grep("HISSP-FGPIA2052", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGPIA2052"
dd$Species[grep("HISSP-FGPIE2055", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGPIE2055"
dd$Species[grep("HISSP-FGBON2001", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGBON2001"
dd$Species[grep("HISSP-FGBIK2051", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGBIK2051"
dd$Species[grep("HISSP-FGJOS2044", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGJOS2044"
dd$Species[grep("JB_062632-Hc_062632", dd$new_seq2)] = "Histoplasma capsulatum JB_062632-Hc_062632" 
dd$Species[grep("NACVFR_Histo_HC1070058_2", dd$new_seq2)] = "Histoplasma capsulatum NACVFR_Histo_HC1070058"
dd$Species[grep("SECH_103-Nam2_3_11G", dd$new_seq2)] = "Histoplasma suramericanum SECH_103-Nam2_3_11G"
dd$Species[grep("SECH_105", dd$new_seq2)] = "Histoplasma suramericanum SECH_105-Nam2_21_14"
dd$Species[grep("HISSP-FGPERS2034", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGPERS2034"
dd$Species[grep("HISSP-B05821", dd$new_seq2)] = "Histoplasma capsulatum HISSP-B05821"
dd$Species[grep("JB_031837-Hc_031837", dd$new_seq2)] = "Histoplasma capsulatum JB_031837-Hc_031837"
dd$Species[grep("JB_042430-Hc_042430", dd$new_seq2)] = "Histoplasma capsulatum JB_042430-Hc_042430"
dd$Species[grep("JB_062775-Hc_062775", dd$new_seq2)] = "Histoplasma capsulatum JB_062775-Hc_062775"
dd$Species[grep("JB_083285_2-Hc_083285_2", dd$new_seq2)] = "Histoplasma capsulatum JB_083285_2-Hc_083285_2"
dd$Species[grep("CB053-2", dd$new_seq2)] = "Histoplasma capsulatum CB053" 
dd$Species[grep("CB055-2", dd$new_seq2)] = "Histoplasma capsulatum CB055"
dd$Species[grep("G184A", dd$new_seq2)] = "Histoplasma capsulatum SECH_101-Nam2_G184A"
dd$Species[grep("HISSP-CM6408", dd$new_seq2)] = "Histoplasma capsulatum HISSP-CM6408"
dd$Species[grep("HISSP-1014-Belem3", dd$new_seq2)] = "Histoplasma capsulatum HISSP-1014-Belem3"
dd$Species[grep("ES2_83Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_83Z"
dd$Species[grep("ES2_90Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_90Z"
dd$Species[grep("ES2_86Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_86Z"
dd$Species[grep("ES2_88Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_88Z"
dd$Species[grep("HCH143", dd$new_seq2)] = "Histoplasma capsulatum HCH143"
dd$Species[grep("duboisii", dd$new_seq2)] = "Histoplasma capsulatum SECH_107-mis_Hc_duboisii-B" 
dd$Species[grep("ES2_85Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_85Z"
dd$Species[grep("HISSP-CM7256", dd$new_seq2)] = "Histoplasma capsulatum HISSP-CM7256"
dd$Species[grep("ES2_89Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_89Z"
dd$Species[grep("SA15", dd$new_seq2)] = "Histoplasma capsulatum SA15"
dd$Species[grep("Histoplasma_ohiense_SECH_98.Nam2_CI_30", dd$new_seq2)] = "Histoplasma ohiense SECH_98_Nam2_CI_30"
dd$Species[grep("SECH_110", dd$new_seq2)] = "Histoplasma capsulatum SECH_110"
dd$Species[grep("G217B", dd$new_seq2)] = "Histoplasma ohiense SECH_91-Nam2_G217B"
dd$Species[grep("SECH_92-Nam2_G222B", dd$new_seq2)] = "Histoplasma ohiense SECH_92-Nam2_G222B"
dd$Species[grep("Hc1986", dd$new_seq2)] = "Histoplasma ohiense Hc1986"
dd$Species[grep("SECH_82-Nam1_CI_4", dd$new_seq2)] = "Histoplasma ohiense SECH_82-Nam1_CI_4"
dd$Species[grep("Histoplasma_ohiense_SECH_95.Nam2_CI_10", dd$new_seq2)] = "Histoplasma ohiense SECH_95_Nam2_CI_10"
dd$Species[grep("SECH_94-Nam2_CI_9", dd$new_seq2)] = "Histoplasma ohiense SECH_94-Nam2_CI_9"
dd$Species[grep("Histoplasma_ohiense_SECH_96.Nam2_CI_17", dd$new_seq2)] = "Histoplasma ohiense SECH_96_Nam2_CI_17"
dd$Species[grep("Histoplasma_ohiense_SECH_93.Nam2_CI_6", dd$new_seq2)] = "Histoplasma ohiense SECH_93_Nam2_CI_6"
dd$Species[grep("SECH_97-Nam2_CI_18", dd$new_seq2)] = "Histoplasma ohiense SECH_97-Nam2_CI_18"
dd$Species[grep("SECH_99-Nam2_CI_35", dd$new_seq2)] = "Histoplasma ohiense SECH_99_Nam2_CI_35"
dd$Species[grep("SECH_104-Nam2_27_14", dd$new_seq2)] = "Histoplasma suramericanum SECH_104-Nam2_27_14" 
dd$Species[grep("CB174-2", dd$new_seq2)] = "Histoplasma capsulatum CB174"
dd$Species[grep("SECH_109", dd$new_seq2)] = "Histoplasma capsulatum SECH_109"
dd$Species[grep("CB062-2", dd$new_seq2)] = "Histoplasma capsulatum CB062"
dd$Species[grep("CB065-2", dd$new_seq2)] = "Histoplasma capsulatum CB065"
dd$Species[grep("CB192-2", dd$new_seq2)] = "Histoplasma capsulatum CB192"
dd$Species[grep("CB064-2", dd$new_seq2)] = "Histoplasma capsulatum CB064"
dd$Species[grep("CB180-2", dd$new_seq2)] = "Histoplasma capsulatum CB180"
dd$Species[grep("CB063-2", dd$new_seq2)] = "Histoplasma capsulatum CB063"
dd$Species[grep("CB186-2", dd$new_seq2)] = "Histoplasma capsulatum CB186" 
dd$Species[grep("SECH_81", dd$new_seq2)] = "Histoplasma mississippiense II SECH_81-Nam1_WU24"
dd$Species[grep("CB066-2", dd$new_seq2)] = "Histoplasma capsulatum CB066"
dd$Species[grep("SECH_88-Nam1_CI_43", dd$new_seq2)] = "Histoplasma mississippiense II SECH_88-Nam1_CI_43"
dd$Species[grep("SECH_86-Nam1_CI_24", dd$new_seq2)] = "Histoplasma mississippiense II SECH_86-Nam1_CI_24"
dd$Species[grep("CB168-2", dd$new_seq2)] = "Histoplasma capsulatum CB168"
dd$Species[grep("SECH_85-Nam1_CI_22", dd$new_seq2)] = "Histoplasma mississippiense II SECH_85-Nam1_CI_22"
dd$Species[grep("SECH_83-Nam1_CI_7", dd$new_seq2)] = "Histoplasma mississippiense II SECH_83-Nam1_CI_7" 
dd$Species[grep("SECH_87-Nam1_CI_42", dd$new_seq2)] = "Histoplasma mississippiense II SECH_87-Nam1_CI_42"
dd$Species[grep("SECH_90-Nam1_DOWNS", dd$new_seq2)] = "Histoplasma mississippiense II SECH_90-Nam1_DOWNS"
dd$Species[grep("SECH_102", dd$new_seq2)] = "Histoplasma mississippiense II SECH_102_Nam2_505"
dd$Species[grep("SECH_100", dd$new_seq2)] = "Histoplasma capsulatum SECH_100-Nam2_G186A"
dd$Species[grep("SECH_84-Nam1_CI_19", dd$new_seq2)] = "Histoplasma mississippiense II SECH_84-Nam1_CI_19"
dd$Species[grep("SECH_89-Nam1_UCLA-531", dd$new_seq2)] = "Histoplasma mississippiense II SECH_89-Nam1_UCLA-531"
dd$Species[grep("JB_01752-Hc_01752", dd$new_seq2)] = "Histoplasma capsulatum JB_01752-Hc_01752"
dd$Species[grep("104_p_06", dd$new_seq2)] = "Histoplasma capsulatum 104_p_06"
dd$Species[grep("1517_p_17", dd$new_seq2)] = "Histoplasma capsulatum 1517_p_17"
dd$Species[grep("117_p_12", dd$new_seq2)] = "Histoplasma capsulatum 117_p_12"
dd$Species[grep("388_p_11", dd$new_seq2)] = "Histoplasma capsulatum 388_p_11" 
dd$Species[grep("144_p_08", dd$new_seq2)] = "Histoplasma capsulatum 144_p_08"
dd$Species[grep("327_P_12", dd$new_seq2)] = "Histoplasma capsulatum 327_P_12"
dd$Species[grep("316_p_10", dd$new_seq2)] = "Histoplasma capsulatum 316_p_10"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S14", dd$new_seq2)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14"
dd$Species[grep("Histo-485P20", dd$new_seq2)] = "Histoplasma capsulatum Histo-485P20"
dd$Species[grep("136_P_07", dd$new_seq2)] = "Histoplasma capsulatum 136_P_07"
dd$Species[grep("343_p_18", dd$new_seq2)] = "Histoplasma capsulatum 343_p_18"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S11", dd$new_seq2)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S16", dd$new_seq2)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16"
dd$Species[grep("122_p_10_B", dd$new_seq2)] = "Histoplasma capsulatum 122_p_10_B"
dd$Species[grep("256_P_18", dd$new_seq2)] = "Histoplasma capsulatum 256_P_18" 
dd$Species[grep("JB_021091-Hc_021091", dd$new_seq2)] = "Histoplasma capsulatum JB_021091-Hc_021091"
dd$Species[grep("JB_073129-Hc_073129", dd$new_seq2)] = "Histoplasma capsulatum JB_073129-Hc_073129"
dd$Species[grep("HISSP-FGFAR0189", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGFAR0189"
dd$Species[grep("HISSP-FGFIN2028", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGFIN2028"
dd$Species[grep("HISSP-FGPSO2043", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGPSO2043"
dd$Species[grep("HISSP-FGMAR2044", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGMAR2044"
dd$Species[grep("HISSP-FGTRO0285", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGTRO0285"
tree_AA$tip.label <- dd$Species

# Visualize the tree
p_AA <- ggtree(tree_AA, layout = "rectangular") +
  geom_tiplab(size = 3) +
  theme_tree2() + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.008) +
  theme(legend.position = "none")
print(p_AA)

# Save the tree as a PDF
ggsave("AA_ultrametric_tree.pdf", p_AA, width = 20, height = 5, dpi = 150)
ggsave(filename = "species_tree_with_numbers.png", plot = p,width = 9, height = 6, units = c("in"), dpi = 300, scale =2, device = "png")

p_AA <-  ggtree(tree_AA, layout = "rectangular") +
  theme_tree2() + geom_tiplab(align=TRUE)

print(p_AA)

#pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_Stacked_barplot_AA2.pdf', width=15, height=15)

p_AA1 <- p_AA + geom_facet(panel = "TE Abundance", data = PERC_SUM_melt2, geom = geom_col, 
                     aes(x = value, color = class, 
                         fill = class), orientation = 'y', width = .6) + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.05) +
  theme_tree2(legend.position=c(.95, .25)) + xlim_expand(c(0,0.8), 'TE Abundance') + xlim_expand(c(0, 2.5), 'Tree')

#dev.off()

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_Stacked_barplotAA2.pdf', width=15, height=15)

facet_widths(p_AA3, widths = c(1, 0.40))

dev.off()

######### Testing for phylogenetic signal ###########
#AA Tree
#tree_AA

colnames(Percentage_new5) = gsub("[.]", "_", colnames(Percentage_new5))
tree_AA$tip.label = gsub("[-]", "_", tree_AA$tip.label)
all(colnames(Percentage_new5)  %in% tree_AA$tip.label)
colnames(Percentage_new5)[which(colnames(Percentage_new5)  %in% tree_AA$tip.label==F)]
#Just first column doesn't agree which is class

#Test TE for phylogenetic signal - but phylosig requires named vector for data
row = which(Percentage_new5$class=="DNA TE")
data.vec = as.numeric(Percentage_new5[row, ]); names(data.vec) = colnames(Percentage_new5)
# remove the first column, which contains no data
data.vec=data.vec[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam = phylosig(tree=tree_AA, x=data.vec, method="lambda", test=TRUE)
sig.lam

#Phylogenetic signal lambda : 0.924593 
#logL(lambda) : 252.703 
#LR(lambda=0) : 34.2999 
#P-value (based on LR test) : 4.72403e-09

sig.k = (phylosig(tree=tree_AA, x=data.vec, method="K", test=TRUE))
sig.k

#Phylogenetic signal K : 0.780382 
#P-value (based on 1000 randomizations) : 0.017 


#Redo AA trees with Blastomyces_parvus_UAMH130 and Emmonsia_crescens_UAMH4076

##Emmonsia crescens UAMH4076
#Add table to the phylotree
library(ggplot2)
library(ggtree)
library(dplyr)
library(ape)
library(phytools)
# AA TREE
tree_AAE<-read.tree("/nas/longleaf/home/taniak/taniak/Tree_Construction_PHYling/Trial5.AA/combined_top50.partition.ULTRA.contree.nw")
tree_AAE
print(tree_AAE)
ggtree(tree_AAE)
tree_AAE <- root(tree_AAE, outgroup="Emmonsia_crescens_UAMH4076", resolve.root=TRUE)

Percentage_Emm <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Percentage_Emm.tsv")
Percentage_Emm<- as.data.frame(Percentage_Emm)
Percentage_Emm[is.na(Percentage_Emm)] <- 0

PERC_SUM_melt_Emm<- melt(Percentage_Emm, id.vars="class", 
                     measure.vars=c("Histoplasma capsulatum HISSP-CM7256","Histoplasma capsulatum Histo-485P20",
                                    "Histoplasma capsulatum JB_01752-Hc_01752","Histoplasma capsulatum JB_021091-Hc_021091",
                                    "Histoplasma capsulatum JB_031837-Hc_031837","Histoplasma capsulatum JB_042430-Hc_042430",
                                    "Histoplasma capsulatum JB_062632-Hc_062632","Histoplasma capsulatum JB_062775-Hc_062775",
                                    "Histoplasma capsulatum JB_073129-Hc_073129","Histoplasma capsulatum JB_083285_2-Hc_083285_2",
                                    "Histoplasma capsulatum SECH_101-Nam2_G184A","Histoplasma mississippiense II SECH_102_Nam2_505",
                                    "Histoplasma suramericanum SECH_103-Nam2_3_11G","Histoplasma suramericanum SECH_104-Nam2_27_14",
                                    "Histoplasma capsulatum SECH_107-mis_Hc_duboisii-B","Histoplasma capsulatum SECH_109",
                                    "Histoplasma capsulatum SECH_110","Histoplasma mississippiense II SECH_81-Nam1_WU24",
                                    "Histoplasma ohiense SECH_82-Nam1_CI_4","Histoplasma mississippiense II SECH_83-Nam1_CI_7",
                                    "Histoplasma mississippiense II SECH_84-Nam1_CI_19","Histoplasma mississippiense II SECH_85-Nam1_CI_22",
                                    "Histoplasma mississippiense II SECH_86-Nam1_CI_24","Histoplasma mississippiense II SECH_87-Nam1_CI_42",
                                    "Histoplasma mississippiense II SECH_88-Nam1_CI_43","Histoplasma mississippiense II SECH_89-Nam1_UCLA-531",
                                    "Histoplasma mississippiense II SECH_90-Nam1_DOWNS","Histoplasma ohiense SECH_91-Nam2_G217B",
                                    "Histoplasma ohiense SECH_92-Nam2_G222B","Histoplasma ohiense SECH_93.Nam2_CI_6","Histoplasma ohiense SECH_94-Nam2_CI_9",
                                    "Histoplasma ohiense SECH_95.Nam2_CI_10","Histoplasma ohiense SECH_96.Nam2_CI_17","Histoplasma ohiense SECH_97-Nam2_CI_18",
                                    "Histoplasma ohiense SECH_98.Nam2_CI_30","Histoplasma ohiense SECH_99.Nam2_CI_35","Histoplasma capsulatum HCH143",
                                    "Histoplasma capsulatum NACVFR_Histo_HC1070058","Histoplasma capsulatum HISSP-FGTRO0285",
                                    "Histoplasma capsulatum HISSP-FGPSO2043","Histoplasma capsulatum HISSP-FGPIE2055","Histoplasma capsulatum HISSP-FGPIA2052",
                                    "Histoplasma capsulatum HISSP-FGPERS2034","Histoplasma capsulatum HISSP-FGMAR2044","Histoplasma capsulatum HISSP-FGLIN2055",
                                    "Histoplasma capsulatum HISSP-FGJOS2044","Histoplasma capsulatum HISSP-FGFIN2028","Histoplasma capsulatum HISSP-FGFER2036",
                                    "Histoplasma capsulatum HISSP-FGFAR0189","Histoplasma capsulatum HISSP-FGFAN2059","Histoplasma capsulatum HISSP-FGBON2001",
                                    "Histoplasma capsulatum HISSP-FGBIK2051","Histoplasma capsulatum HISSP-FGAMA2041","Histoplasma capsulatum HISSP-CM6408",
                                    "Histoplasma capsulatum HISSP-CM6015","Histoplasma capsulatum HISSP-B05821","Histoplasma capsulatum HISSP-11571-Belem1",
                                    "Histoplasma capsulatum HISSP-1014-Belem3","Histoplasma capsulatum 07_12-RJ","Histoplasma capsulatum 104_p_06",
                                    "Histoplasma capsulatum 117_p_12","Histoplasma capsulatum 122_p_10_B","Histoplasma capsulatum 136_P_07",
                                    "Histoplasma capsulatum 144_p_08","Histoplasma capsulatum 1517_p_17","Histoplasma capsulatum 256_P_18",
                                    "Histoplasma capsulatum 316_p_10","Histoplasma capsulatum 327_P_12","Histoplasma capsulatum 343_p_18",
                                    "Histoplasma capsulatum 388_p_11","Histoplasma capsulatum ES2_83Z","Histoplasma capsulatum ES2_85Z",
                                    "Histoplasma capsulatum ES2_86Z","Histoplasma capsulatum ES2_88Z","Histoplasma capsulatum ES2_89Z",
                                    "Histoplasma capsulatum ES2_90Z","Histoplasma capsulatum SA15","Emmonsia crescens UAMH4076","Histoplasma capsulatum CB053",
                                    "Histoplasma capsulatum CB055","Histoplasma capsulatum CB062","Histoplasma capsulatum CB063","Histoplasma capsulatum CB064",
                                    "Histoplasma capsulatum CB065","Histoplasma capsulatum CB066","Histoplasma capsulatum CB168","Histoplasma capsulatum CB174",
                                    "Histoplasma capsulatum CB180","Histoplasma capsulatum CB186","Histoplasma capsulatum CB192",
                                    "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11","Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14",
                                    "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16","Histoplasma capsulatum SECH_100-Nam2_G186A",
                                    "Histoplasma capsulatum 19VMG15","Histoplasma suramericanum SECH_105-Nam2_21_14"))



new_seq2 <- tree_AAE$tip.label
dd = data.frame(new_seq2)
dd$Species = 0
dd$Species[grep("UAMH4076", dd$new_seq2)] = "Emmonsia crescens UAMH4076"
dd$Species[grep("19VMG-15", dd$new_seq2)] = "Histoplasma capsulatum 19VMG15"
dd$Species[grep("07_12-RJ", dd$new_seq2)] = "Histoplasma capsulatum 07_12-RJ"
dd$Species[grep("HISSP-11571-Belem1", dd$new_seq2)] = "Histoplasma capsulatum HISSP-11571-Belem1"
dd$Species[grep("HISSP-FGFER2036", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGFER2036"
dd$Species[grep("HISSP-CM6015", dd$new_seq2)] = "Histoplasma capsulatum HISSP-CM6015"
dd$Species[grep("HISSP-FGAMA2041", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGAMA2041"
dd$Species[grep("HISSP-FGLIN2055", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGLIN2055"
dd$Species[grep("HISSP-FGFAN2059", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGFAN2059"
dd$Species[grep("HISSP-FGPIA2052", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGPIA2052"
dd$Species[grep("HISSP-FGPIE2055", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGPIE2055"
dd$Species[grep("HISSP-FGBON2001", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGBON2001"
dd$Species[grep("HISSP-FGBIK2051", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGBIK2051"
dd$Species[grep("HISSP-FGJOS2044", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGJOS2044"
dd$Species[grep("JB_062632-Hc_062632", dd$new_seq2)] = "Histoplasma capsulatum JB_062632-Hc_062632" 
dd$Species[grep("NACVFR_Histo_HC1070058_2", dd$new_seq2)] = "Histoplasma capsulatum NACVFR_Histo_HC1070058"
dd$Species[grep("SECH_103-Nam2_3_11G", dd$new_seq2)] = "Histoplasma suramericanum SECH_103-Nam2_3_11G"
dd$Species[grep("SECH_105", dd$new_seq2)] = "Histoplasma suramericanum SECH_105-Nam2_21_14"
dd$Species[grep("HISSP-FGPERS2034", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGPERS2034"
dd$Species[grep("HISSP-B05821", dd$new_seq2)] = "Histoplasma capsulatum HISSP-B05821"
dd$Species[grep("JB_031837-Hc_031837", dd$new_seq2)] = "Histoplasma capsulatum JB_031837-Hc_031837"
dd$Species[grep("JB_042430-Hc_042430", dd$new_seq2)] = "Histoplasma capsulatum JB_042430-Hc_042430"
dd$Species[grep("JB_062775-Hc_062775", dd$new_seq2)] = "Histoplasma capsulatum JB_062775-Hc_062775"
dd$Species[grep("JB_083285_2-Hc_083285_2", dd$new_seq2)] = "Histoplasma capsulatum JB_083285_2-Hc_083285_2"
dd$Species[grep("CB053-2", dd$new_seq2)] = "Histoplasma capsulatum CB053" 
dd$Species[grep("CB055-2", dd$new_seq2)] = "Histoplasma capsulatum CB055"
dd$Species[grep("G184A", dd$new_seq2)] = "Histoplasma capsulatum SECH_101-Nam2_G184A"
dd$Species[grep("HISSP-CM6408", dd$new_seq2)] = "Histoplasma capsulatum HISSP-CM6408"
dd$Species[grep("HISSP-1014-Belem3", dd$new_seq2)] = "Histoplasma capsulatum HISSP-1014-Belem3"
dd$Species[grep("ES2_83Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_83Z"
dd$Species[grep("ES2_90Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_90Z"
dd$Species[grep("ES2_86Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_86Z"
dd$Species[grep("ES2_88Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_88Z"
dd$Species[grep("HCH143", dd$new_seq2)] = "Histoplasma capsulatum HCH143"
dd$Species[grep("duboisii", dd$new_seq2)] = "Histoplasma capsulatum SECH_107-mis_Hc_duboisii-B" 
dd$Species[grep("ES2_85Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_85Z"
dd$Species[grep("HISSP-CM7256", dd$new_seq2)] = "Histoplasma capsulatum HISSP-CM7256"
dd$Species[grep("ES2_89Z", dd$new_seq2)] = "Histoplasma capsulatum ES2_89Z"
dd$Species[grep("SA15", dd$new_seq2)] = "Histoplasma capsulatum SA15"
dd$Species[grep("Histoplasma_ohiense_SECH_98.Nam2_CI_30", dd$new_seq2)] = "Histoplasma ohiense SECH_98.Nam2_CI_30"
dd$Species[grep("SECH_110", dd$new_seq2)] = "Histoplasma capsulatum SECH_110"
dd$Species[grep("G217B", dd$new_seq2)] = "Histoplasma ohiense SECH_91-Nam2_G217B"
dd$Species[grep("SECH_92-Nam2_G222B", dd$new_seq2)] = "Histoplasma ohiense SECH_92-Nam2_G222B"
dd$Species[grep("Hc1986", dd$new_seq2)] = "Histoplasma ohiense Hc1986"
dd$Species[grep("SECH_82-Nam1_CI_4", dd$new_seq2)] = "Histoplasma ohiense SECH_82-Nam1_CI_4"
dd$Species[grep("Histoplasma_ohiense_SECH_95.Nam2_CI_10", dd$new_seq2)] = "Histoplasma ohiense SECH_95.Nam2_CI_10"
dd$Species[grep("SECH_94-Nam2_CI_9", dd$new_seq2)] = "Histoplasma ohiense SECH_94-Nam2_CI_9"
dd$Species[grep("Histoplasma_ohiense_SECH_96.Nam2_CI_17", dd$new_seq2)] = "Histoplasma ohiense SECH_96.Nam2_CI_17"
dd$Species[grep("Histoplasma_ohiense_SECH_93.Nam2_CI_6", dd$new_seq2)] = "Histoplasma ohiense SECH_93.Nam2_CI_6"
dd$Species[grep("SECH_97-Nam2_CI_18", dd$new_seq2)] = "Histoplasma ohiense SECH_97-Nam2_CI_18"
dd$Species[grep("SECH_99-Nam2_CI_35", dd$new_seq2)] = "Histoplasma ohiense SECH_99.Nam2_CI_35"
dd$Species[grep("SECH_104-Nam2_27_14", dd$new_seq2)] = "Histoplasma suramericanum SECH_104-Nam2_27_14" 
dd$Species[grep("CB174-2", dd$new_seq2)] = "Histoplasma capsulatum CB174"
dd$Species[grep("SECH_109", dd$new_seq2)] = "Histoplasma capsulatum SECH_109"
dd$Species[grep("CB062-2", dd$new_seq2)] = "Histoplasma capsulatum CB062"
dd$Species[grep("CB065-2", dd$new_seq2)] = "Histoplasma capsulatum CB065"
dd$Species[grep("CB192-2", dd$new_seq2)] = "Histoplasma capsulatum CB192"
dd$Species[grep("CB064-2", dd$new_seq2)] = "Histoplasma capsulatum CB064"
dd$Species[grep("CB180-2", dd$new_seq2)] = "Histoplasma capsulatum CB180"
dd$Species[grep("CB063-2", dd$new_seq2)] = "Histoplasma capsulatum CB063"
dd$Species[grep("CB186-2", dd$new_seq2)] = "Histoplasma capsulatum CB186" 
dd$Species[grep("SECH_81", dd$new_seq2)] = "Histoplasma mississippiense II SECH_81-Nam1_WU24"
dd$Species[grep("CB066-2", dd$new_seq2)] = "Histoplasma capsulatum CB066"
dd$Species[grep("SECH_88-Nam1_CI_43", dd$new_seq2)] = "Histoplasma mississippiense II SECH_88-Nam1_CI_43"
dd$Species[grep("SECH_86-Nam1_CI_24", dd$new_seq2)] = "Histoplasma mississippiense II SECH_86-Nam1_CI_24"
dd$Species[grep("CB168-2", dd$new_seq2)] = "Histoplasma capsulatum CB168"
dd$Species[grep("SECH_85-Nam1_CI_22", dd$new_seq2)] = "Histoplasma mississippiense II SECH_85-Nam1_CI_22"
dd$Species[grep("SECH_83-Nam1_CI_7", dd$new_seq2)] = "Histoplasma mississippiense II SECH_83-Nam1_CI_7" 
dd$Species[grep("SECH_87-Nam1_CI_42", dd$new_seq2)] = "Histoplasma mississippiense II SECH_87-Nam1_CI_42"
dd$Species[grep("SECH_90-Nam1_DOWNS", dd$new_seq2)] = "Histoplasma mississippiense II SECH_90-Nam1_DOWNS"
dd$Species[grep("SECH_102", dd$new_seq2)] = "Histoplasma mississippiense II SECH_102_Nam2_505"
dd$Species[grep("SECH_100", dd$new_seq2)] = "Histoplasma capsulatum SECH_100-Nam2_G186A"
dd$Species[grep("SECH_84-Nam1_CI_19", dd$new_seq2)] = "Histoplasma mississippiense II SECH_84-Nam1_CI_19"
dd$Species[grep("SECH_89-Nam1_UCLA-531", dd$new_seq2)] = "Histoplasma mississippiense II SECH_89-Nam1_UCLA-531"
dd$Species[grep("JB_01752-Hc_01752", dd$new_seq2)] = "Histoplasma capsulatum JB_01752-Hc_01752"
dd$Species[grep("104_p_06", dd$new_seq2)] = "Histoplasma capsulatum 104_p_06"
dd$Species[grep("1517_p_17", dd$new_seq2)] = "Histoplasma capsulatum 1517_p_17"
dd$Species[grep("117_p_12", dd$new_seq2)] = "Histoplasma capsulatum 117_p_12"
dd$Species[grep("388_p_11", dd$new_seq2)] = "Histoplasma capsulatum 388_p_11" 
dd$Species[grep("144_p_08", dd$new_seq2)] = "Histoplasma capsulatum 144_p_08"
dd$Species[grep("327_P_12", dd$new_seq2)] = "Histoplasma capsulatum 327_P_12"
dd$Species[grep("316_p_10", dd$new_seq2)] = "Histoplasma capsulatum 316_p_10"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S14", dd$new_seq2)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14"
dd$Species[grep("Histo-485P20", dd$new_seq2)] = "Histoplasma capsulatum Histo-485P20"
dd$Species[grep("136_P_07", dd$new_seq2)] = "Histoplasma capsulatum 136_P_07"
dd$Species[grep("343_p_18", dd$new_seq2)] = "Histoplasma capsulatum 343_p_18"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S11", dd$new_seq2)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11"
dd$Species[grep("Dr_Anuradha_Fungal_WGS_S16", dd$new_seq2)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16"
dd$Species[grep("122_p_10_B", dd$new_seq2)] = "Histoplasma capsulatum 122_p_10_B"
dd$Species[grep("256_P_18", dd$new_seq2)] = "Histoplasma capsulatum 256_P_18" 
dd$Species[grep("JB_021091-Hc_021091", dd$new_seq2)] = "Histoplasma capsulatum JB_021091-Hc_021091"
dd$Species[grep("JB_073129-Hc_073129", dd$new_seq2)] = "Histoplasma capsulatum JB_073129-Hc_073129"
dd$Species[grep("HISSP-FGFAR0189", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGFAR0189"
dd$Species[grep("HISSP-FGFIN2028", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGFIN2028"
dd$Species[grep("HISSP-FGPSO2043", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGPSO2043"
dd$Species[grep("HISSP-FGMAR2044", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGMAR2044"
dd$Species[grep("HISSP-FGTRO0285", dd$new_seq2)] = "Histoplasma capsulatum HISSP-FGTRO0285"
tree_AAE$tip.label <- dd$Species

# Visualize the tree
p_AAE <- ggtree(tree_AAE, layout = "rectangular") +
  geom_tiplab(size = 3) +
  theme_tree2() + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.008) +
  theme(legend.position = "none")
print(p_AAE)
# Save the tree as a PDF
ggsave("AAE_ultrametric_tree.pdf", p_AAE, width = 20, height = 5, dpi = 150)
ggsave(filename = "species_tree_with_numbersE.png", plot = p_AAE,width = 9, height = 6, units = c("in"), dpi = 300, scale =2, device = "png")

p_AAE <-  ggtree(tree_AAE, layout = "rectangular") +
  theme_tree2() + geom_tiplab(align=TRUE)

print(p_AAE)
#Sanity check!
colnames(Percentage_Emm)[which(colnames(Percentage_Emm)  %in% tree_AAE$tip.label==F)]

colnames(PERC_SUM_melt_Emm)
PERC_SUM_melt_Emm2 <- PERC_SUM_melt_Emm[,c(2,1,3)]

#pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_Stacked_barplot_AAE.pdf', width=15, height=15)

p_AAE1 <- p_AAE + geom_facet(panel = "TE Abundance", data = PERC_SUM_melt_Emm2, geom = geom_col, 
                           aes(x = value, color = class, 
                               fill = class), orientation = 'y', width = .6) + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.05) +
  theme_tree2(legend.position=c(.95, .25)) + xlim_expand(c(0,0.8), 'TE Abundance') + xlim_expand(c(0, 2.5), 'Tree')

#dev.off()

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_Stacked_barplotAAE.pdf', width=15, height=15)

facet_widths(p_AAE1, widths = c(1, 0.40))

dev.off()

#### BLASTOMYCES

#Add table to the phylotree
library(ggplot2)
library(ggtree)
library(dplyr)
library(ape)
library(phytools)
library(picante)
library(readr)
library(reshape2)

# AA TREE
tree_AAB<-read.tree("/nas/longleaf/home/taniak/taniak/Tree_Construction_PHYling/Trial5.AA-Blast/combined_top50.partition.ULTRA.contree.nw")
tree_AAB
is.ultrametric(tree_AAB)
ggtree(tree_AAB) + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.008)



#Percentage_B <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Percentage_BB.tsv")
Percentage_B <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Percentage_BBS.tsv")
Percentage_B<- as.data.frame(Percentage_B)
Percentage_B[is.na(Percentage_B)] <- 0

PERC_SUM_melt_B<- melt(Percentage_B, id.vars="class", 
                         measure.vars=c("Histoplasma capsulatum HISSP-CM7256","Histoplasma capsulatum Histo-485P20",
                                        "Histoplasma capsulatum JB_01752-Hc_01752","Histoplasma capsulatum JB_021091-Hc_021091",
                                        "Histoplasma capsulatum JB_031837-Hc_031837","Histoplasma capsulatum JB_042430-Hc_042430",
                                        "Histoplasma capsulatum JB_062632-Hc_062632","Histoplasma capsulatum JB_062775-Hc_062775",
                                        "Histoplasma capsulatum JB_073129-Hc_073129","Histoplasma capsulatum JB_083285_2-Hc_083285_2",
                                        "Histoplasma capsulatum senso_stricto SECH_101_H184AR","Histoplasma mississippiense SECH_102_Nam2_505",
                                        "Histoplasma suramericanum SECH_103_3_11G","Histoplasma suramericanum SECH_104_27_14",
                                        "Histoplasma capsulatum Africa SECH_107-duboisii-B-H88","Histoplasma capsulatum SECH_109",
                                        "Histoplasma capsulatum SECH_110","Histoplasma mississippiense SECH_81-Nam1_WU24",
                                        "Histoplasma ohiense SECH_82-Nam1_CI_4","Histoplasma mississippiense SECH_83-Nam1_CI_7",
                                        "Histoplasma mississippiense SECH_84-Nam1_CI_19","Histoplasma mississippiense SECH_85-Nam1_CI_22",
                                        "Histoplasma mississippiense SECH_86-Nam1_CI_24","Histoplasma mississippiense SECH_87-Nam1_CI_42",
                                        "Histoplasma mississippiense SECH_88-Nam1_CI_43","Histoplasma mississippiense SECH_89-Nam1_UCLA-531",
                                        "Histoplasma mississippiense SECH_90-Nam1_DOWNS","Histoplasma ohiense SECH_91-Nam2_G217B",
                                        "Histoplasma ohiense SECH_92-Nam2_G222B","Histoplasma ohiense SECH_93.Nam2_CI_6","Histoplasma ohiense SECH_94-Nam2_CI_9",
                                        "Histoplasma ohiense SECH_95.Nam2_CI_10","Histoplasma ohiense SECH_96.Nam2_CI_17","Histoplasma ohiense SECH_97-Nam2_CI_18",
                                        "Histoplasma ohiense SECH_98.Nam2_CI_30","Histoplasma ohiense SECH_99.Nam2_CI_35","Histoplasma capsulatum HCH143",
                                        "Histoplasma capsulatum NACVFR_Histo_HC1070058","Histoplasma capsulatum HISSP-FGTRO0285",
                                        "Histoplasma capsulatum HISSP-FGPSO2043","Histoplasma capsulatum HISSP-FGPIE2055","Histoplasma capsulatum HISSP-FGPIA2052",
                                        "Histoplasma capsulatum HISSP-FGPERS2034","Histoplasma capsulatum HISSP-FGMAR2044","Histoplasma capsulatum HISSP-FGLIN2055",
                                        "Histoplasma capsulatum HISSP-FGJOS2044","Histoplasma capsulatum HISSP-FGFIN2028","Histoplasma capsulatum HISSP-FGFER2036",
                                        "Histoplasma capsulatum HISSP-FGFAR0189","Histoplasma capsulatum HISSP-FGFAN2059","Histoplasma capsulatum HISSP-FGBON2001",
                                        "Histoplasma capsulatum HISSP-FGBIK2051","Histoplasma capsulatum HISSP-FGAMA2041","Histoplasma capsulatum HISSP-CM6408",
                                        "Histoplasma capsulatum HISSP-CM6015","Histoplasma capsulatum HISSP-B05821","Histoplasma capsulatum HISSP-11571-Belem1",
                                        "Histoplasma capsulatum HISSP-1014-Belem3","Histoplasma capsulatum 07_12-RJ","Histoplasma capsulatum 104_p_06",
                                        "Histoplasma capsulatum 117_p_12","Histoplasma capsulatum 122_p_10_B","Histoplasma capsulatum 136_P_07",
                                        "Histoplasma capsulatum 144_p_08","Histoplasma capsulatum 1517_p_17","Histoplasma capsulatum 256_P_18",
                                        "Histoplasma capsulatum 316_p_10","Histoplasma capsulatum 327_P_12","Histoplasma capsulatum 343_p_18",
                                        "Histoplasma capsulatum 388_p_11","Histoplasma capsulatum ES2_83Z","Histoplasma capsulatum ES2_85Z",
                                        "Histoplasma capsulatum ES2_86Z","Histoplasma capsulatum ES2_88Z","Histoplasma capsulatum ES2_89Z",
                                        "Histoplasma capsulatum ES2_90Z","Histoplasma capsulatum SA15","Blastomyces parvus UAMH130","Histoplasma capsulatum CB053",
                                        "Histoplasma capsulatum CB055","Histoplasma capsulatum CB062","Histoplasma capsulatum CB063","Histoplasma capsulatum CB064",
                                        "Histoplasma capsulatum CB065","Histoplasma capsulatum CB066","Histoplasma capsulatum CB168","Histoplasma capsulatum CB174",
                                        "Histoplasma capsulatum CB180","Histoplasma capsulatum CB186","Histoplasma capsulatum CB192",
                                        "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11","Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14",
                                        "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16","Histoplasma capsulatum senso_stricto SECH_100_H186AR",
                                        "Histoplasma capsulatum 19VMG15","Histoplasma suramericanum SECH_105_21_14"))


pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Stacked_barplot_Blast_Star.pdf', width=15, height=20)

ggplot(PERC_SUM_melt_B, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of TE per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired")

dev.off()

#Genome length + TE composition next to stacked barplot
#also need a more complex detailed stacked barplot
library(reshape2) 

Summ_GL_RE <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/SUM_Star_GL_TE.txt")
Summ_GL_RE_melt <- melt(Summ_GL_RE, id.vars="values", 
                           measure.vars=c("Histoplasma capsulatum 19VMG15","Histoplasma capsulatum HISSP-CM7256",
                                          "Histoplasma capsulatum Histo-485P20","Histoplasma capsulatum JB_01752-Hc_01752",
                                          "Histoplasma capsulatum JB_021091-Hc_021091","Histoplasma capsulatum JB_031837-Hc_031837",
                                          "Histoplasma capsulatum JB_042430-Hc_042430","Histoplasma capsulatum JB_062632-Hc_062632",
                                          "Histoplasma capsulatum JB_062775-Hc_062775","Histoplasma capsulatum JB_073129-Hc_073129",
                                          "Histoplasma capsulatum JB_083285_2-Hc_083285_2","Histoplasma capsulatum senso_stricto SECH_100_H186AR",
                                          "Histoplasma capsulatum senso_stricto SECH_101_H184AR","Histoplasma mississippiense SECH_102_Nam2_505",
                                          "Histoplasma suramericanum SECH_103_3_11G","Histoplasma suramericanum SECH_104_27_14",
                                          "Histoplasma suramericanum SECH_105_21_14","Histoplasma capsulatum Africa SECH_107-duboisii-B-H88",
                                          "Histoplasma capsulatum SECH_109","Histoplasma capsulatum SECH_110",
                                          "Histoplasma mississippiense SECH_81-Nam1_WU24","Histoplasma ohiense SECH_82-Nam1_CI_4",
                                          "Histoplasma mississippiense SECH_83-Nam1_CI_7","Histoplasma mississippiense SECH_84-Nam1_CI_19",
                                          "Histoplasma mississippiense SECH_85-Nam1_CI_22","Histoplasma mississippiense SECH_86-Nam1_CI_24",
                                          "Histoplasma mississippiense SECH_87-Nam1_CI_42","Histoplasma mississippiense SECH_88-Nam1_CI_43",
                                          "Histoplasma mississippiense SECH_89-Nam1_UCLA-531","Histoplasma mississippiense SECH_90-Nam1_DOWNS",
                                          "Histoplasma ohiense SECH_91-Nam2_G217B","Histoplasma ohiense SECH_92-Nam2_G222B",
                                          "Histoplasma ohiense SECH_93.Nam2_CI_6","Histoplasma ohiense SECH_94-Nam2_CI_9",
                                          "Histoplasma ohiense SECH_95.Nam2_CI_10","Histoplasma ohiense SECH_96.Nam2_CI_17",
                                          "Histoplasma ohiense SECH_97-Nam2_CI_18","Histoplasma ohiense SECH_98.Nam2_CI_30",
                                          "Histoplasma ohiense SECH_99.Nam2_CI_35","Histoplasma capsulatum HCH143",
                                          "Histoplasma capsulatum NACVFR_Histo_HC1070058","Histoplasma capsulatum HISSP-FGTRO0285",
                                          "Histoplasma capsulatum HISSP-FGPSO2043","Histoplasma capsulatum HISSP-FGPIE2055",
                                          "Histoplasma capsulatum HISSP-FGPIA2052","Histoplasma capsulatum HISSP-FGPERS2034",
                                          "Histoplasma capsulatum HISSP-FGMAR2044","Histoplasma capsulatum HISSP-FGLIN2055",
                                          "Histoplasma capsulatum HISSP-FGJOS2044","Histoplasma capsulatum HISSP-FGFIN2028",
                                          "Histoplasma capsulatum HISSP-FGFER2036","Histoplasma capsulatum HISSP-FGFAR0189",
                                          "Histoplasma capsulatum HISSP-FGFAN2059","Histoplasma capsulatum HISSP-FGBON2001",
                                          "Histoplasma capsulatum HISSP-FGBIK2051","Histoplasma capsulatum HISSP-FGAMA2041",
                                          "Histoplasma capsulatum HISSP-CM6408","Histoplasma capsulatum HISSP-CM6015",
                                          "Histoplasma capsulatum HISSP-B05821","Histoplasma capsulatum HISSP-11571-Belem1",
                                          "Histoplasma capsulatum HISSP-1014-Belem3","Histoplasma capsulatum 07_12-RJ",
                                          "Histoplasma capsulatum 104_p_06","Histoplasma capsulatum 117_p_12","Histoplasma capsulatum 122_p_10_B",
                                          "Histoplasma capsulatum 136_P_07","Histoplasma capsulatum 144_p_08","Histoplasma capsulatum 1517_p_17",
                                          "Histoplasma capsulatum 256_P_18","Histoplasma capsulatum 316_p_10","Histoplasma capsulatum 327_P_12",
                                          "Histoplasma capsulatum 343_p_18","Histoplasma capsulatum 388_p_11","Histoplasma capsulatum ES2_83Z",
                                          "Histoplasma capsulatum ES2_85Z","Histoplasma capsulatum ES2_86Z","Histoplasma capsulatum ES2_88Z",
                                          "Histoplasma capsulatum ES2_89Z","Histoplasma capsulatum ES2_90Z","Histoplasma capsulatum SA15",
                                          "Histoplasma capsulatum CB053","Histoplasma capsulatum CB055","Histoplasma capsulatum CB062",
                                          "Histoplasma capsulatum CB063","Histoplasma capsulatum CB064","Histoplasma capsulatum CB065",
                                          "Histoplasma capsulatum CB066","Histoplasma capsulatum CB168","Histoplasma capsulatum CB174",
                                          "Histoplasma capsulatum CB180","Histoplasma capsulatum CB186","Histoplasma capsulatum CB192",
                                          "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11","Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14",
                                          "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16"))

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Star_TE_GL.pdf', width=15, height=20)
ggplot(Summ_GL_RE_melt, (aes(x=variable, y=value, fill=values)))+ geom_bar(position="stack", stat="identity") + 
  coord_flip() + labs(title = "Histoplasma Repetative Elements & Starships vs Genome Length", y = "Genome Assembly (Mb)", x = "Histoplasma Isolates") + 
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired")
dev.off()

Stack_GL_TE

library("ggplot2")
library("ggtree")

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Star_TE_GL.pdf', width=15, height=20)
Stack_GL_TE + geom_facet(panel = "TE Abundance", data = Summ_GL_RE_melt, geom = geom_col, aes(x = value, color = value, fill = value), orientation = 'y', width = .6) 
dev.off()


new_seq3 <- tree_AAB$tip.label
d = data.frame(new_seq3)
d$Species = 0
d$Species[grep("UAMH130", d$new_seq3)] = "Blastomyces parvus UAMH130"
d$Species[grep("19VMG-15", d$new_seq3)] = "Histoplasma capsulatum 19VMG15"
d$Species[grep("07_12-RJ", d$new_seq3)] = "Histoplasma capsulatum 07_12-RJ"
d$Species[grep("HISSP-11571-Belem1", d$new_seq3)] = "Histoplasma capsulatum HISSP-11571-Belem1"
d$Species[grep("HISSP-FGFER2036", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGFER2036"
d$Species[grep("HISSP-CM6015", d$new_seq3)] = "Histoplasma capsulatum HISSP-CM6015"
d$Species[grep("HISSP-FGAMA2041", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGAMA2041"
d$Species[grep("HISSP-FGLIN2055", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGLIN2055"
d$Species[grep("HISSP-FGFAN2059", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGFAN2059"
d$Species[grep("HISSP-FGPIA2052", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGPIA2052"
d$Species[grep("HISSP-FGPIE2055", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGPIE2055"
d$Species[grep("HISSP-FGBON2001", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGBON2001"
d$Species[grep("HISSP-FGBIK2051", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGBIK2051"
d$Species[grep("HISSP-FGJOS2044", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGJOS2044"
d$Species[grep("JB_062632-Hc_062632", d$new_seq3)] = "Histoplasma capsulatum JB_062632-Hc_062632" 
d$Species[grep("NACVFR_Histo_HC1070058_2", d$new_seq3)] = "Histoplasma capsulatum NACVFR_Histo_HC1070058"
d$Species[grep("SECH_103", d$new_seq3)] = "Histoplasma suramericanum SECH_103_3_11G"
d$Species[grep("SECH_105", d$new_seq3)] = "Histoplasma suramericanum SECH_105_21_14"
d$Species[grep("HISSP-FGPERS2034", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGPERS2034"
d$Species[grep("HISSP-B05821", d$new_seq3)] = "Histoplasma capsulatum HISSP-B05821"
d$Species[grep("JB_031837-Hc_031837", d$new_seq3)] = "Histoplasma capsulatum JB_031837-Hc_031837"
d$Species[grep("JB_042430-Hc_042430", d$new_seq3)] = "Histoplasma capsulatum JB_042430-Hc_042430"
d$Species[grep("JB_062775-Hc_062775", d$new_seq3)] = "Histoplasma capsulatum JB_062775-Hc_062775"
d$Species[grep("JB_083285_2-Hc_083285_2", d$new_seq3)] = "Histoplasma capsulatum JB_083285_2-Hc_083285_2"
d$Species[grep("CB053-2", d$new_seq3)] = "Histoplasma capsulatum CB053" 
d$Species[grep("CB055-2", d$new_seq3)] = "Histoplasma capsulatum CB055"
d$Species[grep("SECH_101", d$new_seq3)] = "Histoplasma capsulatum senso_stricto SECH_101_H184AR"
d$Species[grep("HISSP-CM6408", d$new_seq3)] = "Histoplasma capsulatum HISSP-CM6408"
d$Species[grep("HISSP-1014-Belem3", d$new_seq3)] = "Histoplasma capsulatum HISSP-1014-Belem3"
d$Species[grep("ES2_83Z", d$new_seq3)] = "Histoplasma capsulatum ES2_83Z"
d$Species[grep("ES2_90Z", d$new_seq3)] = "Histoplasma capsulatum ES2_90Z"
d$Species[grep("ES2_86Z", d$new_seq3)] = "Histoplasma capsulatum ES2_86Z"
d$Species[grep("ES2_88Z", d$new_seq3)] = "Histoplasma capsulatum ES2_88Z"
d$Species[grep("HCH143", d$new_seq3)] = "Histoplasma capsulatum HCH143"
d$Species[grep("SECH_107", d$new_seq3)] = "Histoplasma capsulatum Africa SECH_107-duboisii-B-H88" 
d$Species[grep("ES2_85Z", d$new_seq3)] = "Histoplasma capsulatum ES2_85Z"
d$Species[grep("HISSP-CM7256", d$new_seq3)] = "Histoplasma capsulatum HISSP-CM7256"
d$Species[grep("ES2_89Z", d$new_seq3)] = "Histoplasma capsulatum ES2_89Z"
d$Species[grep("SA15", d$new_seq3)] = "Histoplasma capsulatum SA15"
d$Species[grep("SECH_98", d$new_seq3)] = "Histoplasma ohiense SECH_98.Nam2_CI_30"
d$Species[grep("SECH_110", d$new_seq3)] = "Histoplasma capsulatum SECH_110"
d$Species[grep("G217B", d$new_seq3)] = "Histoplasma ohiense SECH_91-Nam2_G217B"
d$Species[grep("SECH_92", d$new_seq3)] = "Histoplasma ohiense SECH_92-Nam2_G222B"
d$Species[grep("Hc1986", d$new_seq3)] = "Histoplasma ohiense Hc1986"
d$Species[grep("SECH_82", d$new_seq3)] = "Histoplasma ohiense SECH_82-Nam1_CI_4"
d$Species[grep("SECH_95", d$new_seq3)] = "Histoplasma ohiense SECH_95.Nam2_CI_10"
d$Species[grep("SECH_94", d$new_seq3)] = "Histoplasma ohiense SECH_94-Nam2_CI_9"
d$Species[grep("SECH_96", d$new_seq3)] = "Histoplasma ohiense SECH_96.Nam2_CI_17"
d$Species[grep("SECH_93", d$new_seq3)] = "Histoplasma ohiense SECH_93.Nam2_CI_6"
d$Species[grep("SECH_97", d$new_seq3)] = "Histoplasma ohiense SECH_97-Nam2_CI_18"
d$Species[grep("SECH_99", d$new_seq3)] = "Histoplasma ohiense SECH_99.Nam2_CI_35"
d$Species[grep("SECH_104", d$new_seq3)] = "Histoplasma suramericanum SECH_104_27_14" 
d$Species[grep("CB174-2", d$new_seq3)] = "Histoplasma capsulatum CB174"
d$Species[grep("SECH_109", d$new_seq3)] = "Histoplasma capsulatum SECH_109"
d$Species[grep("CB062-2", d$new_seq3)] = "Histoplasma capsulatum CB062"
d$Species[grep("CB065-2", d$new_seq3)] = "Histoplasma capsulatum CB065"
d$Species[grep("CB192-2", d$new_seq3)] = "Histoplasma capsulatum CB192"
d$Species[grep("CB064-2", d$new_seq3)] = "Histoplasma capsulatum CB064"
d$Species[grep("CB180-2", d$new_seq3)] = "Histoplasma capsulatum CB180"
d$Species[grep("CB063-2", d$new_seq3)] = "Histoplasma capsulatum CB063"
d$Species[grep("CB186-2", d$new_seq3)] = "Histoplasma capsulatum CB186" 
d$Species[grep("SECH_81", d$new_seq3)] = "Histoplasma mississippiense SECH_81-Nam1_WU24"
d$Species[grep("CB066-2", d$new_seq3)] = "Histoplasma capsulatum CB066"
d$Species[grep("SECH_88", d$new_seq3)] = "Histoplasma mississippiense SECH_88-Nam1_CI_43"
d$Species[grep("SECH_86", d$new_seq3)] = "Histoplasma mississippiense SECH_86-Nam1_CI_24"
d$Species[grep("CB168-2", d$new_seq3)] = "Histoplasma capsulatum CB168"
d$Species[grep("SECH_85", d$new_seq3)] = "Histoplasma mississippiense SECH_85-Nam1_CI_22"
d$Species[grep("SECH_83", d$new_seq3)] = "Histoplasma mississippiense SECH_83-Nam1_CI_7" 
d$Species[grep("SECH_87", d$new_seq3)] = "Histoplasma mississippiense SECH_87-Nam1_CI_42"
d$Species[grep("SECH_90", d$new_seq3)] = "Histoplasma mississippiense SECH_90-Nam1_DOWNS"
d$Species[grep("SECH_102", d$new_seq3)] = "Histoplasma mississippiense SECH_102_Nam2_505"
d$Species[grep("SECH_100", d$new_seq3)] = "Histoplasma capsulatum senso_stricto SECH_100_H186AR"
d$Species[grep("SECH_84", d$new_seq3)] = "Histoplasma mississippiense SECH_84-Nam1_CI_19"
d$Species[grep("SECH_89", d$new_seq3)] = "Histoplasma mississippiense SECH_89-Nam1_UCLA-531"
d$Species[grep("JB_01752-Hc_01752", d$new_seq3)] = "Histoplasma capsulatum JB_01752-Hc_01752"
d$Species[grep("104_p_06", d$new_seq3)] = "Histoplasma capsulatum 104_p_06"
d$Species[grep("1517_p_17", d$new_seq3)] = "Histoplasma capsulatum 1517_p_17"
d$Species[grep("117_p_12", d$new_seq3)] = "Histoplasma capsulatum 117_p_12"
d$Species[grep("388_p_11", d$new_seq3)] = "Histoplasma capsulatum 388_p_11" 
d$Species[grep("144_p_08", d$new_seq3)] = "Histoplasma capsulatum 144_p_08"
d$Species[grep("327_P_12", d$new_seq3)] = "Histoplasma capsulatum 327_P_12"
d$Species[grep("316_p_10", d$new_seq3)] = "Histoplasma capsulatum 316_p_10"
d$Species[grep("Dr_Anuradha_Fungal_WGS_S14", d$new_seq3)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14"
d$Species[grep("Histo-485P20", d$new_seq3)] = "Histoplasma capsulatum Histo-485P20"
d$Species[grep("136_P_07", d$new_seq3)] = "Histoplasma capsulatum 136_P_07"
d$Species[grep("343_p_18", d$new_seq3)] = "Histoplasma capsulatum 343_p_18"
d$Species[grep("Dr_Anuradha_Fungal_WGS_S11", d$new_seq3)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S11"
d$Species[grep("Dr_Anuradha_Fungal_WGS_S16", d$new_seq3)] = "Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16"
d$Species[grep("122_p_10_B", d$new_seq3)] = "Histoplasma capsulatum 122_p_10_B"
d$Species[grep("256_P_18", d$new_seq3)] = "Histoplasma capsulatum 256_P_18" 
d$Species[grep("JB_021091-Hc_021091", d$new_seq3)] = "Histoplasma capsulatum JB_021091-Hc_021091"
d$Species[grep("JB_073129-Hc_073129", d$new_seq3)] = "Histoplasma capsulatum JB_073129-Hc_073129"
d$Species[grep("HISSP-FGFAR0189", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGFAR0189"
d$Species[grep("HISSP-FGFIN2028", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGFIN2028"
d$Species[grep("HISSP-FGPSO2043", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGPSO2043"
d$Species[grep("HISSP-FGMAR2044", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGMAR2044"
d$Species[grep("HISSP-FGTRO0285", d$new_seq3)] = "Histoplasma capsulatum HISSP-FGTRO0285"
tree_AAB$tip.label <- d$Species

ggtree(tree_AAB, layout = "rectangular")
ggtree(tree_AAB, layout = "fan")

# Visualize the tree
p_AABBB1 <- ggtree(tree_AAB, layout = "rectangular") + geom_tiplab(size = 3, fontface = 3, align=TRUE) +
  theme_tree2() + theme(legend.position = "none") 
print(p_AABBB1)
# Save the tree as a PDF
ggsave("AAB_ultrametric_tree.pdf", p_AAB, width = 20, height = 5, dpi = 150)
ggsave(filename = "species_tree_with_numbersB.png", plot = p_AAB,width = 9, height = 6, units = c("in"), dpi = 300, scale =2, device = "png")

#Sanity check!
colnames(Percentage_B)[which(colnames(Percentage_B)  %in% tree_AAB$tip.label==F)]

colnames(PERC_SUM_melt_B)
PERC_SUM_melt_B2 <- PERC_SUM_melt_B[,c(2,1,3)]

#pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_fan_Stacked_barplot_AAB-Star.pdf', width=15, height=15)

p_AABBB2 <- p_AABBB1 + geom_facet(panel = "TE Abundance", data = PERC_SUM_melt_B2, geom = geom_col, 
                                aes(x = value, color = class, 
                                    fill = class), orientation = 'y', width = .6) + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.05) +
  theme_tree2(legend.position=c(.95, .25)) + xlim_expand(c(0,0.8), 'TE Abundance') + xlim_expand(c(0, 4), 'Tree') + geom_cladelabel(node=140, label="NAm1", color="red", offset = 1.5) +
  geom_cladelabel(node=160, label="NAm2", color="#808000", offset = 1.5) + geom_cladelabel(node=172, label="India", color="blue", offset = 1.5) + geom_cladelabel(node=113, label="Africa", color="#FFA500", offset = 1.5) +
  geom_cladelabel(node=105, label="LAmA", color="purple", offset = 1.5) + geom_cladelabel(node=121, label="LAmA", color="purple", offset = 1.5) + geom_cladelabel(node=187, label="LAmA", color="purple", offset = 1.5) +
  geom_cladelabel(node=191, label="Africa", color="#FFA500", offset = 1.5) + geom_cladelabel(node=137, label="LAmA", color="purple", offset = 1.5) + geom_cladelabel(node=1, label="Root", color="green", offset = 1.5)

#dev.off()

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_Stacked_barplotAAB-LINES_STAR.pdf', width=15, height=15)

facet_widths(p_AABBB2, widths = c(1, 0.40))

dev.off()

png(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_fan_Stacked_barplotAAB-LINES_STAR.png', width=1700, height=1500, res=100)

facet_widths(p_AABBB2, widths = c(1, 0.40))

dev.off()

#if you need to figure out what the nodes are in your tree
NODES_TREE<-ggtree(tree_AAB) + geom_text(aes(label=node), hjust=-.3)

## Starship Data Visualization

StarHisto <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/starship_count2.tsv")

#melt data
Star_melt<- melt(StarHisto, id.vars="class", 
                      measure.vars=c("Histoplasma_capsulatum_19VMG-15",
                                     "Histoplasma_capsulatum_HISSP-CM7256",
                                     "Histoplasma_capsulatum_Histo-485P20",
                                     "Histoplasma_capsulatum_JB_01752-Hc_01752",
                                     "Histoplasma_capsulatum_JB_021091-Hc_021091",
                                     "Histoplasma_capsulatum_JB_031837-Hc_031837",
                                     "Histoplasma_capsulatum_JB_042430-Hc_042430",
                                     "Histoplasma_capsulatum_JB_062632-Hc_062632",
                                     "Histoplasma_capsulatum_JB_062775-Hc_062775",
                                     "Histoplasma_capsulatum_JB_073129-Hc_073129",
                                     "Histoplasma_capsulatum_JB_083285_2-Hc_083285_2",
                                     "Histoplasma_capsulatum_SECH_100-Nam2_G186A",
                                     "Histoplasma_capsulatum_SECH_101-Nam2_G184A",
                                     "Histoplasma_mississippiense_SECH_102Nam2_505",
                                     "Histoplasma_suramericanum_SECH_103-Nam2_3_11G",
                                     "Histoplasma_suramericanum_SECH_104-Nam2_27_14",
                                     "Histoplasma_suramericanum_SECH_105-Nam2_21_14",
                                     "Histoplasma_capsulatum_SECH_107-mis_Hc_duboisii-B",
                                     "Histoplasma_capsulatum_SECH_109",
                                     "Histoplasma_capsulatum_SECH_110",
                                     "Histoplasma_mississippiense_SECH_81-Nam1_WU24",
                                     "Histoplasma_ohiense_SECH_82-Nam1_CI_4",
                                     "Histoplasma_mississippiense_SECH_83-Nam1_CI_7",
                                     "Histoplasma_mississippiense_SECH_84-Nam1_CI_19",
                                     "Histoplasma_mississippiense_SECH_85-Nam1_CI_22",
                                     "Histoplasma_mississippiense_SECH_86-Nam1_CI_24",
                                     "Histoplasma_mississippiense_SECH_87-Nam1_CI_42",
                                     "Histoplasma_mississippiense_SECH_88-Nam1_CI_43",
                                     "Histoplasma_mississippiense_SECH_89-Nam1_UCLA-531",
                                     "Histoplasma_mississippiense_SECH_90-Nam1_DOWNS",
                                     "Histoplasma_ohiense_SECH_91-Nam2_G217B",
                                     "Histoplasma_ohiense_SECH_92-Nam2_G222B",
                                     "Histoplasma_ohiense_SECH_93.Nam2_CI_6",
                                     "Histoplasma_ohiense_SECH_94-Nam2_CI_9",
                                     "Histoplasma_ohiense_SECH_95.Nam2_CI_10",
                                     "Histoplasma_ohiense_SECH_96.Nam2_CI_17",
                                     "Histoplasma_ohiense_SECH_97-Nam2_CI_18",
                                     "Histoplasma_ohiense_SECH_98.Nam2_CI_30",
                                     "Histoplasma_ohiense_SECH_99-Nam2_CI_35",
                                     "Histoplasma_capsulatum_HCH143",
                                     "Histoplasma_capsulatum_NACVFR_Histo_HC1070058",
                                     "Histoplasma_capsulatum_HISSP-FGTRO0285",
                                     "Histoplasma_capsulatum_HISSP-FGPSO2043",
                                     "Histoplasma_capsulatum_HISSP-FGPIE2055",
                                     "Histoplasma_capsulatum_HISSP-FGPIA2052",
                                     "Histoplasma_capsulatum_HISSP-FGPERS2034",
                                     "Histoplasma_capsulatum_HISSP-FGMAR2044",
                                     "Histoplasma_capsulatum_HISSP-FGLIN2055",
                                     "Histoplasma_capsulatum_HISSP-FGJOS2044",
                                     "Histoplasma_capsulatum_HISSP-FGFIN2028",
                                     "Histoplasma_capsulatum_HISSP-FGFER2036",
                                     "Histoplasma_capsulatum_HISSP-FGFAR0189",
                                     "Histoplasma_capsulatum_HISSP-FGFAN2059",
                                     "Histoplasma_capsulatum_HISSP-FGBON2001",
                                     "Histoplasma_capsulatum_HISSP-FGBIK2051",
                                     "Histoplasma_capsulatum_HISSP-FGAMA2041",
                                     "Histoplasma_capsulatum_HISSP-CM6408",
                                     "Histoplasma_capsulatum_HISSP-CM6015",
                                     "Histoplasma_capsulatum_HISSP-B05821",
                                     "Histoplasma_capsulatum_HISSP-11571-Belem1",
                                     "Histoplasma_capsulatum_HISSP-1014-Belem3",
                                     "Histoplasma_capsulatum_07_12-RJ",
                                     "Histoplasma_capsulatum_104_p_06",
                                     "Histoplasma_capsulatum_117_p_12",
                                     "Histoplasma_capsulatum_122_p_10_B",
                                     "Histoplasma_capsulatum_136_P_07",
                                     "Histoplasma_capsulatum_144_p_08",
                                     "Histoplasma_capsulatum_1517_p_17",
                                     "Histoplasma_capsulatum_256_P_18",
                                     "Histoplasma_capsulatum_316_p_10",
                                     "Histoplasma_capsulatum_327_P_12",
                                     "Histoplasma_capsulatum_343_p_18",
                                     "Histoplasma_capsulatum_388_p_11",
                                     "Histoplasma_capsulatum_ES2_83Z",
                                     "Histoplasma_capsulatum_ES2_85Z",
                                     "Histoplasma_capsulatum_ES2_86Z",
                                     "Histoplasma_capsulatum_ES2_88Z",
                                     "Histoplasma_capsulatum_ES2_89Z",
                                     "Histoplasma_capsulatum_ES2_90Z",
                                     "Histoplasma_capsulatum_SA15",
                                     "Histoplasma_capsulatum_CB053",
                                     "Histoplasma_capsulatum_CB055",
                                     "Histoplasma_capsulatum_CB062",
                                     "Histoplasma_capsulatum_CB063",
                                     "Histoplasma_capsulatum_CB064",
                                     "Histoplasma_capsulatum_CB065",
                                     "Histoplasma_capsulatum_CB066",
                                     "Histoplasma_capsulatum_CB168",
                                     "Histoplasma_capsulatum_CB174",
                                     "Histoplasma_capsulatum_CB180",
                                     "Histoplasma_capsulatum_CB186",
                                     "Histoplasma_capsulatum_CB192",
                                     "Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11",
                                     "Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14",
                                     "Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16",
                                     "Histoplasma_mississippiense_WU24",
                                     "Histoplasma_ohiense_G217B",
                                     "Histoplasma_capsulatum_G184AR",
                                     "Histoplasma_capsulatum_G186AR",
                                     "Histoplasma_capsulatum_H88"))



library("ggplot2")


#Saving results to a PDF file
pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/Starship_genome.pdf', width=8, height=8)

ggplot(Star_melt, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of Starship per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired")

dev.off()


#RE-RUN significance tests on tree_AAB

colnames(Percentage_B) = gsub("[.]", "_", colnames(Percentage_B))
tree_AAB$tip.label = gsub("[-]", "_", tree_AAB$tip.label)
all(colnames(Percentage_B)  %in% tree_AAB$tip.label)
colnames(Percentage_B)[which(colnames(Percentage_B)  %in% tree_AAB$tip.label==F)]
#Just first column doesn't agree which is class

#Test TE for phylogenetic signal - but phylosig requires named vector for data
row_DNA = which(Percentage_B$class=="DNA TE")
data.vec_DNA = as.numeric(Percentage_B[row_DNA, ]); names(data.vec_DNA) = colnames(Percentage_B)

# remove the first column, which contains no data
data.vec_DNA=data.vec_DNA[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.DNA = phylosig(tree=tree_AAB, x=data.vec_DNA, method="lambda", test=TRUE)
sig.lam.DNA
plot.phylosig(sig.lam.DNA)


#Phylogenetic signal lambda : 7.33137e-05 
#logL(lambda) : 235.542
#LR(lambda=0) : -0.000662666  
#P-value (based on LR test) : 1

sig.k.DNA = (phylosig(tree=tree_AAB, x=data.vec_DNA, method="K", test=TRUE))
sig.k.DNA
plot.phylosig(sig.k.DNA)


#Phylogenetic signal K : 0.359745
#P-value (based on 1000 randomizations) : 0.464

### RNA-TE
row_RNA = which(Percentage_B$class=="RNA TE")
data.vec_RNA = as.numeric(Percentage_B[row_RNA, ]); names(data.vec_RNA) = colnames(Percentage_B)

# remove the first column, which contains no data
data.vec_RNA=data.vec_RNA[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.RNA = phylosig(tree=tree_AAB, x=data.vec_RNA, method="lambda", test=TRUE)
sig.lam.RNA
plot.phylosig(sig.lam.RNA)

#Phylogenetic signal lambda : 0.805879
#logL(lambda) : 114.562
#LR(lambda=0) : 41.2796 
#P-value (based on LR test) : 1.31942e-10 

sig.k.RNA = (phylosig(tree=tree_AAB, x=data.vec_RNA, method="K", test=TRUE))
sig.k.RNA
plot.phylosig(sig.k.RNA)
#Phylogenetic signal K : 0.717039 
#P-value (based on 1000 randomizations) : 0.001

### LC-TE
row_LC = which(Percentage_B$class=="Low_complexity")
data.vec_LC = as.numeric(Percentage_B[row_LC, ]); names(data.vec_LC) = colnames(Percentage_B)

# remove the first column, which contains no data
data.vec_LC=data.vec_LC[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.LC = phylosig(tree=tree_AAB, x=data.vec_LC, method="lambda", test=TRUE)
sig.lam.LC
plot.phylosig(sig.lam.LC)


#Phylogenetic signal lambda : 0.80018
#logL(lambda) : 547.847 
#LR(lambda=0) : 69.2144 
#P-value (based on LR test) : 1.31942e-10 

sig.k.LC = (phylosig(tree=tree_AAB, x=data.vec_LC, method="K", test=TRUE))
sig.k.LC
plot.phylosig(sig.k.LC)

#Phylogenetic signal K : 0.994923 
#P-value (based on 1000 randomizations) : 0.001 

### MITE
row_MITE = which(Percentage_B$class=="MITE")
data.vec_MITE = as.numeric(Percentage_B[row_MITE, ]); names(data.vec_MITE) = colnames(Percentage_B)

# remove the first column, which contains no data
data.vec_MITE=data.vec_MITE[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.MITE = phylosig(tree=tree_AAB, x=data.vec_MITE, method="lambda", test=TRUE)
sig.lam.MITE
plot.phylosig(sig.lam.MITE)

#Phylogenetic signal lambda : 7.33137e-05 
#logL(lambda) : 447.566 
#LR(lambda=0) : -0.00197191 
#P-value (based on LR test) : 1

sig.k.MITE = (phylosig(tree=tree_AAB, x=data.vec_MITE, method="K", test=TRUE))
sig.k.MITE
plot.phylosig(sig.k.MITE)

#Phylogenetic signal K : 0.351299 
#P-value (based on 1000 randomizations) : 0.513 

### Satellite
row_SAT = which(Percentage_B$class=="Satellite")
data.vec_SAT = as.numeric(Percentage_B[row_SAT, ]); names(data.vec_SAT) = colnames(Percentage_B)

# remove the first column, which contains no data
data.vec_SAT=data.vec_SAT[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.SAT = phylosig(tree=tree_AAB, x=data.vec_SAT, method="lambda", test=TRUE)
sig.lam.SAT
plot.phylosig(sig.lam.SAT)

#Phylogenetic signal lambda : 7.33137e-05 
#logL(lambda) : 1137.67 
#LR(lambda=0) : -0.00217452 
#P-value (based on LR test) : 1 

sig.k.SAT = (phylosig(tree=tree_AAB, x=data.vec_SAT, method="K", test=TRUE))
sig.k.SAT
plot.phylosig(sig.k.SAT)

#Phylogenetic signal K : 0.379067 
#P-value (based on 1000 randomizations) : 0.441 

### Simple Repeat
row_SR = which(Percentage_B$class=="Simple_repeat")
data.vec_SR = as.numeric(Percentage_B[row_SR, ]); names(data.vec_SR) = colnames(Percentage_B)

# remove the first column, which contains no data
data.vec_SR=data.vec_SR[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.SR = phylosig(tree=tree_AAB, x=data.vec_SR, method="lambda", test=TRUE)
sig.lam.SR
plot.phylosig(sig.lam.SR)

#Phylogenetic signal lambda : 0.861457 
#logL(lambda) : 422.439 
#LR(lambda=0) : 82.4716 
#P-value (based on LR test) : 1.072e-19

sig.k.SR = (phylosig(tree=tree_AAB, x=data.vec_SR, method="K", test=TRUE))
sig.k.SR
plot.phylosig(sig.k.SR)

#Phylogenetic signal K : 1.18526 
#P-value (based on 1000 randomizations) : 0.001 

## Starships
row_Star = which(Percentage_B$class=="Starship")
data.vec_Star = as.numeric(Percentage_B[row_Star, ]); names(data.vec_Star) = colnames(Percentage_B)

# remove the first column, which contains no data
data.vec_Star=data.vec_Star[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.Star = phylosig(tree=tree_AAB, x=data.vec_Star, method="lambda", test=TRUE)
sig.lam.Star
plot.phylosig(sig.lam.Star)

#Phylogenetic signal lambda : 0.911629 
#logL(lambda) : 321.212 
#LR(lambda=0) : 112.234 
#P-value (based on LR test) : 3.17469e-26  

sig.k.Star = (phylosig(tree=tree_AAB, x=data.vec_Star, method="K", test=TRUE))
sig.k.Star
plot.phylosig(sig.k.Star)

#Phylogenetic signal K : 1.49134 
#P-value (based on 1000 randomizations) : 0.001 
