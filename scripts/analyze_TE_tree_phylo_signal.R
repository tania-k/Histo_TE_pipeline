#Large repetative and messy script which will help you collect and assess your TE results to a stacked barplot visual.
#First calculate TEs, collect all varieties per isolate, combine into large table
#I worked the final table on my computer, had the script read that in to do final TE analysis.
#Trees are built here too, but had two CDS and AA trees, went with the AA tree.
#Created trees using Aspergillus fumigatus, Emmonsia and Blastomyces as the root species. Final tree uses Blastomyces as root.
#Also added Phylogenetic signal to each rooted species. Make sure to run the phylogenetic signal analysis for each TE type.
#I had 7 groups done, RNA TEs, DNA TEs, Simple Repeats, Low Complexity, MITEs, Satellites, and Unknowns. I've seen many disregard Low complexity, simple repeats and unknowns. Up to your discretion. 


library(tidyverse)
library(readr)
library(dplyr)
library(plyr)
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
write_tsv(Histoplasma_capsulatum_WU24class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_WU24class.tsv")
Histoplasma_capsulatum_WU24class <- Histoplasma_capsulatum_WU24class %>% mutate('Histoplasma capsulatum WU24'= percentage)

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
write_tsv(Histoplasma_ohiense_G217Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_G217Bclass.tsv")
Histoplasma_ohiense_G217Bclass <- Histoplasma_ohiense_G217Bclass %>% mutate('Histoplasma ohiense G217B'= percentage)
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
write_tsv(Histoplasma_capsulatum_var.duboisii_H88class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_var.duboisii_H88class.tsv")
Histoplasma_capsulatum_var.duboisii_H88class <- Histoplasma_capsulatum_var.duboisii_H88class %>% mutate('Histoplasma capsulatum var.duboisii H88'= percentage)
Histoplasma_capsulatum_var.duboisii_H88class1 <- Histoplasma_capsulatum_var.duboisii_H88class %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_G186AR <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_G186AR.RM/Histoplasma_capsulatum_G186AR.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_G186AR) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_G186AR  <- Histoplasma_capsulatum_G186AR  %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_G186AR,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_G186AR.tsv")

Histoplasma_capsulatum_G186ARclass <- Histoplasma_capsulatum_G186AR %>% group_by(class) %>% 
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/31111494)
write_tsv(Histoplasma_capsulatum_G186ARclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_G186ARclass.tsv")
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
write_tsv(Histoplasma_capsulatum_G184ARclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_G184ARclass.tsv")
Histoplasma_capsulatum_G184ARclass <- Histoplasma_capsulatum_G184ARclass %>% mutate('Histoplasma capsulatum G184AR'= percentage)
Histoplasma_capsulatum_G184ARclass1 <- Histoplasma_capsulatum_G184ARclass %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_HCCI_17 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_HCCI_17.RM/Histoplasma_ohiense_HCCI_17.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_HCCI_17) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_HCCI_17 <- Histoplasma_ohiense_HCCI_17 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_HCCI_17,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCCI_17_version2.tsv")

Histoplasma_ohiense_HCCI_17class<- Histoplasma_ohiense_HCCI_17 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/34006939)
write_tsv(Histoplasma_ohiense_HCCI_17class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCCI_17class_version2.tsv")
Histoplasma_ohiense_HCCI_17class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCCI_17class_version2.tsv")
Histoplasma_ohiense_HCCI_17class <- Histoplasma_ohiense_HCCI_17class %>% mutate('Histoplasma ohiense HCCI_17'= percentage)
Histoplasma_ohiense_HCCI_17class1 <- Histoplasma_ohiense_HCCI_17class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_HCCI_6 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_HCCI_6.RM/Histoplasma_ohiense_HCCI_6.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_HCCI_6) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_HCCI_6 <- Histoplasma_ohiense_HCCI_6 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_HCCI_6,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCCI_6.tsv")

Histoplasma_ohiense_HCCI_6class<- Histoplasma_ohiense_HCCI_6 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28445673)
write_tsv(Histoplasma_ohiense_HCCI_6class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCCI_6class.tsv")
Histoplasma_ohiense_HCCI_6class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCCI_6class.tsv")
Histoplasma_ohiense_HCCI_6class <- Histoplasma_ohiense_HCCI_6class %>% mutate('Histoplasma ohiense HCCI_6'= percentage)
Histoplasma_ohiense_HCCI_6class1 <- Histoplasma_ohiense_HCCI_6class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_HCG217B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_HCG217B.RM/Histoplasma_ohiense_HCG217B.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_HCG217B) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_HCG217B <- Histoplasma_ohiense_HCG217B %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_HCG217B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCG217B.tsv")

Histoplasma_ohiense_HCG217Bclass<- Histoplasma_ohiense_HCG217B %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/33872679)
write_tsv(Histoplasma_ohiense_HCG217Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCG217Bclass.tsv")
Histoplasma_ohiense_HCG217Bclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_HCG217Bclass.tsv")
Histoplasma_ohiense_HCG217Bclass <- Histoplasma_ohiense_HCG217Bclass %>% mutate('Histoplasma ohiense HCG217B'= percentage)
Histoplasma_ohiense_HCG217Bclass1 <- Histoplasma_ohiense_HCG217Bclass %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_82.Nam1_CI_4 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_82-Nam1_CI_4.RM/Histoplasma_ohiense_SECH_82-Nam1_CI_4.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_82.Nam1_CI_4) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_82.Nam1_CI_4 <- Histoplasma_ohiense_SECH_82.Nam1_CI_4 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_82.Nam1_CI_4,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_82-Nam1_CI_4.tsv")

Histoplasma_ohiense_SECH_82.Nam1_CI_4class<- Histoplasma_ohiense_SECH_82.Nam1_CI_4 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/30880349)
write_tsv(Histoplasma_ohiense_SECH_82.Nam1_CI_4class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_82-Nam1_CI_4class.tsv")
Histoplasma_ohiense_SECH_82.Nam1_CI_4class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_82-Nam1_CI_4class.tsv")
Histoplasma_ohiense_SECH_82.Nam1_CI_4class <- Histoplasma_ohiense_SECH_82.Nam1_CI_4class %>% mutate('Histoplasma ohiense SECH_82-Nam1_CI_4'= percentage)
Histoplasma_ohiense_SECH_82.Nam1_CI_4class1 <- Histoplasma_ohiense_SECH_82.Nam1_CI_4class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_91.Nam2_G217B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_91-Nam2_G217B.RM/Histoplasma_ohiense_SECH_91-Nam2_G217B.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_91.Nam2_G217B) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_91.Nam2_G217B <- Histoplasma_ohiense_SECH_91.Nam2_G217B %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_91.Nam2_G217B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_91-Nam2_G217B.tsv")

Histoplasma_ohiense_SECH_91.Nam2_G217Bclass<- Histoplasma_ohiense_SECH_91.Nam2_G217B %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29816969)
write_tsv(Histoplasma_ohiense_SECH_91.Nam2_G217Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_91-Nam2_G217Bclass.tsv")
Histoplasma_ohiense_SECH_91.Nam2_G217Bclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_91-Nam2_G217Bclass.tsv")
Histoplasma_ohiense_SECH_91.Nam2_G217Bclass <- Histoplasma_ohiense_SECH_91.Nam2_G217Bclass %>% mutate('Histoplasma ohiense SECH_91-Nam2_G217B'= percentage)
Histoplasma_ohiense_SECH_91.Nam2_G217Bclass1 <- Histoplasma_ohiense_SECH_91.Nam2_G217Bclass %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_92.Nam2_G222B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_92-Nam2_G222B.RM/Histoplasma_ohiense_SECH_92-Nam2_G222B.scaffolds.fa.out",skip=3,col_names = F)

colnames(Histoplasma_ohiense_SECH_92.Nam2_G222B) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_92.Nam2_G222B <- Histoplasma_ohiense_SECH_92.Nam2_G222B %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_92.Nam2_G222B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_92-Nam2_G222B.tsv")

Histoplasma_ohiense_SECH_92.Nam2_G222Bclass<- Histoplasma_ohiense_SECH_92.Nam2_G222B %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29639711)
write_tsv(Histoplasma_ohiense_SECH_92.Nam2_G222Bclass,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_92-Nam2_G222Bclass.tsv")
Histoplasma_ohiense_SECH_92.Nam2_G222Bclass <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_92-Nam2_G222Bclass.tsv")
Histoplasma_ohiense_SECH_92.Nam2_G222Bclass <- Histoplasma_ohiense_SECH_92.Nam2_G222Bclass %>% mutate('Histoplasma ohiense SECH_92-Nam2_G222B'= percentage)
Histoplasma_ohiense_SECH_92.Nam2_G222Bclass1 <- Histoplasma_ohiense_SECH_92.Nam2_G222Bclass %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_93.Nam2_CI_6 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_93-Nam2_CI_6.RM/Histoplasma_ohiense_SECH_93-Nam2_CI_6.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_93.Nam2_CI_6) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_93.Nam2_CI_6 <- Histoplasma_ohiense_SECH_93.Nam2_CI_6 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_93.Nam2_CI_6,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_93-Nam2_CI_6.tsv")

Histoplasma_ohiense_SECH_93.Nam2_CI_6class<- Histoplasma_ohiense_SECH_93.Nam2_CI_6 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23646057)
write_tsv(Histoplasma_ohiense_SECH_93.Nam2_CI_6class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_93-Nam2_CI_6class.tsv")
Histoplasma_ohiense_SECH_93.Nam2_CI_6class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_93-Nam2_CI_6class.tsv")
Histoplasma_ohiense_SECH_93.Nam2_CI_6class <- Histoplasma_ohiense_SECH_93.Nam2_CI_6class %>% mutate('Histoplasma ohiense SECH_93.Nam2_CI_6'= percentage)
Histoplasma_ohiense_SECH_93.Nam2_CI_6class1 <- Histoplasma_ohiense_SECH_93.Nam2_CI_6class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_94.Nam2_CI_9 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_94-Nam2_CI_9.RM/Histoplasma_ohiense_SECH_94-Nam2_CI_9.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_94.Nam2_CI_9) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_94.Nam2_CI_9 <- Histoplasma_ohiense_SECH_94.Nam2_CI_9 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_94.Nam2_CI_9,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_94.Nam2_CI_9.tsv")

Histoplasma_ohiense_SECH_94.Nam2_CI_9class<- Histoplasma_ohiense_SECH_94.Nam2_CI_9 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/20233006)
write_tsv(Histoplasma_ohiense_SECH_94.Nam2_CI_9class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_94.Nam2_CI_9class.tsv")
Histoplasma_ohiense_SECH_94.Nam2_CI_9class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_94.Nam2_CI_9class.tsv")
Histoplasma_ohiense_SECH_94.Nam2_CI_9class <- Histoplasma_ohiense_SECH_94.Nam2_CI_9class %>% mutate('Histoplasma ohiense SECH_94-Nam2_CI_9'= percentage)
Histoplasma_ohiense_SECH_94.Nam2_CI_9class1 <- Histoplasma_ohiense_SECH_94.Nam2_CI_9class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_95.Nam2_CI_10 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_95-Nam2_CI_10.RM/Histoplasma_ohiense_SECH_95-Nam2_CI_10.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_95.Nam2_CI_10) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_95.Nam2_CI_10 <- Histoplasma_ohiense_SECH_95.Nam2_CI_10 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_95.Nam2_CI_10,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_95.Nam2_CI_10.tsv")

Histoplasma_ohiense_SECH_95.Nam2_CI_10class<- Histoplasma_ohiense_SECH_95.Nam2_CI_10 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25510614)
write_tsv(Histoplasma_ohiense_SECH_95.Nam2_CI_10class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_95.Nam2_CI_10class.tsv")
Histoplasma_ohiense_SECH_95.Nam2_CI_10class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_95.Nam2_CI_10class.tsv")
Histoplasma_ohiense_SECH_95.Nam2_CI_10class <- Histoplasma_ohiense_SECH_95.Nam2_CI_10class %>% mutate('Histoplasma ohiense SECH_95.Nam2_CI_10'= percentage)
Histoplasma_ohiense_SECH_95.Nam2_CI_10class1 <- Histoplasma_ohiense_SECH_95.Nam2_CI_10class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_96.Nam2_CI_17 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_96-Nam2_CI_17.RM/Histoplasma_ohiense_SECH_96-Nam2_CI_17.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_96.Nam2_CI_17) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_96.Nam2_CI_17 <- Histoplasma_ohiense_SECH_96.Nam2_CI_17 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_96.Nam2_CI_17,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_96.Nam2_CI_17.tsv")

Histoplasma_ohiense_SECH_96.Nam2_CI_17class<- Histoplasma_ohiense_SECH_96.Nam2_CI_17 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25893578)
write_tsv(Histoplasma_ohiense_SECH_96.Nam2_CI_17class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_96.Nam2_CI_17class.tsv")
Histoplasma_ohiense_SECH_96.Nam2_CI_17class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_96.Nam2_CI_17class.tsv")
Histoplasma_ohiense_SECH_96.Nam2_CI_17class <- Histoplasma_ohiense_SECH_96.Nam2_CI_17class %>% mutate('Histoplasma ohiense SECH_96.Nam2_CI_17'= percentage)
Histoplasma_ohiense_SECH_96.Nam2_CI_17class1 <- Histoplasma_ohiense_SECH_96.Nam2_CI_17class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_97.Nam2_CI_18 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_97-Nam2_CI_18.RM/Histoplasma_ohiense_SECH_97-Nam2_CI_18.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_97.Nam2_CI_18) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_97.Nam2_CI_18 <- Histoplasma_ohiense_SECH_97.Nam2_CI_18 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_97.Nam2_CI_18,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_97-Nam2_CI_18.tsv")

Histoplasma_ohiense_SECH_97.Nam2_CI_18class<- Histoplasma_ohiense_SECH_97.Nam2_CI_18 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28259788)
write_tsv(Histoplasma_ohiense_SECH_97.Nam2_CI_18class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_97-Nam2_CI_18class.tsv")
Histoplasma_ohiense_SECH_97.Nam2_CI_18class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_97-Nam2_CI_18class.tsv")
Histoplasma_ohiense_SECH_97.Nam2_CI_18class <- Histoplasma_ohiense_SECH_97.Nam2_CI_18class %>% mutate('Histoplasma ohiense SECH_97-Nam2_CI_18'= percentage)
Histoplasma_ohiense_SECH_97.Nam2_CI_18class1 <- Histoplasma_ohiense_SECH_97.Nam2_CI_18class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_98.Nam2_CI_30 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_98-Nam2_CI_30.RM/Histoplasma_ohiense_SECH_98-Nam2_CI_30.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_98.Nam2_CI_30) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_98.Nam2_CI_30 <- Histoplasma_ohiense_SECH_98.Nam2_CI_30 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_98.Nam2_CI_30,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_98.Nam2_CI_30.tsv")

Histoplasma_ohiense_SECH_98.Nam2_CI_30class<- Histoplasma_ohiense_SECH_98.Nam2_CI_30 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24259982)
write_tsv(Histoplasma_ohiense_SECH_98.Nam2_CI_30class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_98.Nam2_CI_30class.tsv")
Histoplasma_ohiense_SECH_98.Nam2_CI_30class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_98.Nam2_CI_30class.tsv")
Histoplasma_ohiense_SECH_98.Nam2_CI_30class <- Histoplasma_ohiense_SECH_98.Nam2_CI_30class %>% mutate('Histoplasma ohiense SECH_98.Nam2_CI_30'= percentage)
Histoplasma_ohiense_SECH_98.Nam2_CI_30class1 <- Histoplasma_ohiense_SECH_98.Nam2_CI_30class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiense_SECH_99.Nam2_CI_35 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiense_SECH_99-Nam2_CI_35.RM/Histoplasma_ohiense_SECH_99-Nam2_CI_35.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiense_SECH_99.Nam2_CI_35) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiense_SECH_99.Nam2_CI_35 <- Histoplasma_ohiense_SECH_99.Nam2_CI_35 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiense_SECH_99.Nam2_CI_35,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_99.Nam2_CI_35.tsv")

Histoplasma_ohiense_SECH_99.Nam2_CI_35class<- Histoplasma_ohiense_SECH_99.Nam2_CI_35 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29591009)
write_tsv(Histoplasma_ohiense_SECH_99.Nam2_CI_35class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_99.Nam2_CI_35class.tsv")
Histoplasma_ohiense_SECH_99.Nam2_CI_35class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_99.Nam2_CI_35class.tsv")
Histoplasma_ohiense_SECH_99.Nam2_CI_35class <- Histoplasma_ohiense_SECH_99.Nam2_CI_35class %>% mutate('Histoplasma ohiense SECH_99.Nam2_CI_35'= percentage)
Histoplasma_ohiense_SECH_99.Nam2_CI_35class1 <- Histoplasma_ohiense_SECH_99.Nam2_CI_35class %>% select(-c(n,total,percentage))

##
Histoplasma_ohiensis_Hc1986 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_ohiensis_Hc1986.RM/Histoplasma_ohiensis_Hc1986.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_ohiensis_Hc1986) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_ohiensis_Hc1986 <- Histoplasma_ohiensis_Hc1986 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_ohiensis_Hc1986,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiensis_Hc1986.tsv")

Histoplasma_ohiensis_Hc1986class<- Histoplasma_ohiensis_Hc1986 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/13148448)
write_tsv(Histoplasma_ohiensis_Hc1986class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiensis_Hc1986class.tsv")
Histoplasma_ohiensis_Hc1986class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiensis_Hc1986class.tsv")
Histoplasma_ohiensis_Hc1986class <- Histoplasma_ohiensis_Hc1986class %>% mutate('Histoplasma ohiensis Hc1986'= percentage)
Histoplasma_ohiensis_Hc1986class1 <- Histoplasma_ohiensis_Hc1986class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass <- Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass %>% mutate('Histoplasma capsulatum HISSP-CM7256'= percentage)
Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass1 <- Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_JB_01752Hc_01752class1 <- Histoplasma_capsulatum_JB_01752Hc_01752class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_JB_062775Hc_062775class1 <- Histoplasma_capsulatum_JB_062775Hc_062775class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_JB_073129Hc_073129class <- Histoplasma_capsulatum_JB_073129Hc_073129class %>% mutate('Histoplasma capsulatum JB_073129-Hc_073129'= percentage)
Histoplasma_capsulatum_JB_073129Hc_073129class1 <- Histoplasma_capsulatum_JB_073129Hc_073129class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_SECH_101Nam2_G184Aclass <- Histoplasma_capsulatum_SECH_101Nam2_G184Aclass %>% mutate('Histoplasma capsulatum senso_stricto SECH_101_G184AR'= percentage)
Histoplasma_capsulatum_SECH_101Nam2_G184Aclass1 <- Histoplasma_capsulatum_SECH_101Nam2_G184Aclass %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_SECH_109class1 <- Histoplasma_capsulatum_SECH_109class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_SECH_110class1 <- Histoplasma_capsulatum_SECH_110class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_HISSPFGTRO0285class <- Histoplasma_capsulatum_HISSPFGTRO0285class %>% mutate('Histoplasma capsulatum HISSP-FGTRO0285'= percentage)
Histoplasma_capsulatum_HISSPFGTRO0285class1 <- Histoplasma_capsulatum_HISSPFGTRO0285class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_HISSPFGAMA2041class1 <- Histoplasma_capsulatum_HISSPFGAMA2041class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_HISSPCM6408class1 <- Histoplasma_capsulatum_HISSPCM6408class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_HISSPB05821class1 <- Histoplasma_capsulatum_HISSPB05821class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_HISSP11571Belem1class1 <- Histoplasma_capsulatum_HISSP11571Belem1class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_HISSP1014Belem3class1 <- Histoplasma_capsulatum_HISSP1014Belem3class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_1517_p_17class1 <- Histoplasma_capsulatum_1517_p_17class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_256_P_18class <- Histoplasma_capsulatum_256_P_18class %>% mutate('Histoplasma capsulatum 256_P_18'= percentage)
Histoplasma_capsulatum_256_P_18class1 <- Histoplasma_capsulatum_256_P_18class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_316_p_10class1 <- Histoplasma_capsulatum_316_p_10class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_327_P_12class1 <- Histoplasma_capsulatum_327_P_12class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_343_p_18class <- Histoplasma_capsulatum_343_p_18class %>% mutate('Histoplasma capsulatum 343_p_18'= percentage)
Histoplasma_capsulatum_343_p_18class1 <- Histoplasma_capsulatum_343_p_18class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_388_p_11class <- Histoplasma_capsulatum_388_p_11class %>% mutate('Histoplasma capsulatum 388_p_11'= percentage)
Histoplasma_capsulatum_388_p_11class1 <- Histoplasma_capsulatum_388_p_11class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_ES2_83Zclass1 <- Histoplasma_capsulatum_ES2_83Zclass %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_ES2_85Zclass <- Histoplasma_capsulatum_ES2_85Zclass %>% mutate('Histoplasma capsulatum ES2_85Z'= percentage)
Histoplasma_capsulatum_ES2_85Zclass1 <- Histoplasma_capsulatum_ES2_85Zclass %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_ES2_86Zclass <- Histoplasma_capsulatum_ES2_86Zclass %>% mutate('Histoplasma capsulatum ES2_86Z'= percentage)
Histoplasma_capsulatum_ES2_86Zclass1 <- Histoplasma_capsulatum_ES2_86Zclass %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_ES2_87class <- Histoplasma_capsulatum_ES2_87class %>% mutate('Histoplasma capsulatum ES2_87'= percentage)
Histoplasma_capsulatum_ES2_87class1 <- Histoplasma_capsulatum_ES2_87class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_ES2_88Zclass <- Histoplasma_capsulatum_ES2_88Zclass %>% mutate('Histoplasma capsulatum ES2_88Z'= percentage)
Histoplasma_capsulatum_ES2_88Zclass1 <- Histoplasma_capsulatum_ES2_88Zclass %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_ES2_89Zclass <- Histoplasma_capsulatum_ES2_89Zclass %>% mutate('Histoplasma capsulatum ES2_89Z'= percentage)
Histoplasma_capsulatum_ES2_89Zclass1 <- Histoplasma_capsulatum_ES2_89Zclass %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_ES2_90Zclass <- Histoplasma_capsulatum_ES2_90Zclass %>% mutate('Histoplasma capsulatum ES2_90Z'= percentage)
Histoplasma_capsulatum_ES2_90Zclass1 <- Histoplasma_capsulatum_ES2_90Zclass %>% select(-c(n,total,percentage))

##
Histoplasma_capsulatum_SA15 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_EDTA/Histoplasma_capsulatum_SA15.RM/Histoplasma_capsulatum_SA15.scaffolds.fa.out",skip=3,col_names = F)
colnames(Histoplasma_capsulatum_SA15) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
Histoplasma_capsulatum_SA15 <- Histoplasma_capsulatum_SA15 %>% mutate(length=endq-beginq +1) 
write_tsv(Histoplasma_capsulatum_SA15,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SA15.tsv")

Histoplasma_capsulatum_SA15class<- Histoplasma_capsulatum_SA15 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28298332)
write_tsv(Histoplasma_capsulatum_SA15class,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SA15class.tsv")
Histoplasma_capsulatum_SA15class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SA15class.tsv")
Histoplasma_capsulatum_SA15class <- Histoplasma_capsulatum_SA15class %>% mutate('Histoplasma capsulatum SA15'= percentage)
Histoplasma_capsulatum_SA15class1 <- Histoplasma_capsulatum_SA15class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB053class <- Histoplasma_capsulatum_CB053class %>% mutate('Histoplasma capsulatum CB053'= percentage)
Histoplasma_capsulatum_CB053class1 <- Histoplasma_capsulatum_CB053class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB062class <- Histoplasma_capsulatum_CB062class %>% mutate('Histoplasma capsulatum CB062'= percentage)
Histoplasma_capsulatum_CB062class1 <- Histoplasma_capsulatum_CB062class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB063class1 <- Histoplasma_capsulatum_CB063class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB066class <- Histoplasma_capsulatum_CB066class %>% mutate('Histoplasma capsulatum CB066'= percentage)
Histoplasma_capsulatum_CB066class1 <- Histoplasma_capsulatum_CB066class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB180class <- Histoplasma_capsulatum_CB180class %>% mutate('Histoplasma capsulatum CB180'= percentage)
Histoplasma_capsulatum_CB180class1 <- Histoplasma_capsulatum_CB180class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB0522class <- Histoplasma_capsulatum_CB0522class %>% mutate('Histoplasma capsulatum CB052-2'= percentage)
Histoplasma_capsulatum_CB0522class1 <- Histoplasma_capsulatum_CB0522class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB0532class <- Histoplasma_capsulatum_CB0532class %>% mutate('Histoplasma capsulatum CB053-2'= percentage)
Histoplasma_capsulatum_CB0532class1 <- Histoplasma_capsulatum_CB0532class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB0552class1 <- Histoplasma_capsulatum_CB0552class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB0632class <- Histoplasma_capsulatum_CB0632class %>% mutate('Histoplasma capsulatum CB063-2'= percentage)
Histoplasma_capsulatum_CB0632class1 <- Histoplasma_capsulatum_CB0632class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB0642class <- Histoplasma_capsulatum_CB0642class %>% mutate('Histoplasma capsulatum CB064-2'= percentage)
Histoplasma_capsulatum_CB0642class1 <- Histoplasma_capsulatum_CB0642class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB0652class <- Histoplasma_capsulatum_CB0652class %>% mutate('Histoplasma capsulatum CB065-2'= percentage)
Histoplasma_capsulatum_CB0652class1 <- Histoplasma_capsulatum_CB0652class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB0662class <- Histoplasma_capsulatum_CB0662class %>% mutate('Histoplasma capsulatum CB066-2'= percentage)
Histoplasma_capsulatum_CB0662class1 <- Histoplasma_capsulatum_CB0662class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB1682class <- Histoplasma_capsulatum_CB1682class %>% mutate('Histoplasma capsulatum CB168-2'= percentage)
Histoplasma_capsulatum_CB1682class1 <- Histoplasma_capsulatum_CB1682class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB1742class <- Histoplasma_capsulatum_CB1742class %>% mutate('Histoplasma capsulatum CB174-2'= percentage)
Histoplasma_capsulatum_CB1742class1 <- Histoplasma_capsulatum_CB1742class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB1802class <- Histoplasma_capsulatum_CB1802class %>% mutate('Histoplasma capsulatum CB180-2'= percentage)
Histoplasma_capsulatum_CB1802class1 <- Histoplasma_capsulatum_CB1802class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB1862class <- Histoplasma_capsulatum_CB1862class %>% mutate('Histoplasma capsulatum CB186-2'= percentage)
Histoplasma_capsulatum_CB1862class1 <- Histoplasma_capsulatum_CB1862class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_CB1922class <- Histoplasma_capsulatum_CB1922class %>% mutate('Histoplasma capsulatum CB192-2'= percentage)
Histoplasma_capsulatum_CB1922class1 <- Histoplasma_capsulatum_CB1922class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class %>% mutate('Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S14'= percentage)
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class1 <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class %>% select(-c(n,total,percentage))

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
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class %>% mutate('Histoplasma capsulatum Dr_Anuradha_Fungal_WGS_S16'= percentage)
Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class1 <- Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_102Nam2_505class <- Histoplasma_mississippienseII_SECH_102Nam2_505class %>% mutate('Histoplasma mississippiense II SECH_102Nam2_505'= percentage)
Histoplasma_mississippienseII_SECH_102Nam2_505class1 <- Histoplasma_mississippienseII_SECH_102Nam2_505class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_81Nam1_WU24class <- Histoplasma_mississippienseII_SECH_81Nam1_WU24class %>% mutate('Histoplasma mississippiense II SECH_81-Nam1_WU24'= percentage)
Histoplasma_mississippienseII_SECH_81Nam1_WU24class1 <- Histoplasma_mississippienseII_SECH_81Nam1_WU24class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_83Nam1_CI_7class <- Histoplasma_mississippienseII_SECH_83Nam1_CI_7class %>% mutate('Histoplasma mississippiense II SECH_83-Nam1_CI_7'= percentage)
Histoplasma_mississippienseII_SECH_83Nam1_CI_7class1 <- Histoplasma_mississippienseII_SECH_83Nam1_CI_7class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_84Nam1_CI_19class <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19class %>% mutate('Histoplasma mississippiense II SECH_84-Nam1_CI_19'= percentage)
Histoplasma_mississippienseII_SECH_84Nam1_CI_19class1 <- Histoplasma_mississippienseII_SECH_84Nam1_CI_19class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_85Nam1_CI_22class <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_85-Nam1_CI_22class.tsv")
Histoplasma_mississippienseII_SECH_85Nam1_CI_22class <- Histoplasma_mississippienseII_SECH_85Nam1_CI_22class %>% mutate('Histoplasma mississippiense II SECH_85-Nam1_CI_22'= percentage)
Histoplasma_mississippienseII_SECH_85Nam1_CI_22class1 <- Histoplasma_mississippienseII_SECH_85Nam1_CI_22class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_86Nam1_CI_24class <- Histoplasma_mississippienseII_SECH_86Nam1_CI_24class %>% mutate('Histoplasma mississippiense II SECH_86-Nam1_CI_24'= percentage)
Histoplasma_mississippienseII_SECH_86Nam1_CI_24class1 <- Histoplasma_mississippienseII_SECH_86Nam1_CI_24class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_87Nam1_CI_42class <- Histoplasma_mississippienseII_SECH_87Nam1_CI_42class %>% mutate('Histoplasma mississippiense II SECH_87-Nam1_CI_42'= percentage)
Histoplasma_mississippienseII_SECH_87Nam1_CI_42class1 <- Histoplasma_mississippienseII_SECH_87Nam1_CI_42class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_88Nam1_CI_43class <- Histoplasma_mississippienseII_SECH_88Nam1_CI_43class %>% mutate('Histoplasma mississippiense II SECH_88-Nam1_CI_43'= percentage)
Histoplasma_mississippienseII_SECH_88Nam1_CI_43class1 <- Histoplasma_mississippienseII_SECH_88Nam1_CI_43class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class <- Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class %>% mutate('Histoplasma mississippiense II SECH_89-Nam1_UCLA-531'= percentage)
Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class1 <- Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_HCCI_43class <- Histoplasma_mississippienseII_HCCI_43class %>% mutate('Histoplasma mississippiense II HCCI_43'= percentage)
Histoplasma_mississippienseII_HCCI_43class1 <- Histoplasma_mississippienseII_HCCI_43class %>% select(-c(n,total,percentage))

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
Histoplasma_mississippienseII_HCCI_19class <- Histoplasma_mississippienseII_HCCI_19class %>% mutate('Histoplasma mississippiense II HCCI_19'= percentage)
Histoplasma_mississippienseII_HCCI_19class1 <- Histoplasma_mississippienseII_HCCI_19class %>% select(-c(n,total,percentage))

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
Histoplasma_suramericanum_SECH_104Nam2_27_14class <- Histoplasma_suramericanum_SECH_104Nam2_27_14class %>% mutate('Histoplasma suramericanum SECH_104-Nam2_27_14'= percentage)
Histoplasma_suramericanum_SECH_104Nam2_27_14class1 <- Histoplasma_suramericanum_SECH_104Nam2_27_14class %>% select(-c(n,total,percentage))

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
Histoplasma_suramericanum_SECH_105Nam2_21_14class <- Histoplasma_suramericanum_SECH_105Nam2_21_14class %>% mutate('Histoplasma suramericanum SECH_105-Nam2_21_14'= percentage)
Histoplasma_suramericanum_SECH_105Nam2_21_14class1 <- Histoplasma_suramericanum_SECH_105Nam2_21_14class %>% select(-c(n,total,percentage))

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
Histoplasma_suramericanum_HC27_14class <- Histoplasma_suramericanum_HC27_14class %>% mutate('Histoplasma suramericanum HC27_14'= percentage)
Histoplasma_suramericanum_HC27_14class1 <- Histoplasma_suramericanum_HC27_14class %>% select(-c(n,total,percentage))

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

library(plyr)
library(dplyr)

#Combine results into one large table
Total_table<-join_all(list(Histoplasma_ohiense_HCCI_17class1, Histoplasma_ohiense_HCCI_6class1, Histoplasma_ohiense_HCG217Bclass1, Histoplasma_ohiense_SECH_82.Nam1_CI_4class1, Histoplasma_ohiense_SECH_91.Nam2_G217Bclass1, Histoplasma_ohiense_SECH_92.Nam2_G222Bclass1, Histoplasma_ohiense_SECH_93.Nam2_CI_6class1, Histoplasma_ohiense_SECH_94.Nam2_CI_9class1, Histoplasma_ohiense_SECH_95.Nam2_CI_10class1, Histoplasma_ohiense_SECH_96.Nam2_CI_17class1, Histoplasma_ohiense_SECH_97.Nam2_CI_18class1, Histoplasma_ohiense_SECH_98.Nam2_CI_30class1, Histoplasma_ohiense_SECH_99.Nam2_CI_35class1, Histoplasma_ohiensis_Hc1986class1, Histoplasma_capsulatum_19VMG15class1, Histoplasma_capsulatum_HISSPCM7256xxCLSENxxxx036BBclass1, Histoplasma_capsulatum_Histo485P20class1, Histoplasma_capsulatum_JB_01752Hc_01752class1, Histoplasma_capsulatum_JB_021091Hc_021091class1, Histoplasma_capsulatum_JB_031837Hc_031837class1, Histoplasma_capsulatum_JB_042430Hc_042430class1, Histoplasma_capsulatum_JB_062632Hc_062632class1, Histoplasma_capsulatum_JB_062775Hc_062775class1, Histoplasma_capsulatum_JB_073129Hc_073129class1, Histoplasma_capsulatum_JB_083285_2Hc_083285_2class1, Histoplasma_capsulatum_SECH_101Nam2_G184Aclass1, Histoplasma_capsulatum_SECH_107mis_Hc_duboisiiBclass1, Histoplasma_capsulatum_SECH_109class1, Histoplasma_capsulatum_SECH_110class1, Histoplasma_capsulatum_HCH143class1, Histoplasma_capsulatum_HCG186Aclass1, Histoplasma_capsulatum_NACVFR_Histo_HC1070058_2class1, Histoplasma_capsulatum_HISSPFGTRO0285class1, Histoplasma_capsulatum_HISSPFGPSO2043class1, Histoplasma_capsulatum_HISSPFGPIE2055class1, Histoplasma_capsulatum_HISSPFGPIA2052class1,Histoplasma_capsulatum_HISSPFGPERS2034class1, Histoplasma_capsulatum_HISSPFGMAR2044class1, Histoplasma_capsulatum_HISSPFGLIN2055class1, Histoplasma_capsulatum_HISSPFGJOS2044class1, Histoplasma_capsulatum_HISSPFGGRE2022class1, Histoplasma_capsulatum_HISSPFGFIN2028class1, Histoplasma_capsulatum_HISSPFGFER2036class1, Histoplasma_capsulatum_HISSPFGFAR0189class1, Histoplasma_capsulatum_HISSPFGFAN2059class1, Histoplasma_capsulatum_HISSPFGBON2001class1, Histoplasma_capsulatum_HISSPFGBIK2051class1, Histoplasma_capsulatum_HISSPFGAMA2041class1, Histoplasma_capsulatum_HISSPCM6408class1, Histoplasma_capsulatum_HISSPCM6015class1, Histoplasma_capsulatum_HISSPB05821class1, Histoplasma_capsulatum_HISSP11571Belem1class1, Histoplasma_capsulatum_HISSP1014Belem3class1, Histoplasma_capsulatum_07_12RJclass1, Histoplasma_capsulatum_104_p_06class1, Histoplasma_capsulatum_104_P_19class1, Histoplasma_capsulatum_117_p_12class1, Histoplasma_capsulatum_122_p_10_Bclass1, Histoplasma_capsulatum_136_P_07class1, Histoplasma_capsulatum_144_p_08class1, Histoplasma_capsulatum_1517_p_17class1, Histoplasma_capsulatum_256_P_18class1, Histoplasma_capsulatum_CB063class1, Histoplasma_capsulatum_CB066class1, Histoplasma_capsulatum_CB180class1, Histoplasma_capsulatum_CB0522class1, Histoplasma_capsulatum_CB0532class1, Histoplasma_capsulatum_CB0552class1, Histoplasma_capsulatum_CB0622class1, Histoplasma_capsulatum_CB0632class1, Histoplasma_capsulatum_CB0642class1, Histoplasma_capsulatum_CB0652class1, Histoplasma_capsulatum_CB0662class1, Histoplasma_capsulatum_CB1682class1, Histoplasma_capsulatum_CB1742class1, Histoplasma_capsulatum_CB1802class1, Histoplasma_capsulatum_CB1862class1, Histoplasma_capsulatum_CB1922class1, Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class1, Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class1, Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class1, Histoplasma_mississippienseII_SECH_102Nam2_505class1, Histoplasma_mississippienseII_SECH_81Nam1_WU24class1, Histoplasma_mississippienseII_SECH_83Nam1_CI_7class1, Histoplasma_mississippienseII_SECH_84Nam1_CI_19class1, Histoplasma_mississippienseII_SECH_85Nam1_CI_22class1, Histoplasma_mississippienseII_SECH_86Nam1_CI_24class1, Histoplasma_mississippienseII_SECH_87Nam1_CI_42class1, Histoplasma_mississippienseII_SECH_88Nam1_CI_43class1),by="class", type="full")
Total_table2<-join_all(list(Histoplasma_capsulatum_SECH_100Nam2_G186Aclass1,Histoplasma_mississippienseII_SECH_89Nam1_UCLA531class1,Histoplasma_mississippienseII_SECH_90Nam1_DOWNSclass1, Histoplasma_mississippienseII_HCWU24class1, Histoplasma_mississippienseII_HCCI_43class1, Histoplasma_mississippienseII_HCCI_19class1, Histoplasma_suramericanum_SECH_103Nam2_3_11Gclass1, Histoplasma_suramericanum_SECH_104Nam2_27_14class1, Histoplasma_suramericanum_SECH_105Nam2_21_14class1, Histoplasma_suramericanum_HC27_14class1, Histoplasma_suramericanum_HC21_14class1, Histoplasma_capsulatum_CB053class1, Histoplasma_capsulatum_CB062class1, Histoplasma_capsulatum_316_p_10class1, Histoplasma_capsulatum_327_P_12class1, Histoplasma_capsulatum_343_p_18class1, Histoplasma_capsulatum_388_p_11class1, Histoplasma_capsulatum_ES2_83Zclass1, Histoplasma_capsulatum_ES2_85Zclass1,Histoplasma_capsulatum_ES2_86Zclass1,Histoplasma_capsulatum_ES2_87class1, Histoplasma_capsulatum_ES2_88Zclass1,Histoplasma_capsulatum_ES2_89Zclass1,Histoplasma_capsulatum_ES2_90Zclass1, Histoplasma_capsulatum_SA15class1, Blastomyces_parvus_UAMH130class1),by="class", type="full")
Table_ALL<-join_all(list(Total_table, Total_table2), by="class", type="full")
write_tsv(Table_ALL,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/All_Histo.tsv")

#Manipulated the table outputs to seven types of TEs instead of the 34 initially described, then loaded these back in.
Complete_TE_complex <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Percentage_95_complex.tsv")
Complete_TE_simple <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Percentage_95_clean.tsv")
##Building the stacked barplot
#Found this description to follow https://zebrabi.com/guide/how-to-customize-stacked-bar-chart-in-r-ggplot2/#:~:text=Creating%20a%20basic%20stacked%20bar%20chart%20using%20ggplot2%20in%20R&text=For%20instance%2C%20we%20can%20change,bars%20using%20the%20width%20parameter

library(dplyr)
library(tidyr)
library(stringr)
library(reshape2) 
library(reshape)
library(RColorBrewer)
library(ggplot2)


Complete_TE_complex <- as.data.frame(Complete_TE_complex)
Complete_TE_simple <- as.data.frame(Complete_TE_simple)
Complete_TE_complex[is.na(Complete_TE_complex)] <- 0
Complete_TE_simple[is.na(Complete_TE_simple)] <- 0

Complete_TE_melt_complex<- melt(Complete_TE_complex, id.vars="class", 
                            measure.vars=c("Histoplasma capsulatum 19VMG15","Histoplasma capsulatum HISSP-CM7256",
                                           "Histoplasma capsulatum Histo-485P20","Histoplasma capsulatum JB_01752-Hc_01752",
                                           "Histoplasma capsulatum JB_021091-Hc_021091","Histoplasma capsulatum JB_031837-Hc_031837",
                                           "Histoplasma capsulatum JB_042430-Hc_042430","Histoplasma capsulatum JB_062632-Hc_062632",
                                           "Histoplasma capsulatum JB_062775-Hc_062775","Histoplasma capsulatum JB_073129-Hc_073129",
                                           "Histoplasma capsulatum JB_083285_2-Hc_083285_2","Histoplasma capsulatum SECH_100-G186A",
                                           "Histoplasma capsulatum senso_stricto SECH_101_G184AR","Histoplasma mississippiense SECH_102-505",
                                           "Histoplasma capsulatum SECH_103-3_11G","Histoplasma capsulatum SECH_104-27_14",
                                           "Histoplasma capsulatumSECH_105-21_14","Histoplasma capsulatum SECH_107-duboisii-B",
                                           "Histoplasma capsulatum SECH_109","Histoplasma capsulatum SECH_110",
                                           "Histoplasma mississippiense SECH_81-WU24","Histoplasma mississippiense SECH_83-CI_7",
                                           "Histoplasma mississippiense SECH_84-CI_19","Histoplasma mississippiense SECH_85-CI_22",
                                           "Histoplasma mississippiense SECH_86-CI_24","Histoplasma mississippiense SECH_87-CI_42",
                                           "Histoplasma mississippiense SECH_88-CI_43","Histoplasma mississippiense SECH_89-UCLA-531",
                                           "Histoplasma mississippiense SECH_90-DOWNS","Histoplasma ohiense SECH_82-CI_4",
                                           "Histoplasma ohiense SECH_91-G217B","Histoplasma ohiense SECH_92-G222B","Histoplasma ohiense SECH_93-CI_6",
                                           "Histoplasma ohiense SECH_94-CI_9","Histoplasma ohiense SECH_95-CI_10","Histoplasma ohiense SECH_96-CI_17",
                                           "Histoplasma ohiense SECH_97-CI_18","Histoplasma ohiense SECH_98-CI_30","Histoplasma ohiense SECH_99-CI_35",
                                           "Histoplasma capsulatum HCH143","Histoplasma capsulatum NACVFR_Histo_HC1070058_2",
                                           "Histoplasma capsulatum HISSP-FGTRO0285","Histoplasma capsulatum HISSP-FGPSO2043",
                                           "Histoplasma capsulatum HISSP-FGPIE2055","Histoplasma capsulatum HISSP-FGPIA2052",
                                           "Histoplasma capsulatum HISSP-FGPERS2034","Histoplasma capsulatum HISSP-FGMAR2044",
                                           "Histoplasma capsulatum HISSP-FGLIN2055","Histoplasma capsulatum HISSP-FGJOS2044",
                                           "Histoplasma capsulatum HISSP-FGFIN2028","Histoplasma capsulatum HHISSP-FGFER2036",
                                           "Histoplasma capsulatum HISSP-FGFAR0189","Histoplasma capsulatum HISSP-FGFAN2059",
                                           "Histoplasma capsulatum HISSP-FGBON2001","Histoplasma capsulatum HISSP-FGBIK205",
                                           "Histoplasma capsulatum HISSP-FGAMA2041","Histoplasma capsulatum HISSP-CM6408",
                                           "Histoplasma capsulatum HISSP-CM6015","Histoplasma capsulatum HISSP-B05821",
                                           "Histoplasma capsulatum HISSP-11571-Belem1","Histoplasma capsulatum HISSP-1014-Belem3",
                                           "Histoplasma capsulatum 07_12-RJ","Histoplasma capsulatum 104_p_06","Histoplasma capsulatum 117_p_12",
                                           "Histoplasma capsulatum 122_p_10_B","Histoplasma capsulatum 136_P_07","Histoplasma capsulatum 144_p_08",
                                           "Histoplasma capsulatum 1517_p_17","Histoplasma capsulatum 256_P_18","Histoplasma capsulatum 316_p_10",
                                           "Histoplasma capsulatum 327_P_12","Histoplasma capsulatum 343_p_18","Histoplasma capsulatum 388_p_11",
                                           "Histoplasma capsulatum ES2_83Z","Histoplasma capsulatum ES2_85Z","Histoplasma capsulatum ES2_86Z",
                                           "Histoplasma capsulatum ES2_88Z","Histoplasma capsulatum ES2_89Z","Histoplasma capsulatum ES2_90Z",
                                           "Histoplasma capsulatum SA15","Histoplasma capsulatum CB053-2","Histoplasma capsulatum CB055-2",
                                           "Histoplasma capsulatum CB062-2","Histoplasma capsulatum CB063-2","Histoplasma capsulatum CB064-2",
                                           "Histoplasma capsulatum CB065-2","Histoplasma capsulatum CB066-2","Histoplasma capsulatum CB168-2",
                                           "Histoplasma capsulatum CB174-2","Histoplasma capsulatum CB180-2","Histoplasma capsulatum CB186-2",
                                           "Histoplasma capsulatum CB192-2","Histoplasma capsulatum WGS_S11","Histoplasma capsulatum WGS_S14",
                                           "Histoplasma capsulatum WGS_S16","Blastomyces parvus UAMH130"))

write_tsv(Complete_TE_melt_complex,"/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/FULL_Table.tsv")

#Saving results to a PDF file
pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/Stacked_barplot_COMPLEX.pdf', width=15, height=20)

ggplot(Complete_TE_melt_complex, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip(ylim=c()) + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of Repetative Elements per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') 

dev.off()

png(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/Stacked_barplot_COMPLEX.png',width=1500, height=1900, res=100)

ggplot(Complete_TE_melt_complex, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of TE per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') 

dev.off()


Complete_TE_melt_short<- melt(Complete_TE_simple, id.vars="class", 
                              measure.vars=c("Histoplasma capsulatum 19VMG15","Histoplasma capsulatum HISSP-CM7256",
                                          "Histoplasma capsulatum Histo-485P20","Histoplasma capsulatum JB_01752-Hc_01752",
                                          "Histoplasma capsulatum JB_021091-Hc_021091","Histoplasma capsulatum JB_031837-Hc_031837",
                                          "Histoplasma capsulatum JB_042430-Hc_042430","Histoplasma capsulatum JB_062632-Hc_062632",
                                          "Histoplasma capsulatum JB_062775-Hc_062775","Histoplasma capsulatum JB_073129-Hc_073129",
                                          "Histoplasma capsulatum JB_083285_2-Hc_083285_2","Histoplasma capsulatum SECH_100-G186A",
                                          "Histoplasma capsulatum senso_stricto SECH_101_G184AR","Histoplasma mississippiense SECH_102-505",
                                          "Histoplasma capsulatum SECH_103-3_11G","Histoplasma capsulatum SECH_104-27_14",
                                          "Histoplasma capsulatum SECH_105-21_14","Histoplasma capsulatum SECH_107-duboisii-B",
                                          "Histoplasma capsulatum SECH_109","Histoplasma capsulatum SECH_110",
                                          "Histoplasma mississippiense SECH_81-WU24","Histoplasma mississippiense SECH_83-CI_7",
                                          "Histoplasma mississippiense SECH_84-CI_19","Histoplasma mississippiense SECH_85-CI_22",
                                          "Histoplasma mississippiense SECH_86-CI_24","Histoplasma mississippiense SECH_87-CI_42",
                                          "Histoplasma mississippiense SECH_88-CI_43","Histoplasma mississippiense SECH_89-UCLA-531",
                                          "Histoplasma mississippiense SECH_90-DOWNS","Histoplasma ohiense SECH_82-CI_4",
                                          "Histoplasma ohiense SECH_91-G217B","Histoplasma ohiense SECH_92-G222B","Histoplasma ohiense SECH_93-CI_6",
                                          "Histoplasma ohiense SECH_94-CI_9","Histoplasma ohiense SECH_95-CI_10","Histoplasma ohiense SECH_96-CI_17",
                                          "Histoplasma ohiense SECH_97-CI_18","Histoplasma ohiense SECH_98-CI_30","Histoplasma ohiense SECH_99-CI_35",
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
                                          "Histoplasma capsulatum 07_12-RJ","Histoplasma capsulatum 104_p_06","Histoplasma capsulatum 117_p_12",
                                          "Histoplasma capsulatum 122_p_10_B","Histoplasma capsulatum 136_P_07","Histoplasma capsulatum 144_p_08",
                                          "Histoplasma capsulatum 1517_p_17","Histoplasma capsulatum 256_P_18","Histoplasma capsulatum 316_p_10",
                                          "Histoplasma capsulatum 327_P_12","Histoplasma capsulatum 343_p_18","Histoplasma capsulatum 388_p_11",
                                          "Histoplasma capsulatum ES2_83Z","Histoplasma capsulatum ES2_85Z","Histoplasma capsulatum ES2_86Z",
                                          "Histoplasma capsulatum ES2_88Z","Histoplasma capsulatum ES2_89Z","Histoplasma capsulatum ES2_90Z",
                                          "Histoplasma capsulatum SA15","Histoplasma capsulatum CB053","Histoplasma capsulatum CB055",
                                          "Histoplasma capsulatum CB062","Histoplasma capsulatum CB063","Histoplasma capsulatum CB064",
                                          "Histoplasma capsulatum CB065","Histoplasma capsulatum CB066","Histoplasma capsulatum CB168",
                                          "Histoplasma capsulatum CB174","Histoplasma capsulatum CB180","Histoplasma capsulatum CB186",
                                          "Histoplasma capsulatum CB192","Histoplasma capsulatum WGS_S11","Histoplasma capsulatum WGS_S14",
                                          "Histoplasma capsulatum WGS_S16","Blastomyces parvus UAMH130"))

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/Stacked_barplot_short.pdf', width=15, height=20)

ggplot(Complete_TE_melt_short, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of TE per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired") 

dev.off()

png(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/Stacked_barplot_short.png',width=1500, height=1900, res=100)

ggplot(Complete_TE_melt_short, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of TE per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired") 

dev.off()


############~AMINO ACID TREE~###################################

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
p_AAB <- ggtree(tree_AAB) + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.008)
ggsave(filename = "/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/species_tree_with_numbers.png", plot = p_AAB,width = 9, height = 6, units = c("in"), dpi = 300, scale =2, device = "png")

new_seq <- tree_AAB$tip.label
d = data.frame(new_seq)
d$Species = 0
d$Species[grep("UAMH130", d$new_seq)] = "Blastomyces parvus UAMH130"
d$Species[grep("19VMG-15", d$new_seq)] = "Histoplasma capsulatum 19VMG15"
d$Species[grep("07_12-RJ", d$new_seq)] = "Histoplasma capsulatum 07_12-RJ"
d$Species[grep("HISSP-11571-Belem1", d$new_seq)] = "Histoplasma capsulatum HISSP-11571-Belem1"
d$Species[grep("HISSP-FGFER2036", d$new_seq)] = "Histoplasma capsulatum HISSP-FGFER2036"
d$Species[grep("HISSP-CM6015", d$new_seq)] = "Histoplasma capsulatum HISSP-CM6015"
d$Species[grep("HISSP-FGAMA2041", d$new_seq)] = "Histoplasma capsulatum HISSP-FGAMA2041"
d$Species[grep("HISSP-FGLIN2055", d$new_seq)] = "Histoplasma capsulatum HISSP-FGLIN2055"
d$Species[grep("HISSP-FGFAN2059", d$new_seq)] = "Histoplasma capsulatum HISSP-FGFAN2059"
d$Species[grep("HISSP-FGPIA2052", d$new_seq)] = "Histoplasma capsulatum HISSP-FGPIA2052"
d$Species[grep("HISSP-FGPIE2055", d$new_seq)] = "Histoplasma capsulatum HISSP-FGPIE2055"
d$Species[grep("HISSP-FGBON2001", d$new_seq)] = "Histoplasma capsulatum HISSP-FGBON2001"
d$Species[grep("HISSP-FGBIK2051", d$new_seq)] = "Histoplasma capsulatum HISSP-FGBIK2051"
d$Species[grep("HISSP-FGJOS2044", d$new_seq)] = "Histoplasma capsulatum HISSP-FGJOS2044"
d$Species[grep("JB_062632-Hc_062632", d$new_seq)] = "Histoplasma capsulatum JB_062632-Hc_062632" 
d$Species[grep("NACVFR_Histo_HC1070058_2", d$new_seq)] = "Histoplasma capsulatum NACVFR_Histo_HC1070058"
d$Species[grep("SECH_103", d$new_seq)] = "Histoplasma capsulatum SECH_103_3_11G"
d$Species[grep("SECH_105", d$new_seq)] = "Histoplasma capsulatum SECH_105_21_14"
d$Species[grep("HISSP-FGPERS2034", d$new_seq)] = "Histoplasma capsulatum HISSP-FGPERS2034"
d$Species[grep("HISSP-B05821", d$new_seq)] = "Histoplasma capsulatum HISSP-B05821"
d$Species[grep("JB_031837-Hc_031837", d$new_seq)] = "Histoplasma capsulatum JB_031837-Hc_031837"
d$Species[grep("JB_042430-Hc_042430", d$new_seq)] = "Histoplasma capsulatum JB_042430-Hc_042430"
d$Species[grep("JB_062775-Hc_062775", d$new_seq)] = "Histoplasma capsulatum JB_062775-Hc_062775"
d$Species[grep("JB_083285_2-Hc_083285_2", d$new_seq)] = "Histoplasma capsulatum JB_083285_2-Hc_083285_2"
d$Species[grep("CB053-2", d$new_seq)] = "Histoplasma capsulatum CB053" 
d$Species[grep("CB055-2", d$new_seq)] = "Histoplasma capsulatum CB055"
d$Species[grep("SECH_101", d$new_seq)] = "Histoplasma capsulatum senso_stricto SECH_101_G184AR"
d$Species[grep("HISSP-CM6408", d$new_seq)] = "Histoplasma capsulatum HISSP-CM6408"
d$Species[grep("HISSP-1014-Belem3", d$new_seq)] = "Histoplasma capsulatum HISSP-1014-Belem3"
d$Species[grep("ES2_83Z", d$new_seq)] = "Histoplasma capsulatum ES2_83Z"
d$Species[grep("ES2_90Z", d$new_seq)] = "Histoplasma capsulatum ES2_90Z"
d$Species[grep("ES2_86Z", d$new_seq)] = "Histoplasma capsulatum ES2_86Z"
d$Species[grep("ES2_88Z", d$new_seq)] = "Histoplasma capsulatum ES2_88Z"
d$Species[grep("HCH143", d$new_seq)] = "Histoplasma capsulatum HCH143"
d$Species[grep("SECH_107", d$new_seq)] = "Histoplasma capsulatum SECH_107-duboisii-B" 
d$Species[grep("ES2_85Z", d$new_seq)] = "Histoplasma capsulatum ES2_85Z"
d$Species[grep("HISSP-CM7256", d$new_seq)] = "Histoplasma capsulatum HISSP-CM7256"
d$Species[grep("ES2_89Z", d$new_seq)] = "Histoplasma capsulatum ES2_89Z"
d$Species[grep("SA15", d$new_seq)] = "Histoplasma capsulatum SA15"
d$Species[grep("SECH_98", d$new_seq)] = "Histoplasma ohiense SECH_98-CI_30"
d$Species[grep("SECH_110", d$new_seq)] = "Histoplasma capsulatum SECH_110"
d$Species[grep("G217B", d$new_seq)] = "Histoplasma ohiense SECH_91-G217B"
d$Species[grep("SECH_92", d$new_seq)] = "Histoplasma ohiense SECH_92-G222B"
d$Species[grep("Hc1986", d$new_seq)] = "Histoplasma ohiense Hc1986"
d$Species[grep("SECH_82", d$new_seq)] = "Histoplasma ohiense SECH_82-CI_4"
d$Species[grep("SECH_95", d$new_seq)] = "Histoplasma ohiense SECH_95-CI_10"
d$Species[grep("SECH_94", d$new_seq)] = "Histoplasma ohiense SECH_94-CI_9"
d$Species[grep("SECH_96", d$new_seq)] = "Histoplasma ohiense SECH_96-CI_17"
d$Species[grep("SECH_93", d$new_seq)] = "Histoplasma ohiense SECH_93-CI_6"
d$Species[grep("SECH_97", d$new_seq)] = "Histoplasma ohiense SECH_97-CI_18"
d$Species[grep("SECH_99", d$new_seq)] = "Histoplasma ohiense SECH_99-CI_35"
d$Species[grep("SECH_104", d$new_seq)] = "Histoplasma capsulatum SECH_104-27_14" 
d$Species[grep("CB174-2", d$new_seq)] = "Histoplasma capsulatum CB174"
d$Species[grep("SECH_109", d$new_seq)] = "Histoplasma capsulatum SECH_109"
d$Species[grep("CB062-2", d$new_seq)] = "Histoplasma capsulatum CB062"
d$Species[grep("CB065-2", d$new_seq)] = "Histoplasma capsulatum CB065"
d$Species[grep("CB192-2", d$new_seq)] = "Histoplasma capsulatum CB192"
d$Species[grep("CB064-2", d$new_seq)] = "Histoplasma capsulatum CB064"
d$Species[grep("CB180-2", d$new_seq)] = "Histoplasma capsulatum CB180"
d$Species[grep("CB063-2", d$new_seq)] = "Histoplasma capsulatum CB063"
d$Species[grep("CB186-2", d$new_seq)] = "Histoplasma capsulatum CB186" 
d$Species[grep("SECH_81", d$new_seq)] = "Histoplasma mississippiense SECH_81-WU24"
d$Species[grep("CB066-2", d$new_seq)] = "Histoplasma capsulatum CB066"
d$Species[grep("SECH_88", d$new_seq)] = "Histoplasma mississippiense SECH_88-CI_43"
d$Species[grep("SECH_86", d$new_seq)] = "Histoplasma mississippiense SECH_86-CI_24"
d$Species[grep("CB168-2", d$new_seq)] = "Histoplasma capsulatum CB168"
d$Species[grep("SECH_85", d$new_seq)] = "Histoplasma mississippiense SECH_85-CI_22"
d$Species[grep("SECH_83", d$new_seq)] = "Histoplasma mississippiense SECH_83-CI_7" 
d$Species[grep("SECH_87", d$new_seq)] = "Histoplasma mississippiense SECH_87-CI_42"
d$Species[grep("SECH_90", d$new_seq)] = "Histoplasma mississippiense SECH_90-DOWNS"
d$Species[grep("SECH_102", d$new_seq)] = "Histoplasma mississippiense SECH_102_505"
d$Species[grep("SECH_100", d$new_seq)] = "Histoplasma capsulatum SECH_100-G186A"
d$Species[grep("SECH_84", d$new_seq)] = "Histoplasma mississippiense SECH_84-CI_19"
d$Species[grep("SECH_89", d$new_seq)] = "Histoplasma mississippiense SECH_89-UCLA-531"
d$Species[grep("JB_01752-Hc_01752", d$new_seq)] = "Histoplasma capsulatum JB_01752-Hc_01752"
d$Species[grep("104_p_06", d$new_seq)] = "Histoplasma capsulatum 104_p_06"
d$Species[grep("1517_p_17", d$new_seq)] = "Histoplasma capsulatum 1517_p_17"
d$Species[grep("117_p_12", d$new_seq)] = "Histoplasma capsulatum 117_p_12"
d$Species[grep("388_p_11", d$new_seq)] = "Histoplasma capsulatum 388_p_11" 
d$Species[grep("144_p_08", d$new_seq)] = "Histoplasma capsulatum 144_p_08"
d$Species[grep("327_P_12", d$new_seq)] = "Histoplasma capsulatum 327_P_12"
d$Species[grep("316_p_10", d$new_seq)] = "Histoplasma capsulatum 316_p_10"
d$Species[grep("Dr_Anuradha_Fungal_WGS_S14", d$new_seq)] = "Histoplasma capsulatum WGS_S14"
d$Species[grep("Histo-485P20", d$new_seq)] = "Histoplasma capsulatum Histo-485P20"
d$Species[grep("136_P_07", d$new_seq)] = "Histoplasma capsulatum 136_P_07"
d$Species[grep("343_p_18", d$new_seq)] = "Histoplasma capsulatum 343_p_18"
d$Species[grep("Dr_Anuradha_Fungal_WGS_S11", d$new_seq)] = "Histoplasma capsulatum WGS_S11"
d$Species[grep("Dr_Anuradha_Fungal_WGS_S16", d$new_seq)] = "Histoplasma capsulatum WGS_S16"
d$Species[grep("122_p_10_B", d$new_seq)] = "Histoplasma capsulatum 122_p_10_B"
d$Species[grep("256_P_18", d$new_seq)] = "Histoplasma capsulatum 256_P_18" 
d$Species[grep("JB_021091-Hc_021091", d$new_seq)] = "Histoplasma capsulatum JB_021091-Hc_021091"
d$Species[grep("JB_073129-Hc_073129", d$new_seq)] = "Histoplasma capsulatum JB_073129-Hc_073129"
d$Species[grep("HISSP-FGFAR0189", d$new_seq)] = "Histoplasma capsulatum HISSP-FGFAR0189"
d$Species[grep("HISSP-FGFIN2028", d$new_seq)] = "Histoplasma capsulatum HISSP-FGFIN2028"
d$Species[grep("HISSP-FGPSO2043", d$new_seq)] = "Histoplasma capsulatum HISSP-FGPSO2043"
d$Species[grep("HISSP-FGMAR2044", d$new_seq)] = "Histoplasma capsulatum HISSP-FGMAR2044"
d$Species[grep("HISSP-FGTRO0285", d$new_seq)] = "Histoplasma capsulatum HISSP-FGTRO0285"
tree_AAB$tip.label <- d$Species

# Visualize the tree
p_AABBB1 <- ggtree(tree_AAB, layout = "rectangular") + geom_tiplab(size = 3, fontface = 3, align=TRUE) +
  theme_tree2() + theme(legend.position = "none") 
print(p_AABBB1)
# Save the tree as a PDF
ggsave("/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/AAB_ultrametric_tree.pdf", p_AABBB1, width = 20, height = 5, dpi = 150)

#Sanity check!
colnames(Complete_TE_melt_short)[which(colnames(Complete_TE_melt_short)  %in% tree_AAB$tip.label==F)]

colnames(Complete_TE_melt_short)
Complete_TE_melt_short2 <- Complete_TE_melt_short[,c(2,1,3)]

#pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_fan_Stacked_barplot_AAB-Star.pdf', width=15, height=15)

p_AABBB2 <- p_AABBB1 + geom_facet(panel = "TE Abundance", data = Complete_TE_melt_short2, geom = geom_col, 
                                aes(x = value, color = class, 
                                    fill = class), orientation = 'y', width = .6) + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.05) +
  theme_tree2(legend.position=c(.95, .25)) + xlim_expand(c(0,0.8), 'TE Abundance') + xlim_expand(c(0, 4), 'Tree') + geom_cladelabel(node=140, label="NAm1", color="red", offset = 1.5) +
  geom_cladelabel(node=160, label="NAm2", color="#808000", offset = 1.5) + geom_cladelabel(node=172, label="India", color="blue", offset = 1.5) + geom_cladelabel(node=113, label="Africa", color="#FFA500", offset = 1.5) +
  geom_cladelabel(node=105, label="LAmA", color="purple", offset = 1.5) + geom_cladelabel(node=121, label="LAmA", color="purple", offset = 1.5) + geom_cladelabel(node=187, label="LAmA", color="purple", offset = 1.5) +
  geom_cladelabel(node=191, label="Africa", color="#FFA500", offset = 1.5) + geom_cladelabel(node=137, label="LAmA", color="purple", offset = 1.5) + geom_cladelabel(node=1, label="Root", color="green", offset = 1.5)

#dev.off()

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/Tree_Stacked_barplotAAB-LINES.pdf', width=15, height=15)

facet_widths(p_AABBB2, widths = c(1, 0.40))

dev.off()

png(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/PLOTS/Tree_fan_Stacked_barplotAAB-LINES.png', width=1700, height=1500, res=100)

facet_widths(p_AABBB2, widths = c(1, 0.40))

dev.off()

#if you need to figure out what the nodes are in your tree
NODES_TREE<-ggtree(tree_AAB) + geom_text(aes(label=node), hjust=-.3)

library(phytools)

#RUN significance tests on tree_AAB
colnames(Complete_TE_simple) = gsub("[-]", "_", colnames(Complete_TE_simple))
tree_AAB$tip.label = gsub("[-]", "_", tree_AAB$tip.label)
all(colnames(Complete_TE_simple)  %in% tree_AAB$tip.label)
colnames(Complete_TE_simple)[which(colnames(Complete_TE_simple)  %in% tree_AAB$tip.label==F)]
#Just first column doesn't agree which is class

#Test TE for phylogenetic signal - but phylosig requires named vector for data
row_DNA = which(Complete_TE_simple$class=="DNA")
data.vec_DNA = as.numeric(Complete_TE_simple[row_DNA, ]); names(data.vec_DNA) = colnames(Complete_TE_simple)

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
row_RNA = which(Complete_TE_simple$class=="RNA")
data.vec_RNA = as.numeric(Complete_TE_simple[row_RNA, ]); names(data.vec_RNA) = colnames(Complete_TE_simple)

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
row_LC = which(Complete_TE_simple$class=="Low_complexity")
data.vec_LC = as.numeric(Complete_TE_simple[row_LC, ]); names(data.vec_LC) = colnames(Complete_TE_simple)

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
row_MITE = which(Complete_TE_simple$class=="MITE")
data.vec_MITE = as.numeric(Complete_TE_simple[row_MITE, ]); names(data.vec_MITE) = colnames(Complete_TE_simple)

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
row_SAT = which(Complete_TE_simple$class=="Satellite")
data.vec_SAT = as.numeric(Complete_TE_simple[row_SAT, ]); names(data.vec_SAT) = colnames(Complete_TE_simple)

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
row_SR = which(Complete_TE_simple$class=="Simple_repeat")
data.vec_SR = as.numeric(Complete_TE_simple[row_SR, ]); names(data.vec_SR) = colnames(Complete_TE_simple)

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