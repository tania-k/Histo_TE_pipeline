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

#setwd("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships")

G184AR <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapG184AR.out",col_names = F)
colnames(G184AR) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
G184AR <- G184AR %>% mutate(length=endq-beginq +1)
write_tsv(G184AR,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/G184AR_adjust.tsv")
G184AR <- dplyr::mutate(G184AR, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
G184ARclass <- G184AR %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/30991520)
write_tsv(G184ARclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/G184ARclass_Starship.tsv")
G184ARclass <- G184ARclass %>% mutate('Histoplasma capsulatum G184AR' = percentage)
G184ARclass1 <- G184ARclass %>% select(-c(n,total,percentage))
HiscapG184AR <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_G184ARclass.tsv")
total <- HiscapG184AR %>% left_join(G184ARclass, by = "class")
df <- as.data.frame(total) 
df[is.na(df)] <- 0
df_new <- df %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_G184AR <- df_new %>% select(c(class,n,total,percentage))
write_tsv(df_G184AR, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/Attempt2/Attempt3/jbrowse_visualize/G184ARclass_TOTAL.tsv")
df_G184AR1 <- df_G184AR %>% mutate('Histoplasma capsulatum G184AR' = percentage)
df_G184AR11 <- df_G184AR1 %>% select(-c(n,total,percentage))
##


G186AR <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapG186AR.out",col_names = F)
colnames(G186AR) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
G186AR <- G186AR %>% mutate(length=endq-beginq +1)
write_tsv(G186AR,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/G186AR_adjust.tsv")
G186AR <- dplyr::mutate(G186AR, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
G186ARclass <- G186AR %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/31111494)
write_tsv(G186ARclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/G186ARclass_Starship.tsv")
G186ARclass <- G186ARclass %>% mutate('Histoplasma capsulatum G186AR' = percentage)
G186ARclass1 <- G186ARclass %>% select(-c(n,total,percentage))
HiscapG186AR <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_G186ARclass.tsv")
totalG186AR <- HiscapG186AR %>% left_join(G186ARclass, by = "class")
dfG186AR <- as.data.frame(totalG186AR) 
dfG186AR[is.na(dfG186AR)] <- 0
dfG186AR_new <- dfG186AR %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_G186AR <- dfG186AR_new %>% select(c(class,n,total,percentage))
write_tsv(df_G186AR, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/G186ARclass_TOTAL.tsv")
df_G186AR1 <- df_G186AR %>% mutate('Histoplasma capsulatum G186AR' = percentage)
df_G186AR11 <- df_G186AR1 %>% select(-c(n,total,percentage))

###

G217B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiG217B.out",col_names = F)
colnames(G217B) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
G217B <- G217B %>% mutate(length=endq-beginq +1)
write_tsv(G217B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/G217B_adjust.tsv")
G217B <- dplyr::mutate(G217B, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
G217Bclass <- G217B %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/39447273)
write_tsv(G217Bclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/G217Bclass_Starship.tsv")
G217Bclass <- G217Bclass %>% mutate('Histoplasma ohiense G217B' = percentage)
G217Bclass1 <- G217Bclass %>% select(-c(n,total,percentage))
HiscapG217B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_G217Bclass.tsv")
totalG217B <- HiscapG217B %>% left_join(G217Bclass, by = "class")
dfG217B <- as.data.frame(totalG217B) 
dfG217B[is.na(dfG217B)] <- 0
dfG217B_new <- dfG217B %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_G217B <- dfG217B_new %>% select(c(class,n,total,percentage))
write_tsv(df_G217B, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/G217Bclass_TOTAL.tsv")
df_G217B1 <- df_G186AR %>% mutate('Histoplasma capsulatum G217B' = percentage)
df_G217B11 <- df_G186AR1 %>% select(-c(n,total,percentage))
##

WU24 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissWU24.out",col_names = F)
colnames(WU24) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
WU24 <- WU24 %>% mutate(length=endq-beginq +1)
write_tsv(WU24,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/WU24_adjust.tsv")
WU24 <- dplyr::mutate(WU24, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
WU24class <- WU24 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/32531515)
write_tsv(WU24class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/WU24class_Starship.tsv")
WU24class <- WU24class %>% mutate('Histoplasma mississippiense WU24' = percentage)
WU24class1 <- WU24class %>% select(-c(n,total,percentage))
HistomissWU24 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_WU24class.tsv")
totalWU24 <- HistomissWU24 %>% left_join(WU24class, by = "class")
dfWU24 <- as.data.frame(totalWU24) 
dfWU24[is.na(dfWU24)] <- 0
dfWU24_new <- dfWU24 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_WU24 <- dfWU24_new %>% select(c(class,n,total,percentage))
write_tsv(df_WU24, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/WU24class_TOTAL.tsv")
df_WU241 <- df_WU24 %>% mutate('Histoplasma mississippiense WU24' = percentage)
df_WU2411 <- df_WU241 %>% select(-c(n,total,percentage))
##

H88 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapH88.out",col_names = F)
colnames(H88) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H88 <- H88 %>% mutate(length=endq-beginq +1)
write_tsv(H88,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H88_adjust.tsv")
H88 <- dplyr::mutate(H88, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H88class <- H88 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/37996987)
write_tsv(H88class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H88class_Starship.tsv")
H88class <- H88class %>% mutate('Histoplasma capsulatum H88' = percentage)
H88class1 <- H88class %>% select(-c(n,total,percentage))
HistocapH88 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_var.duboisii_H88class.tsv")
totalH88 <- HistocapH88 %>% left_join(H88class, by = "class")
dfH88 <- as.data.frame(totalH88) 
dfH88[is.na(dfH88)] <- 0
dfH88_new <- dfH88 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H88 <- dfH88_new %>% select(c(class,n,total,percentage))
write_tsv(df_H88, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H88class_TOTAL.tsv")
df_H881 <- df_H88 %>% mutate('Histoplasma capsulatum H88' = percentage)
df_H8811 <- df_H881 %>% select(-c(n,total,percentage))
##

Hcap19VMG15 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap19VMG15.out",col_names = F)
colnames(Hcap19VMG15) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
Hcap19VMG15 <- Hcap19VMG15 %>% mutate(length=endq-beginq +1)
write_tsv(Hcap19VMG15,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Hcap19VMG15_adjust.tsv")
Hcap19VMG15 <- dplyr::mutate(Hcap19VMG15, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
Hcap19VMG15class <- Hcap19VMG15 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26953388)
write_tsv(Hcap19VMG15class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Hcap19VMG15class_Starship.tsv")
Hcap19VMG15class <- Hcap19VMG15class %>% mutate('Histoplasma capsulatum 19VMG-15' = percentage)
Hcap19VMG15class1 <- Hcap19VMG15class %>% select(-c(n,total,percentage))
HHcap19VMG15 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_19VMG15class.tsv")
totalHcap19VMG15  <- HHcap19VMG15  %>% left_join(Hcap19VMG15class, by = "class")
dfHcap19VMG15 <- as.data.frame(totalHcap19VMG15) 
dfHcap19VMG15[is.na(dfHcap19VMG15)] <- 0
dfHcap19VMG15_new <- dfHcap19VMG15 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_Hcap19VMG15 <- dfHcap19VMG15_new %>% select(c(class,n,total,percentage))
write_tsv(df_Hcap19VMG15, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Hcap19VMG15class_TOTAL.tsv")
df_Hcap19VMG151 <- df_Hcap19VMG15 %>% mutate('Histoplasma capsulatum 19VMG-15' = percentage)
df_Hcap19VMG1511 <- df_Hcap19VMG151 %>% select(-c(n,total,percentage))
##

HISSPCM7256 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap19VMG15.out",col_names = F)
colnames(HISSPCM7256) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPCM7256 <- HISSPCM7256 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPCM7256,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPCM7256_adjust.tsv")
HISSPCM7256 <- dplyr::mutate(HISSPCM7256, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPCM7256class <- HISSPCM7256 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/59797550)
write_tsv(HISSPCM7256class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPCM7256class_Starship.tsv")
HISSPCM7256class <- HISSPCM7256class %>% mutate('Histoplasma capsulatum HISSP-CM7256' = percentage)
HISSPCM7256class1 <- HISSPCM7256class %>% select(-c(n,total,percentage))
HcapHISSPCM7256<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM7256-xx-CL-SEN-xxxx-036-BBclass.tsv")
totalHISSPCM7256  <- HcapHISSPCM7256  %>% left_join(HISSPCM7256class, by = "class")
dfHISSPCM7256 <- as.data.frame(totalHISSPCM7256) 
dfHISSPCM7256[is.na(dfHISSPCM7256)] <- 0
dfHISSPCM7256_new <- dfHISSPCM7256 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPCM7256 <- dfHISSPCM7256_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPCM7256, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPCM7256class_TOTAL.tsv")
df_HISSPCM72561 <- df_HISSPCM7256 %>% mutate('Histoplasma capsulatum 19VMG-15' = percentage)
df_HISSPCM725611 <- df_HISSPCM72561 %>% select(-c(n,total,percentage))
##

Histo485P20 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap485P20.out",col_names = F)
colnames(Histo485P20) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
Histo485P20 <- Histo485P20 %>% mutate(length=endq-beginq +1)
write_tsv(Histo485P20,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPCM7256_adjust.tsv")
Histo485P20 <- dplyr::mutate(Histo485P20, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
Histo485P20class <- Histo485P20 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28091669)
write_tsv(Histo485P20class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histo485P20class_Starship.tsv")
Histo485P20class <- Histo485P20class %>% mutate('Histoplasma capsulatum Histo485P20' = percentage)
Histo485P20class1 <- Histo485P20class %>% select(-c(n,total,percentage))
HcapHisto485P20<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Histo-485P20class.tsv")
totalHisto485P20  <- HcapHisto485P20  %>% left_join(Histo485P20class, by = "class")
dfHisto485P20 <- as.data.frame(totalHisto485P20) 
dfHisto485P20[is.na(dfHisto485P20)] <- 0
dfHisto485P20_new <- dfHisto485P20 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_Histo485P20 <- dfHisto485P20_new %>% select(c(class,n,total,percentage))
write_tsv(df_Histo485P20, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histo485P20class_TOTAL.tsv")
df_Histo485P201 <- df_Histo485P20 %>% mutate('Histoplasma capsulatum Histo485P20' = percentage)
df_Histo485P2011 <- df_Histo485P201 %>% select(-c(n,total,percentage))
##

JB_083285 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapJBHc083285.out",col_names = F)
colnames(JB_083285) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
JB_083285 <- JB_083285 %>% mutate(length=endq-beginq +1)
write_tsv(JB_083285,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/JB_083285_adjust.tsv")
JB_083285 <- dplyr::mutate(JB_083285, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
JB_083285class <- JB_083285 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/27324225)
write_tsv(JB_083285class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/JB_083285class_Starship.tsv")
JB_083285class <- JB_083285class %>% mutate('Histoplasma capsulatum JB_083285_2-Hc_083285_2' = percentage)
JB_083285class1 <- JB_083285class %>% select(-c(n,total,percentage))
HcapJB_083285<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_JB_083285_2-Hc_083285_2class.tsv")
totalJB_083285  <- HcapJB_083285  %>% left_join(JB_083285class, by = "class")
dfJB_083285 <- as.data.frame(totalJB_083285) 
dfJB_083285[is.na(dfJB_083285)] <- 0
dfJB_083285_new <- dfJB_083285 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_JB_083285 <- dfJB_083285_new %>% select(c(class,n,total,percentage))
write_tsv(df_JB_083285, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/JB_083285class_TOTAL.tsv")
df_JB_0832851 <- df_JB_083285 %>% mutate('Histoplasma capsulatum JB_083285_2-Hc_083285_2' = percentage)
df_JB_08328511 <- df_JB_0832851 %>% select(-c(n,total,percentage))
##

SECH_100 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapSECH100.out",col_names = F)
colnames(SECH_100) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_100 <- SECH_100 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_100,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_100_adjust.tsv")
SECH_100 <- dplyr::mutate(SECH_100, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_100class <- SECH_100 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/15609890)
write_tsv(SECH_100class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_100class_Starship.tsv")
SECH_100class <- SECH_100class %>% mutate('Histoplasma capsulatum SECH_100-Nam2_G186A' = percentage)
SECH_100class1 <- SECH_100class %>% select(-c(n,total,percentage))
HcapSECH_100<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_100-Nam2_G186Aclass.tsv")
totalSECH_100  <- HcapSECH_100  %>% left_join(SECH_100class, by = "class")
dfSECH_100 <- as.data.frame(totalSECH_100) 
dfSECH_100[is.na(dfSECH_100)] <- 0
dfSECH_100_new <- dfSECH_100 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_100 <- dfSECH_100_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_100, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_100class_TOTAL.tsv")
df_SECH_1001 <- df_SECH_100 %>% mutate('Histoplasma capsulatum SECH_100-G186A' = percentage)
df_SECH_10011 <- df_SECH_1001 %>% select(-c(n,total,percentage))
##

SECH_101 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapSECH101.out",col_names = F)
colnames(SECH_101) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_101 <- SECH_101 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_101,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_101_adjust.tsv")
SECH_101 <- dplyr::mutate(SECH_101, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_101class <- SECH_101 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23224387)
write_tsv(SECH_101class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_101class_Starship.tsv")
SECH_101class <- SECH_101class %>% mutate('Histoplasma capsulatum SECH_101-Nam2_G184A' = percentage)
SECH_101class1 <- SECH_101class %>% select(-c(n,total,percentage))
HcapSECH_101<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_101-Nam2_G184Aclass.tsv")
totalSECH_101  <- HcapSECH_101  %>% left_join(SECH_101class, by = "class")
dfSECH_101 <- as.data.frame(totalSECH_101) 
dfSECH_101[is.na(dfSECH_101)] <- 0
dfSECH_101_new <- dfSECH_101 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_101 <- dfSECH_101_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_101, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_101class_TOTAL.tsv")
df_SECH_1011 <- df_SECH_101 %>% mutate('Histoplasma capsulatum SECH_101-G184A' = percentage)
df_SECH_10111 <- df_SECH_1011 %>% select(-c(n,total,percentage))
##

SECH_102 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH102.out",col_names = F)
colnames(SECH_102) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_102 <- SECH_102 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_102,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_102_adjust.tsv")
SECH_102 <- dplyr::mutate(SECH_102, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_102class <- SECH_102 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/27544444)
write_tsv(SECH_102class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_102class_Starship.tsv")
SECH_102class <- SECH_102class %>% mutate('Histoplasma mississippiense SECH_102 Nam2_505' = percentage)
SECH_102class1 <- SECH_102class %>% select(-c(n,total,percentage))
HcapSECH_102<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_102Nam2_505class.tsv")
totalSECH_102  <- HcapSECH_102  %>% left_join(SECH_102class, by = "class")
dfSECH_102 <- as.data.frame(totalSECH_102) 
dfSECH_102[is.na(dfSECH_102)] <- 0
dfSECH_102_new <- dfSECH_102 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_102 <- dfSECH_102_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_102, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_102class_TOTAL.tsv")
df_SECH_1021 <- df_SECH_102 %>% mutate('Histoplasma capsulatum SECH_102-505' = percentage)
df_SECH_10211 <- df_SECH_1021 %>% select(-c(n,total,percentage))
##

SECH_103 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistosurSECH103.out",col_names = F)
colnames(SECH_103) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_103 <- SECH_103 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_103,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_103_adjust.tsv")
SECH_103 <- dplyr::mutate(SECH_103, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_103class <- SECH_103 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22453624)
write_tsv(SECH_103class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_103class_Starship.tsv")
SECH_103class <- SECH_103class %>% mutate('Histoplasma capsulatum SECH_103 Nam2_3_11G' = percentage)
SECH_103class1 <- SECH_103class %>% select(-c(n,total,percentage))
HcapSECH_103<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_103-Nam2_3_11Gclass.tsv")
totalSECH_103  <- HcapSECH_103  %>% left_join(SECH_103class, by = "class")
dfSECH_103 <- as.data.frame(totalSECH_103) 
dfSECH_103[is.na(dfSECH_103)] <- 0
dfSECH_103_new <- dfSECH_103 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_103 <- dfSECH_103_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_103, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_103class_TOTAL.tsv")
df_SECH_1031 <- df_SECH_103 %>% mutate('Histoplasma capsulatum SECH_103-3_11G' = percentage)
df_SECH_10311 <- df_SECH_1031 %>% select(-c(n,total,percentage))
##

SECH_104 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistosurSECH104.out",col_names = F)
colnames(SECH_104) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_104 <- SECH_104 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_104,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_104_adjust.tsv")
SECH_104 <- dplyr::mutate(SECH_104, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_104class <- SECH_104 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22213937)
write_tsv(SECH_104class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_104class_Starship.tsv")
SECH_104class <- SECH_104class %>% mutate('Histoplasma capsulatum SECH_104 Nam2_27_14' = percentage)
SECH_104class1 <- SECH_104class %>% select(-c(n,total,percentage))
HcapSECH_104<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_104-Nam2_27_14class.tsv")
totalSECH_104  <- HcapSECH_104 %>% left_join(SECH_104class, by = "class")
dfSECH_104 <- as.data.frame(totalSECH_104) 
dfSECH_104[is.na(dfSECH_104)] <- 0
dfSECH_104_new <- dfSECH_104 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_104 <- dfSECH_104_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_104, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_104class_TOTAL.tsv")
df_SECH_1041 <- df_SECH_104 %>% mutate('Histoplasma capsulatum SECH_104-27_14' = percentage)
df_SECH_10411 <- df_SECH_1041 %>% select(-c(n,total,percentage))
##

SECH_105 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistosurSECH105.out",col_names = F)
colnames(SECH_105) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_105 <- SECH_105 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_105,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_105_adjust.tsv")
SECH_105 <- dplyr::mutate(SECH_105, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_105class <- SECH_105 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28048040)
write_tsv(SECH_105class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_105class_Starship.tsv")
SECH_105class <- SECH_105class %>% mutate('Histoplasma capsulatum SECH_105 Nam2_21_14' = percentage)
SECH_105class1 <- SECH_105class %>% select(-c(n,total,percentage))
HcapSECH_105<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_sum/Histoplasma_suramericanum_SECH_105-Nam2_21_14class.tsv")
totalSECH_105  <- HcapSECH_105 %>% left_join(SECH_105class, by = "class")
dfSECH_105 <- as.data.frame(totalSECH_105) 
dfSECH_105[is.na(dfSECH_105)] <- 0
dfSECH_105_new <- dfSECH_105 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_105 <- dfSECH_105_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_105, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_105class_TOTAL.tsv")
df_SECH_1051 <- df_SECH_105 %>% mutate('Histoplasma capsulatum SECH_105-21_14' = percentage)
df_SECH_10511 <- df_SECH_1051 %>% select(-c(n,total,percentage))
##

SECH_107 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapSECH107.out",col_names = F)
colnames(SECH_107) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_107 <- SECH_107 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_107,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_107_adjust.tsv")
SECH_107 <- dplyr::mutate(SECH_107, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_107class <- SECH_107 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22891084)
write_tsv(SECH_107class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_107class_Starship.tsv")
SECH_107class <- SECH_107class %>% mutate('Histoplasma capsulatum SECH_107-mis_Hc_duboisii-B' = percentage)
SECH_107class1 <- SECH_107class %>% select(-c(n,total,percentage))
HcapSECH_107<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_107-mis_Hc_duboisii-Bclass.tsv")
totalSECH_107  <- HcapSECH_107 %>% left_join(SECH_107class, by = "class")
dfSECH_107 <- as.data.frame(totalSECH_107) 
dfSECH_107[is.na(dfSECH_107)] <- 0
dfSECH_107_new <- dfSECH_107 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_107 <- dfSECH_107_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_107, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_107class_TOTAL.tsv")
df_SECH_1071 <- df_SECH_107 %>% mutate('Histoplasma capsulatum SECH_105-Hc_duboisii-B' = percentage)
df_SECH_10711 <- df_SECH_1071 %>% select(-c(n,total,percentage))
##

SECH_110 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapSECH110.out",col_names = F)
colnames(SECH_110) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_110 <- SECH_110 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_110,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_110_adjust.tsv")
SECH_110 <- dplyr::mutate(SECH_110, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_110class <- SECH_110 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22122765)
write_tsv(SECH_110class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_110class_Starship.tsv")
SECH_110class <- SECH_110class %>% mutate('Histoplasma capsulatum SECH_110' = percentage)
SECH_110class1 <- SECH_110class %>% select(-c(n,total,percentage))
HcapSECH_110<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SECH_110class.tsv")
totalSECH_110  <- HcapSECH_110 %>% left_join(SECH_110class, by = "class")
dfSECH_110 <- as.data.frame(totalSECH_110) 
dfSECH_110[is.na(dfSECH_110)] <- 0
dfSECH_110_new <- dfSECH_110 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_110 <- dfSECH_110_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_110, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_110class_TOTAL.tsv")
df_SECH_1101 <- df_SECH_110 %>% mutate('Histoplasma capsulatum SECH_110' = percentage)
df_SECH_11011 <- df_SECH_1101 %>% select(-c(n,total,percentage))
##

SECH_81 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH81.out",col_names = F)
colnames(SECH_81) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_81 <- SECH_81 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_81,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_81_adjust.tsv")
SECH_81 <- dplyr::mutate(SECH_81, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_81class <- SECH_81 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25518040)
write_tsv(SECH_81class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_81class_Starship.tsv")
SECH_81class <- SECH_81class %>% mutate('Histoplasma mississippiense SECH_81 WU24' = percentage)
SECH_81class1 <- SECH_81class %>% select(-c(n,total,percentage))
HcapSECH_81<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_81-Nam1_WU24class.tsv")
totalSECH_81  <- HcapSECH_81 %>% left_join(SECH_81class, by = "class")
dfSECH_81 <- as.data.frame(totalSECH_81) 
dfSECH_81[is.na(dfSECH_81)] <- 0
dfSECH_81_new <- dfSECH_81 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_81 <- dfSECH_81_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_81, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_81class_TOTAL.tsv")
df_SECH_811 <- df_SECH_81 %>% mutate('Histoplasma mississippiense SECH_81-WU24' = percentage)
df_SECH_8111 <- df_SECH_811 %>% select(-c(n,total,percentage))
##

SECH_82 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH82.out",col_names = F)
colnames(SECH_82) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_82 <- SECH_82 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_82,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_82_adjust.tsv")
SECH_82 <- dplyr::mutate(SECH_82, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_82class <- SECH_82 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/30880349)
write_tsv(SECH_82class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_82class_Starship.tsv")
SECH_82class <- SECH_82class %>% mutate('Histoplasma ohiense SECH_82 CI_4' = percentage)
SECH_82class1 <- SECH_82class %>% select(-c(n,total,percentage))
HcapSECH_82<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_82-Nam1_CI_4class.tsv")
totalSECH_82  <- HcapSECH_82 %>% left_join(SECH_82class, by = "class")
dfSECH_82 <- as.data.frame(totalSECH_82) 
dfSECH_82[is.na(dfSECH_82)] <- 0
dfSECH_82_new <- dfSECH_82 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_82 <- dfSECH_82_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_82, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_82class_TOTAL.tsv")
df_SECH_821 <- df_SECH_82 %>% mutate('Histoplasma ohiense SECH_82-CI_4' = percentage)
df_SECH_8211 <- df_SECH_821 %>% select(-c(n,total,percentage))
##

SECH_83 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH83.out",col_names = F)
colnames(SECH_83) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_83 <- SECH_83 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_83,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_83_adjust.tsv")
SECH_83 <- dplyr::mutate(SECH_83, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_83class <- SECH_83 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25776464)
write_tsv(SECH_83class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_83class_Starship.tsv")
SECH_83class <- SECH_83class %>% mutate('Histoplasma mississippiense SECH_83 CI_7' = percentage)
SECH_83class1 <- SECH_83class %>% select(-c(n,total,percentage))
HcapSECH_83<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_83-Nam1_CI_7class.tsv")
totalSECH_83  <- HcapSECH_83 %>% left_join(SECH_83class, by = "class")
dfSECH_83 <- as.data.frame(totalSECH_83) 
dfSECH_83[is.na(dfSECH_83)] <- 0
dfSECH_83_new <- dfSECH_83 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_83 <- dfSECH_83_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_83, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_83class_TOTAL.tsv")
df_SECH_831 <- df_SECH_83 %>% mutate('Histoplasma mississippiense SECH_83-CI_7' = percentage)
df_SECH_8311 <- df_SECH_831 %>% select(-c(n,total,percentage))
##

SECH_84 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH84.out",col_names = F)
colnames(SECH_84) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_84 <- SECH_84 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_84,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_84_adjust.tsv")
SECH_84 <- dplyr::mutate(SECH_84, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_84class <- SECH_84 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26524728)
write_tsv(SECH_84class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_84class_Starship.tsv")
SECH_84class <- SECH_84class %>% mutate('Histoplasma mississippiense SECH_84 CI_19' = percentage)
SECH_84class1 <- SECH_84class %>% select(-c(n,total,percentage))
HcapSECH_84<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_84-Nam1_CI_19class.tsv")
totalSECH_84  <- HcapSECH_84 %>% left_join(SECH_84class, by = "class")
dfSECH_84 <- as.data.frame(totalSECH_84) 
dfSECH_84[is.na(dfSECH_84)] <- 0
dfSECH_84_new <- dfSECH_84 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_84 <- dfSECH_84_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_84, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_84class_TOTAL.tsv")
df_SECH_841 <- df_SECH_84 %>% mutate('Histoplasma mississippiense SECH_84-CI_19' = percentage)
df_SECH_8411 <- df_SECH_841 %>% select(-c(n,total,percentage))
##

SECH_85 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH85.out",col_names = F)
colnames(SECH_85) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_85 <- SECH_85 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_85,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_85_adjust.tsv")
SECH_85 <- dplyr::mutate(SECH_85, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_85class <- SECH_85 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23145789)
write_tsv(SECH_85class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_85class_Starship.tsv")
SECH_85class <- SECH_85class %>% mutate('Histoplasma mississippiense SECH_85 CI_22' = percentage)
SECH_85class1 <- SECH_85class %>% select(-c(n,total,percentage))
HcapSECH_85<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_85-Nam1_CI_22class.tsv")
totalSECH_85  <- HcapSECH_85 %>% left_join(SECH_85class, by = "class")
dfSECH_85 <- as.data.frame(totalSECH_85) 
dfSECH_85[is.na(dfSECH_85)] <- 0
dfSECH_85_new <- dfSECH_85 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_85 <- dfSECH_85_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_85, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_85class_TOTAL.tsv")
df_SECH_851 <- df_SECH_85 %>% mutate('Histoplasma mississippiense SECH_85-CI_22' = percentage)
df_SECH_8511 <- df_SECH_851 %>% select(-c(n,total,percentage))
##

SECH_86 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH86.out",col_names = F)
colnames(SECH_86) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_86 <- SECH_86 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_86,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_86_adjust.tsv")
SECH_86 <- dplyr::mutate(SECH_86, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_86class <- SECH_86 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23809594)
write_tsv(SECH_86class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_86class_Starship.tsv")
SECH_86class <- SECH_86class %>% mutate('Histoplasma mississippiense SECH_86 CI_24' = percentage)
SECH_86class1 <- SECH_86class %>% select(-c(n,total,percentage))
HcapSECH_86<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_86-Nam1_CI_24class.tsv")
totalSECH_86  <- HcapSECH_86 %>% left_join(SECH_86class, by = "class")
dfSECH_86 <- as.data.frame(totalSECH_86) 
dfSECH_86[is.na(dfSECH_86)] <- 0
dfSECH_86_new <- dfSECH_86 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_86 <- dfSECH_86_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_86, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_86class_TOTAL.tsv")
df_SECH_861 <- df_SECH_86 %>% mutate('Histoplasma mississippiense SECH_86-CI_24' = percentage)
df_SECH_8611 <- df_SECH_861 %>% select(-c(n,total,percentage))
##

SECH_87 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH87.out",col_names = F)
colnames(SECH_87) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_87 <- SECH_87 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_87,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_87_adjust.tsv")
SECH_87 <- dplyr::mutate(SECH_87, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_87class <- SECH_87 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26345314)
write_tsv(SECH_87class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_87class_Starship.tsv")
SECH_87class <- SECH_87class %>% mutate('Histoplasma mississippiense SECH_87 CI_42' = percentage)
SECH_87class1 <- SECH_87class %>% select(-c(n,total,percentage))
HcapSECH_87<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_87-Nam1_CI_42class.tsv")
totalSECH_87 <- HcapSECH_87 %>% left_join(SECH_87class, by = "class")
dfSECH_87 <- as.data.frame(totalSECH_87) 
dfSECH_87[is.na(dfSECH_87)] <- 0
dfSECH_87_new <- dfSECH_87 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_87 <- dfSECH_87_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_87, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_87class_TOTAL.tsv")
df_SECH_871 <- df_SECH_87 %>% mutate('Histoplasma mississippiense SECH_87-CI_42' = percentage)
df_SECH_8711 <- df_SECH_871 %>% select(-c(n,total,percentage))
##

SECH_88 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH88.out",col_names = F)
colnames(SECH_88) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_88 <- SECH_88 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_88,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_88_adjust.tsv")
SECH_88 <- dplyr::mutate(SECH_88, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_88class <- SECH_88 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24925558)
write_tsv(SECH_88class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_88class_Starship.tsv")
SECH_88class <- SECH_88class %>% mutate('Histoplasma mississippiense SECH_88 CI_43' = percentage)
SECH_88class1 <- SECH_88class %>% select(-c(n,total,percentage))
HcapSECH_88<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_88-Nam1_CI_43class.tsv")
totalSECH_88 <- HcapSECH_88 %>% left_join(SECH_88class, by = "class")
dfSECH_88 <- as.data.frame(totalSECH_88) 
dfSECH_88[is.na(dfSECH_88)] <- 0
dfSECH_88_new <- dfSECH_88 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_88 <- dfSECH_88_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_88, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_88class_TOTAL.tsv")
df_SECH_881 <- df_SECH_88 %>% mutate('Histoplasma mississippiense SECH_88-CI_43' = percentage)
df_SECH_8811 <- df_SECH_881 %>% select(-c(n,total,percentage))
##

SECH_89 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH89.out",col_names = F)
colnames(SECH_89) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_89 <- SECH_89 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_89,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_89_adjust.tsv")
SECH_89 <- dplyr::mutate(SECH_89, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_89class <- SECH_89 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24430709)
write_tsv(SECH_89class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_89class_Starship.tsv")
SECH_89class <- SECH_89class %>% mutate('Histoplasma mississippiense SECH_89 UCLA531' = percentage)
SECH_89class1 <- SECH_89class %>% select(-c(n,total,percentage))
HcapSECH_89<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_89-Nam1_UCLA-531class.tsv")
totalSECH_89 <- HcapSECH_89 %>% left_join(SECH_89class, by = "class")
dfSECH_89 <- as.data.frame(totalSECH_89) 
dfSECH_89[is.na(dfSECH_89)] <- 0
dfSECH_89_new <- dfSECH_89 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_89 <- dfSECH_89_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_89, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_89class_TOTAL.tsv")
df_SECH_891 <- df_SECH_89 %>% mutate('Histoplasma mississippiense SECH_89-UCLA-531' = percentage)
df_SECH_8911 <- df_SECH_891 %>% select(-c(n,total,percentage))
##

SECH_90 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistomissSECH90.out",col_names = F)
colnames(SECH_90) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_90 <- SECH_90 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_90,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_90_adjust.tsv")
SECH_90 <- dplyr::mutate(SECH_90, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_90class <- SECH_90 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24054288)
write_tsv(SECH_90class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_90class_Starship.tsv")
SECH_90class <- SECH_90class %>% mutate('Histoplasma mississippiense SECH_90 DOWNS' = percentage)
SECH_90class1 <- SECH_90class %>% select(-c(n,total,percentage))
HcapSECH_90<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_mis/Histoplasma_mississippienseII_SECH_90-Nam1_DOWNSclass.tsv")
totalSECH_90 <- HcapSECH_90 %>% left_join(SECH_90class, by = "class")
dfSECH_90 <- as.data.frame(totalSECH_90) 
dfSECH_90[is.na(dfSECH_90)] <- 0
dfSECH_90_new <- dfSECH_90 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_90 <- dfSECH_90_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_90, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_90class_TOTAL.tsv")
df_SECH_901 <- df_SECH_90 %>% mutate('Histoplasma mississippiense SECH_90-DOWNS' = percentage)
df_SECH_9011 <- df_SECH_901 %>% select(-c(n,total,percentage))
##

SECH_91 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH91.out",col_names = F)
colnames(SECH_91) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_91 <- SECH_91 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_91,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_91_adjust.tsv")
SECH_91 <- dplyr::mutate(SECH_91, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_91class <- SECH_91 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29816969)
write_tsv(SECH_91class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_91class_Starship.tsv")
SECH_91class <- SECH_91class %>% mutate('Histoplasma ohiense SECH_91 G217B' = percentage)
SECH_91class1 <- SECH_91class %>% select(-c(n,total,percentage))
HcapSECH_91<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_91-Nam2_G217Bclass.tsv")
totalSECH_91 <- HcapSECH_91 %>% left_join(SECH_91class, by = "class")
dfSECH_91 <- as.data.frame(totalSECH_91) 
dfSECH_91[is.na(dfSECH_91)] <- 0
dfSECH_91_new <- dfSECH_91 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_91 <- dfSECH_91_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_91, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_91class_TOTAL.tsv")
df_SECH_911 <- df_SECH_91 %>% mutate('Histoplasma ohiense SECH_91-G217B' = percentage)
df_SECH_9111 <- df_SECH_911 %>% select(-c(n,total,percentage))
##

SECH_92 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH92.out",col_names = F)
colnames(SECH_92) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_92 <- SECH_92 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_92,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_92_adjust.tsv")
SECH_92 <- dplyr::mutate(SECH_92, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_92class <- SECH_92 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29639711)
write_tsv(SECH_92class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_92class_Starship.tsv")
SECH_92class <- SECH_92class %>% mutate('Histoplasma ohiense SECH_92 G222B' = percentage)
SECH_92class1 <- SECH_92class %>% select(-c(n,total,percentage))
HcapSECH_92<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_92-Nam2_G222Bclass.tsv")
totalSECH_92 <- HcapSECH_92 %>% left_join(SECH_92class, by = "class")
dfSECH_92 <- as.data.frame(totalSECH_92) 
dfSECH_92[is.na(dfSECH_92)] <- 0
dfSECH_92_new <- dfSECH_92 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_92 <- dfSECH_92_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_92, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_92class_TOTAL.tsv")
df_SECH_921 <- df_SECH_92 %>% mutate('Histoplasma ohiense SECH_92-G222B' = percentage)
df_SECH_9211 <- df_SECH_921 %>% select(-c(n,total,percentage))
##

SECH_93 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH93.out",col_names = F)
colnames(SECH_93) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_93 <- SECH_93 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_93,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_93_adjust.tsv")
SECH_93 <- dplyr::mutate(SECH_93, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_93class <- SECH_93 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23646057)
write_tsv(SECH_93class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_93class_Starship.tsv")
SECH_93class <- SECH_93class %>% mutate('Histoplasma ohiense SECH_93 CI_6' = percentage)
SECH_93class1 <- SECH_93class %>% select(-c(n,total,percentage))
HcapSECH_93<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_93-Nam2_CI_6class.tsv")
totalSECH_93 <- HcapSECH_93 %>% left_join(SECH_93class, by = "class")
dfSECH_93 <- as.data.frame(totalSECH_93) 
dfSECH_93[is.na(dfSECH_93)] <- 0
dfSECH_93_new <- dfSECH_93 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_93 <- dfSECH_93_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_93, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_93class_TOTAL.tsv")
df_SECH_931 <- df_SECH_93 %>% mutate('Histoplasma ohiense SECH_93-CI_6' = percentage)
df_SECH_9311 <- df_SECH_931 %>% select(-c(n,total,percentage))
##

SECH_94 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH94.out",col_names = F)
colnames(SECH_94) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_94 <- SECH_94 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_94,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_94_adjust.tsv")
SECH_94 <- dplyr::mutate(SECH_94, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_94class <- SECH_94 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/20233006)
write_tsv(SECH_94class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_94class_Starship.tsv")
SECH_94class <- SECH_94class %>% mutate('Histoplasma ohiense SECH_94 CI_9' = percentage)
SECH_94class1 <- SECH_94class %>% select(-c(n,total,percentage))
HcapSECH_94<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_94.Nam2_CI_9class.tsv")
totalSECH_94 <- HcapSECH_94 %>% left_join(SECH_94class, by = "class")
dfSECH_94 <- as.data.frame(totalSECH_94) 
dfSECH_94[is.na(dfSECH_94)] <- 0
dfSECH_94_new <- dfSECH_94 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_94 <- dfSECH_94_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_94, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_94class_TOTAL.tsv")
df_SECH_941 <- df_SECH_94 %>% mutate('Histoplasma ohiense SECH_94-CI_9' = percentage)
df_SECH_9411 <- df_SECH_941 %>% select(-c(n,total,percentage))
##

SECH_95 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH95.out",col_names = F)
colnames(SECH_95) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_95 <- SECH_95 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_95,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_95_adjust.tsv")
SECH_95 <- dplyr::mutate(SECH_95, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_95class <- SECH_95 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25510614)
write_tsv(SECH_95class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_95class_Starship.tsv")
SECH_95class <- SECH_95class %>% mutate('Histoplasma ohiense SECH_95 CI_10' = percentage)
SECH_95class1 <- SECH_95class %>% select(-c(n,total,percentage))
HcapSECH_95<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_95.Nam2_CI_10class.tsv")
totalSECH_95 <- HcapSECH_95 %>% left_join(SECH_95class, by = "class")
dfSECH_95 <- as.data.frame(totalSECH_95) 
dfSECH_95[is.na(dfSECH_95)] <- 0
dfSECH_95_new <- dfSECH_95 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_95 <- dfSECH_95_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_95, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_95class_TOTAL.tsv")
df_SECH_951 <- df_SECH_95 %>% mutate('Histoplasma ohiense SECH_95-CI_10' = percentage)
df_SECH_9511 <- df_SECH_951 %>% select(-c(n,total,percentage))
##

SECH_96 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH96.out",col_names = F)
colnames(SECH_96) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_96 <- SECH_96 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_96,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_96_adjust.tsv")
SECH_96 <- dplyr::mutate(SECH_96, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_96class <- SECH_96 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25893578)
write_tsv(SECH_96class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_96class_Starship.tsv")
SECH_96class <- SECH_96class %>% mutate('Histoplasma ohiense SECH_96 CI_17' = percentage)
SECH_96class1 <- SECH_96class %>% select(-c(n,total,percentage))
HcapSECH_96<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_96.Nam2_CI_17class.tsv")
totalSECH_96 <- HcapSECH_96 %>% left_join(SECH_96class, by = "class")
dfSECH_96 <- as.data.frame(totalSECH_96) 
dfSECH_96[is.na(dfSECH_96)] <- 0
dfSECH_96_new <- dfSECH_96 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_96 <- dfSECH_96_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_96, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_96class_TOTAL.tsv")
df_SECH_961 <- df_SECH_96 %>% mutate('Histoplasma ohiense SECH_96-CI_17' = percentage)
df_SECH_9611 <- df_SECH_961 %>% select(-c(n,total,percentage))
##

SECH_97 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH97.out",col_names = F)
colnames(SECH_97) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_97 <- SECH_97 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_97,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_97_adjust.tsv")
SECH_97 <- dplyr::mutate(SECH_97, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_97class <- SECH_97 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28259788)
write_tsv(SECH_97class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_97class_Starship.tsv")
SECH_97class <- SECH_97class %>% mutate('Histoplasma ohiense SECH_97 CI_18' = percentage)
SECH_97class1 <- SECH_97class %>% select(-c(n,total,percentage))
HcapSECH_97<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_97-Nam2_CI_18class.tsv")
totalSECH_97 <- HcapSECH_97 %>% left_join(SECH_97class, by = "class")
dfSECH_97 <- as.data.frame(totalSECH_97) 
dfSECH_97[is.na(dfSECH_97)] <- 0
dfSECH_97_new <- dfSECH_97 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_97 <- dfSECH_97_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_97, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_97class_TOTAL.tsv")
df_SECH_971 <- df_SECH_97 %>% mutate('Histoplasma ohiense SECH_97-CI_18' = percentage)
df_SECH_9711 <- df_SECH_971 %>% select(-c(n,total,percentage))
##

SECH_98 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH98.out",col_names = F)
colnames(SECH_98) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_98 <- SECH_98 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_98,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_98_adjust.tsv")
SECH_98 <- dplyr::mutate(SECH_98, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_98class <- SECH_98 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24259982)
write_tsv(SECH_98class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_98class_Starship.tsv")
SECH_98class <- SECH_98class %>% mutate('Histoplasma ohiense SECH_98 CI_30' = percentage)
SECH_98class1 <- SECH_98class %>% select(-c(n,total,percentage))
HcapSECH_98<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_98.Nam2_CI_30class.tsv")
totalSECH_98 <- HcapSECH_98 %>% left_join(SECH_98class, by = "class")
dfSECH_98 <- as.data.frame(totalSECH_98) 
dfSECH_98[is.na(dfSECH_98)] <- 0
dfSECH_98_new <- dfSECH_98 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_98 <- dfSECH_98_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_98, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_98class_TOTAL.tsv")
df_SECH_981 <- df_SECH_98 %>% mutate('Histoplasma ohiense SECH_98-CI_30' = percentage)
df_SECH_9811 <- df_SECH_981 %>% select(-c(n,total,percentage))
##

SECH_99 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistoohiSECH99.out",col_names = F)
colnames(SECH_99) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SECH_99 <- SECH_99 %>% mutate(length=endq-beginq +1)
write_tsv(SECH_99,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_99_adjust.tsv")
SECH_99 <- dplyr::mutate(SECH_99, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SECH_99class <- SECH_99 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29591009)
write_tsv(SECH_99class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_99class_Starship.tsv")
SECH_99class <- SECH_99class %>% mutate('Histoplasma ohiense SECH_99 CI_35' = percentage)
SECH_99class1 <- SECH_99class %>% select(-c(n,total,percentage))
HcapSECH_99<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_ohi/Histoplasma_ohiense_SECH_99.Nam2_CI_35class.tsv")
totalSECH_99 <- HcapSECH_99 %>% left_join(SECH_99class, by = "class")
dfSECH_99 <- as.data.frame(totalSECH_99) 
dfSECH_99[is.na(dfSECH_99)] <- 0
dfSECH_99_new <- dfSECH_99 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SECH_99 <- dfSECH_99_new %>% select(c(class,n,total,percentage))
write_tsv(df_SECH_99, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SECH_99class_TOTAL.tsv")
df_SECH_991 <- df_SECH_99 %>% mutate('Histoplasma ohiense SECH_99-CI_35' = percentage)
df_SECH_9911 <- df_SECH_991 %>% select(-c(n,total,percentage))
##

HCH143 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHCH143.out",col_names = F)
colnames(HCH143) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HCH143 <- HCH143 %>% mutate(length=endq-beginq +1)
write_tsv(HCH143,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HCH143_adjust.tsv")
HCH143 <- dplyr::mutate(HCH143, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HCH143class <- HCH143 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29042875)
write_tsv(HCH143class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HCH143class_Starship.tsv")
HCH143class <- HCH143class %>% mutate('Histoplasma capsulatum HCH143' = percentage)
HCH143class1 <- HCH143class %>% select(-c(n,total,percentage))
HcapHCH143<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HCH143class.tsv")
totalHCH143 <- HcapHCH143 %>% left_join(HCH143class, by = "class")
dfHCH143 <- as.data.frame(totalHCH143) 
dfHCH143[is.na(dfHCH143)] <- 0
dfHCH143_new <- dfHCH143 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HCH143 <- dfHCH143_new %>% select(c(class,n,total,percentage))
write_tsv(df_HCH143, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HCH143class_TOTAL.tsv")
df_HCH1431 <- df_HCH143 %>% mutate('Histoplasma capsulatum HCH143' = percentage)
df_HCH14311 <- df_HCH1431 %>% select(-c(n,total,percentage))
##

HISSPFGPIE2055 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPFGPIE2055.out",col_names = F)
colnames(HISSPFGPIE2055) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPFGPIE2055 <- HISSPFGPIE2055 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPFGPIE2055,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPIE2055_adjust.tsv")
HISSPFGPIE2055 <- dplyr::mutate(HISSPFGPIE2055, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPFGPIE2055class <- HISSPFGPIE2055 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28918902)
write_tsv(HISSPFGPIE2055class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPIE2055class_Starship.tsv")
HISSPFGPIE2055class <- HISSPFGPIE2055class %>% mutate('Histoplasma capsulatum HISSP-FGPIE2055' = percentage)
HISSPFGPIE2055class1 <- HISSPFGPIE2055class %>% select(-c(n,total,percentage))
HcapHISSPFGPIE2055<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPIE2055class.tsv")
totalHISSPFGPIE2055 <- HcapHISSPFGPIE2055 %>% left_join(HISSPFGPIE2055class, by = "class")
dfHISSPFGPIE2055 <- as.data.frame(totalHISSPFGPIE2055) 
dfHISSPFGPIE2055[is.na(dfHISSPFGPIE2055)] <- 0
dfHISSPFGPIE2055_new <- dfHISSPFGPIE2055 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPFGPIE2055 <- dfHISSPFGPIE2055_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPFGPIE2055, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPIE2055class_TOTAL.tsv")
df_HISSPFGPIE20551 <- df_HISSPFGPIE2055 %>% mutate('Histoplasma capsulatum HISSP-FGPIE2055' = percentage)
df_HISSPFGPIE205511 <- df_HISSPFGPIE20551 %>% select(-c(n,total,percentage))
##

HISSPFGPIA2052 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPFGPIA2052.out",col_names = F)
colnames(HISSPFGPIA2052) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPFGPIA2052 <- HISSPFGPIA2052 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPFGPIA2052,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPIA2052_adjust.tsv")
HISSPFGPIA2052 <- dplyr::mutate(HISSPFGPIA2052, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPFGPIA2052class <- HISSPFGPIA2052 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/20859097)
write_tsv(HISSPFGPIA2052class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPIA2052class_Starship.tsv")
HISSPFGPIA2052class <- HISSPFGPIA2052class %>% mutate('Histoplasma capsulatum HISSP-FGPIA2052' = percentage)
HISSPFGPIA2052class1 <- HISSPFGPIA2052class %>% select(-c(n,total,percentage))
HcapHISSPFGPIA2052<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPIA2052class.tsv")
totalHISSPFGPIA2052 <- HcapHISSPFGPIA2052 %>% left_join(HISSPFGPIA2052class, by = "class")
dfHISSPFGPIA2052 <- as.data.frame(totalHISSPFGPIA2052) 
dfHISSPFGPIA2052[is.na(dfHISSPFGPIA2052)] <- 0
dfHISSPFGPIA2052_new <- dfHISSPFGPIA2052 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPFGPIA2052 <- dfHISSPFGPIA2052_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPFGPIA2052, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPIA2052class_TOTAL.tsv")
df_HISSPFGPIA20521 <- df_HISSPFGPIA2052 %>% mutate('Histoplasma capsulatum HISSP-FGPIA2052' = percentage)
df_HISSPFGPIA205211 <- df_HISSPFGPIA20521 %>% select(-c(n,total,percentage))
##

HISSPFGPERS2034 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPFGPERS2034.out",col_names = F)
colnames(HISSPFGPERS2034) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPFGPERS2034 <- HISSPFGPERS2034 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPFGPERS2034,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPERS2034_adjust.tsv")
HISSPFGPERS2034 <- dplyr::mutate(HISSPFGPERS2034, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPFGPERS2034class <- HISSPFGPERS2034 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/24259670)
write_tsv(HISSPFGPERS2034class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPERS2034class_Starship.tsv")
HISSPFGPERS2034class <- HISSPFGPERS2034class %>% mutate('Histoplasma capsulatum HISSP-FGPERS2034' = percentage)
HISSPFGPERS2034class1 <- HISSPFGPERS2034class %>% select(-c(n,total,percentage))
HcapHISSPFGPERS2034<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGPERS2034class.tsv")
totalHISSPFGPERS2034 <- HcapHISSPFGPERS2034 %>% left_join(HISSPFGPERS2034class, by = "class")
dfHISSPFGPERS2034 <- as.data.frame(totalHISSPFGPERS2034) 
dfHISSPFGPERS2034[is.na(dfHISSPFGPERS2034)] <- 0
dfHISSPFGPERS2034_new <- dfHISSPFGPERS2034 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPFGPERS2034 <- dfHISSPFGPERS2034_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPFGPERS2034, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGPERS2034class_TOTAL.tsv")
df_HISSPFGPERS20341 <- df_HISSPFGPERS2034 %>% mutate('Histoplasma capsulatum HISSP-FGPERS2034' = percentage)
df_HISSPFGPERS203411 <- df_HISSPFGPERS20341 %>% select(-c(n,total,percentage))
##

HISSPFGFER2036 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPFGFER2036.out",col_names = F)
colnames(HISSPFGFER2036) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPFGFER2036 <- HISSPFGFER2036 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPFGFER2036,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGFER2036_adjust.tsv")
HISSPFGFER2036 <- dplyr::mutate(HISSPFGFER2036, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPFGFER2036class <- HISSPFGFER2036 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/22918657)
write_tsv(HISSPFGFER2036class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGFER2036class_Starship.tsv")
HISSPFGFER2036class <- HISSPFGFER2036class %>% mutate('Histoplasma capsulatum HISSP-FGFER2036' = percentage)
HISSPFGFER2036class1 <- HISSPFGFER2036class %>% select(-c(n,total,percentage))
HcapHISSPFGFER2036<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFER2036class.tsv")
totalHISSPFGFER2036 <- HcapHISSPFGFER2036 %>% left_join(HISSPFGFER2036class, by = "class")
dfHISSPFGFER2036 <- as.data.frame(totalHISSPFGFER2036) 
dfHISSPFGFER2036[is.na(dfHISSPFGFER2036)] <- 0
dfHISSPFGFER2036_new <- dfHISSPFGFER2036 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPFGFER2036 <- dfHISSPFGFER2036_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPFGFER2036, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGFER2036class_TOTAL.tsv")
df_HISSPFGFER20361 <- df_HISSPFGFER2036 %>% mutate('Histoplasma capsulatum HISSP-FGFER2036' = percentage)
df_HISSPFGFER203611 <- df_HISSPFGFER20361 %>% select(-c(n,total,percentage))
##

HISSPFGFAN2059 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPFGFAN2059.out",col_names = F)
colnames(HISSPFGFAN2059) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPFGFAN2059 <- HISSPFGFAN2059 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPFGFAN2059,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGFAN2059_adjust.tsv")
HISSPFGFAN2059 <- dplyr::mutate(HISSPFGFAN2059, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPFGFAN2059class <- HISSPFGFAN2059 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23710091)
write_tsv(HISSPFGFAN2059class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGFAN2059class_Starship.tsv")
HISSPFGFAN2059class <- HISSPFGFAN2059class %>% mutate('Histoplasma capsulatum HISSP-FGFAN2059' = percentage)
HISSPFGFAN2059class1 <- HISSPFGFAN2059class %>% select(-c(n,total,percentage))
HcapHISSPFGFAN2059<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGFAN2059class.tsv")
totalHISSPFGFAN2059 <- HcapHISSPFGFAN2059 %>% left_join(HISSPFGFAN2059class, by = "class")
dfHISSPFGFAN2059 <- as.data.frame(totalHISSPFGFAN2059) 
dfHISSPFGFAN2059[is.na(dfHISSPFGFAN2059)] <- 0
dfHISSPFGFAN2059_new <- dfHISSPFGFAN2059 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPFGFAN2059 <- dfHISSPFGFAN2059_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPFGFAN2059, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGFAN2059class_TOTAL.tsv")
df_HISSPFGFAN20591 <- df_HISSPFGFAN2059 %>% mutate('Histoplasma capsulatum HISSP-FGFAN2059' = percentage)
df_HISSPFGFAN205911 <- df_HISSPFGFAN20591 %>% select(-c(n,total,percentage))
##

HISSPFGBON2001 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPFGBON2001.out",col_names = F)
colnames(HISSPFGBON2001) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPFGBON2001 <- HISSPFGBON2001 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPFGBON2001,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGBON2001_adjust.tsv")
HISSPFGBON2001 <- dplyr::mutate(HISSPFGBON2001, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPFGBON2001class <- HISSPFGBON2001 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/30753906)
write_tsv(HISSPFGBON2001class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGBON2001class_Starship.tsv")
HISSPFGBON2001class <- HISSPFGBON2001class %>% mutate('Histoplasma capsulatum HISSP-FGBON2001' = percentage)
HISSPFGBON2001class1 <- HISSPFGBON2001class %>% select(-c(n,total,percentage))
HcapHISSPFGBON2001<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGBON2001class.tsv")
totalHISSPFGBON2001 <- HcapHISSPFGBON2001 %>% left_join(HISSPFGBON2001class, by = "class")
dfHISSPFGBON2001 <- as.data.frame(totalHISSPFGBON2001) 
dfHISSPFGBON2001[is.na(dfHISSPFGBON2001)] <- 0
dfHISSPFGBON2001_new <- dfHISSPFGBON2001 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPFGBON2001 <- dfHISSPFGBON2001_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPFGBON2001, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGBON2001class_TOTAL.tsv")
df_HISSPFGBON20011 <- df_HISSPFGBON2001 %>% mutate('Histoplasma capsulatum HISSP-FGBON2001' = percentage)
df_HISSPFGBON200111 <- df_HISSPFGBON20011 %>% select(-c(n,total,percentage))
##

HISSPFGBIK2051 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPFGBIK2051.out",col_names = F)
colnames(HISSPFGBIK2051) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPFGBIK2051 <- HISSPFGBIK2051 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPFGBIK2051,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGBIK2051_adjust.tsv")
HISSPFGBIK2051 <- dplyr::mutate(HISSPFGBIK2051, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPFGBIK2051class <- HISSPFGBIK2051 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23534554)
write_tsv(HISSPFGBIK2051class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGBIK2051class_Starship.tsv")
HISSPFGBIK2051class <- HISSPFGBIK2051class %>% mutate('Histoplasma capsulatum HISSP-FGBIK2051' = percentage)
HISSPFGBIK2051class1 <- HISSPFGBIK2051class %>% select(-c(n,total,percentage))
HcapHISSPFGBIK2051<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-FGBIK2051class.tsv")
totalHISSPFGBIK2051 <- HcapHISSPFGBIK2051 %>% left_join(HISSPFGBIK2051class, by = "class")
dfHISSPFGBIK2051 <- as.data.frame(totalHISSPFGBIK2051) 
dfHISSPFGBIK2051[is.na(dfHISSPFGBIK2051)] <- 0
dfHISSPFGBIK2051_new <- dfHISSPFGBIK2051 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPFGBIK2051 <- dfHISSPFGBIK2051_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPFGBIK2051, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPFGBIK2051class_TOTAL.tsv")
df_HISSPFGBIK20511 <- df_HISSPFGBIK2051 %>% mutate('Histoplasma capsulatum HISSP-FGBIK2051' = percentage)
df_HISSPFGBIK205111 <- df_HISSPFGBIK20511 %>% select(-c(n,total,percentage))
##

HISSPCM6408 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPCM6408.out",col_names = F)
colnames(HISSPCM6408) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPCM6408 <- HISSPCM6408 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPCM6408,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPCM6408_adjust.tsv")
HISSPCM6408 <- dplyr::mutate(HISSPCM6408, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPCM6408class <- HISSPCM6408 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/26008567)
write_tsv(HISSPCM6408class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPCM6408class_Starship.tsv")
HISSPCM6408class <- HISSPCM6408class %>% mutate('Histoplasma capsulatum HISSP-CM6408' = percentage)
HISSPCM6408class1 <- HISSPCM6408class %>% select(-c(n,total,percentage))
HcapHISSPCM6408<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-CM6408class.tsv")
totalHISSPCM6408 <- HcapHISSPCM6408 %>% left_join(HISSPCM6408class, by = "class")
dfHISSPCM6408 <- as.data.frame(totalHISSPCM6408) 
dfHISSPCM6408[is.na(dfHISSPCM6408)] <- 0
dfHISSPCM6408_new <- dfHISSPCM6408 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPCM6408 <- dfHISSPCM6408_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPCM6408, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPCM6408class_TOTAL.tsv")
df_HISSPCM64081 <- df_HISSPCM6408 %>% mutate('Histoplasma capsulatum HISSP-CM6408' = percentage)
df_HISSPCM640811 <- df_HISSPCM64081 %>% select(-c(n,total,percentage))
##

HISSPB05821 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSPB05821.out",col_names = F)
colnames(HISSPB05821) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSPB05821 <- HISSPB05821 %>% mutate(length=endq-beginq +1)
write_tsv(HISSPB05821,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPB05821_adjust.tsv")
HISSPB05821 <- dplyr::mutate(HISSPB05821, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSPB05821class <- HISSPB05821 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/36858421)
write_tsv(HISSPB05821class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPB05821class_Starship.tsv")
HISSPB05821class <- HISSPB05821class %>% mutate('Histoplasma capsulatum HISSP-B05821' = percentage)
HISSPB05821class1 <- HISSPB05821class %>% select(-c(n,total,percentage))
HcapHISSPB05821<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-B05821class.tsv")
totalHISSPB05821 <- HcapHISSPB05821 %>% left_join(HISSPB05821class, by = "class")
dfHISSPB05821 <- as.data.frame(totalHISSPB05821) 
dfHISSPB05821[is.na(dfHISSPB05821)] <- 0
dfHISSPB05821_new <- dfHISSPB05821 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSPB05821 <- dfHISSPB05821_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSPB05821, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSPB05821class_TOTAL.tsv")
df_HISSPB058211 <- df_HISSPB05821 %>% mutate('Histoplasma capsulatum HISSP-B05821' = percentage)
df_HISSPB0582111 <- df_HISSPB058211 %>% select(-c(n,total,percentage))
##

HISSP11571 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSP11571.out",col_names = F)
colnames(HISSP11571) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSP11571 <- HISSP11571 %>% mutate(length=endq-beginq +1)
write_tsv(HISSP11571,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSP11571_adjust.tsv")
HISSP11571 <- dplyr::mutate(HISSP11571, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSP11571class <- HISSP11571 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/21937402)
write_tsv(HISSP11571class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSP11571class_Starship.tsv")
HISSP11571class <- HISSP11571class %>% mutate('Histoplasma capsulatum HISSP-11571' = percentage)
HISSP11571class1 <- HISSP11571class %>% select(-c(n,total,percentage))
HcapHISSP11571<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-11571-Belem1class.tsv")
totalHISSP11571 <- HcapHISSP11571 %>% left_join(HISSP11571class, by = "class")
dfHISSP11571 <- as.data.frame(totalHISSP11571) 
dfHISSP11571[is.na(dfHISSP11571)] <- 0
dfHISSP11571_new <- dfHISSP11571 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSP11571 <- dfHISSP11571_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSP11571, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSP11571class_TOTAL.tsv")
df_HISSP115711 <- df_HISSP11571 %>% mutate('Histoplasma capsulatum HISSP-11571' = percentage)
df_HISSP1157111 <- df_HISSP115711 %>% select(-c(n,total,percentage))
##

HISSP1014 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapHISSP1014.out",col_names = F)
colnames(HISSP1014) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
HISSP1014 <- HISSP1014 %>% mutate(length=endq-beginq +1)
write_tsv(HISSP1014,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSP1014_adjust.tsv")
HISSP1014 <- dplyr::mutate(HISSP1014, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
HISSP1014class <- HISSP1014 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/31744825)
write_tsv(HISSP1014class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSP1014class_Starship.tsv")
HISSP1014class <- HISSP1014class %>% mutate('Histoplasma capsulatum HISSP-1014' = percentage)
HISSP1014class1 <- HISSP1014class %>% select(-c(n,total,percentage))
HcapHISSP1014<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_HISSP-1014-Belem3class.tsv")
totalHISSP1014 <- HcapHISSP1014 %>% left_join(HISSP1014class, by = "class")
dfHISSP1014 <- as.data.frame(totalHISSP1014) 
dfHISSP1014[is.na(dfHISSP1014)] <- 0
dfHISSP1014_new <- dfHISSP1014 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_HISSP1014 <- dfHISSP1014_new %>% select(c(class,n,total,percentage))
write_tsv(df_HISSP1014, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HISSP1014class_TOTAL.tsv")
df_HISSP115711 <- df_HISSP11571 %>% mutate('Histoplasma capsulatum HISSP-11571' = percentage)
df_HISSP1157111 <- df_HISSP115711 %>% select(-c(n,total,percentage))
##

RJ0712 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap0712-RJ.out",col_names = F)
colnames(RJ0712) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
RJ0712 <- RJ0712 %>% mutate(length=endq-beginq +1)
write_tsv(RJ0712,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/RJ0712_adjust.tsv")
RJ0712 <- dplyr::mutate(RJ0712, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
RJ0712class <- RJ0712 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/17244466)
write_tsv(RJ0712class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/RJ0712class_Starship.tsv")
RJ0712class <- RJ0712class %>% mutate('Histoplasma capsulatum 0712-RJ' = percentage)
RJ0712class1 <- RJ0712class %>% select(-c(n,total,percentage))
HcapRJ0712<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_07_12-RJclass.tsv")
totalRJ0712 <- HcapRJ0712 %>% left_join(RJ0712class, by = "class")
dfRJ0712 <- as.data.frame(totalRJ0712) 
dfRJ0712[is.na(dfRJ0712)] <- 0
dfRJ0712_new <- dfRJ0712 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_RJ0712 <- dfRJ0712_new %>% select(c(class,n,total,percentage))
write_tsv(df_RJ0712, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/RJ0712class_TOTAL.tsv")
df_RJ07121 <- df_RJ0712 %>% mutate('Histoplasma capsulatum 07-12-RJ' = percentage)
df_RJ071211 <- df_RJ07121 %>% select(-c(n,total,percentage))
##

H104p06 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap104p06.out",col_names = F)
colnames(H104p06) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H104p06 <- H104p06 %>% mutate(length=endq-beginq +1)
write_tsv(H104p06,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H104p06_adjust.tsv")
H104p06 <- dplyr::mutate(H104p06, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H104p06class <- H104p06 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28442945)
write_tsv(H104p06class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/104p06class_Starship.tsv")
H104p06class <- H104p06class %>% mutate('Histoplasma capsulatum 104_p_06' = percentage)
H104p06class1 <- H104p06class %>% select(-c(n,total,percentage))
HcapH104p06<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_104_p_06class.tsv")
totalH104p06 <- HcapH104p06 %>% left_join(H104p06class, by = "class")
dfH104p06 <- as.data.frame(totalH104p06) 
dfH104p06[is.na(dfH104p06)] <- 0
dfH104p06_new <- dfH104p06 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H104p06 <- dfH104p06_new %>% select(c(class,n,total,percentage))
write_tsv(df_H104p06, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/104p06class_TOTAL.tsv")
df_H104p061 <- df_H104p06 %>% mutate('Histoplasma capsulatum 104_p_06' = percentage)
df_H104p0611 <- df_H104p061 %>% select(-c(n,total,percentage))
##

H117p12 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap117p12.out",col_names = F)
colnames(H117p12) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H117p12 <- H117p12 %>% mutate(length=endq-beginq +1)
write_tsv(H117p12,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H117p12_adjust.tsv")
H117p12 <- dplyr::mutate(H117p12, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H117p12class <- H117p12 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28587451)
write_tsv(H117p12class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H117p12class_Starship.tsv")
H117p12class <- H117p12class %>% mutate('Histoplasma capsulatum 117_p_12' = percentage)
H117p12class1 <- H117p12class %>% select(-c(n,total,percentage))
HcapH117p12<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_117_p_12class.tsv")
totalH117p12 <- HcapH117p12 %>% left_join(H117p12class, by = "class")
dfH117p12 <- as.data.frame(totalH117p12) 
dfH117p12[is.na(dfH117p12)] <- 0
dfH117p12_new <- dfH117p12 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H117p12 <- dfH117p12_new %>% select(c(class,n,total,percentage))
write_tsv(df_H117p12, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H117p12class_TOTAL.tsv")
df_H117p121 <- df_H117p12 %>% mutate('Histoplasma capsulatum 117_p_12' = percentage)
df_H117p1211 <- df_H117p121 %>% select(-c(n,total,percentage))
##

H122p10B <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap122p10B.out",col_names = F)
colnames(H122p10B) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H122p10B <- H122p10B %>% mutate(length=endq-beginq +1)
write_tsv(H122p10B,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H122p10B_adjust.tsv")
H122p10B <- dplyr::mutate(H122p10B, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H122p10Bclass <- H122p10B %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28585056)
write_tsv(H122p10Bclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H122p10Bclass_Starship.tsv")
H122p10Bclass <- H122p10Bclass %>% mutate('Histoplasma capsulatum 122_p_10_B' = percentage)
H122p10Bclass1 <- H122p10Bclass %>% select(-c(n,total,percentage))
HcapH122p10B<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_122_p_10_Bclass.tsv")
totalH122p10B <- HcapH122p10B %>% left_join(H122p10Bclass, by = "class")
dfH122p10B <- as.data.frame(totalH122p10B) 
dfH122p10B[is.na(dfH122p10B)] <- 0
dfH122p10B_new <- dfH122p10B %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H122p10B <- dfH122p10B_new %>% select(c(class,n,total,percentage))
write_tsv(df_H122p10B, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H122p10Bclass_TOTAL.tsv")
df_H122p10B1 <- df_H122p10B %>% mutate('Histoplasma capsulatum 122_p_10B' = percentage)
df_H122p10B11 <- df_H122p10B1 %>% select(-c(n,total,percentage))
##

H136P07 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap136P07.out",col_names = F)
colnames(H136P07) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H136P07 <- H136P07 %>% mutate(length=endq-beginq +1)
write_tsv(H136P07,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H136P07_adjust.tsv")
H136P07 <- dplyr::mutate(H136P07, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H136P07class <- H136P07 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28536980)
write_tsv(H136P07class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H136P07class_Starship.tsv")
H136P07class <- H136P07class %>% mutate('Histoplasma capsulatum 136_P_07' = percentage)
H136P07class1 <- H136P07class %>% select(-c(n,total,percentage))
HcapH136P07<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_136_P_07class.tsv")
totalH136P07 <- HcapH136P07 %>% left_join(H136P07class, by = "class")
dfH136P07 <- as.data.frame(totalH136P07) 
dfH136P07[is.na(dfH136P07)] <- 0
dfH136P07_new <- dfH136P07 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H136P07 <- dfH136P07_new %>% select(c(class,n,total,percentage))
write_tsv(df_H136P07, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H136P07class_TOTAL.tsv")
df_H136P071 <- df_H136P07 %>% mutate('Histoplasma capsulatum 136_P_07' = percentage)
df_H136P0711 <- df_H136P071 %>% select(-c(n,total,percentage))
##

H144p08 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap144p08.out",col_names = F)
colnames(H144p08 ) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H144p08  <- H144p08  %>% mutate(length=endq-beginq +1)
write_tsv(H144p08 ,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H144p08_adjust.tsv")
H144p08 <- dplyr::mutate(H144p08, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H144p08class <- H144p08 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28329803)
write_tsv(H144p08class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H144p08class_Starship.tsv")
H144p08class <- H144p08class %>% mutate('Histoplasma capsulatum 144_p_08' = percentage)
H144p08class1 <- H144p08class %>% select(-c(n,total,percentage))
HcapH144p08<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_144_p_08class.tsv")
totalH144p08 <- HcapH144p08 %>% left_join(H144p08class, by = "class")
dfH144p08 <- as.data.frame(totalH144p08) 
dfH144p08[is.na(dfH144p08)] <- 0
dfH144p08_new <- dfH144p08 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H144p08 <- dfH144p08_new %>% select(c(class,n,total,percentage))
write_tsv(df_H144p08, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H144p08class_TOTAL.tsv")
df_H144p081 <- df_H144p08 %>% mutate('Histoplasma capsulatum 144_p_08' = percentage)
df_H144p0811 <- df_H144p081 %>% select(-c(n,total,percentage))
##

H1517p17 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap1517p17.out",col_names = F)
colnames(H1517p17) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H1517p17 <- H1517p17 %>% mutate(length=endq-beginq +1)
write_tsv(H1517p17 ,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H1517p17_adjust.tsv")
H1517p17 <- dplyr::mutate(H1517p17, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H1517p17class <- H1517p17 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28646603)
write_tsv(H1517p17class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H1517p17class_Starship.tsv")
H1517p17class <- H1517p17class %>% mutate('Histoplasma capsulatum 1517_p_17' = percentage)
H1517p17class1 <- H1517p17class %>% select(-c(n,total,percentage))
HcapH1517p17 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_1517_p_17class.tsv")
totalH1517p17 <- HcapH1517p17 %>% left_join(H1517p17class, by = "class")
dfH1517p17 <- as.data.frame(totalH1517p17) 
dfH1517p17[is.na(dfH1517p17)] <- 0
dfH1517p17_new <- dfH1517p17 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H1517p17 <- dfH1517p17_new %>% select(c(class,n,total,percentage))
write_tsv(df_H1517p17, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H1517p17class_TOTAL.tsv")
df_H1517p171 <- df_H1517p17 %>% mutate('Histoplasma capsulatum 1517_p_17' = percentage)
df_H1517p1711 <- df_H1517p171 %>% select(-c(n,total,percentage))
##

H256P18 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap256P18.out",col_names = F)
colnames(H256P18) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H256P18 <- H256P18 %>% mutate(length=endq-beginq +1)
write_tsv(H256P18,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H256P18_adjust.tsv")
H256P18 <- dplyr::mutate(H256P18, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H256P18class <- H256P18 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28534158)
write_tsv(H256P18class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H256P18class_Starship.tsv")
H256P18class <- H256P18class %>% mutate('Histoplasma capsulatum 256_P_18' = percentage)
H256P18class1 <- H256P18class %>% select(-c(n,total,percentage))
HcapH256P18 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_256_P_18class.tsv")
totalH256P18 <- HcapH256P18 %>% left_join(H256P18class, by = "class")
dfH256P18 <- as.data.frame(totalH256P18) 
dfH256P18[is.na(dfH256P18)] <- 0
dfH256P18_new <- dfH256P18 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H256P18 <- dfH256P18_new %>% select(c(class,n,total,percentage))
write_tsv(df_H256P18, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H256P18class_TOTAL.tsv")
df_H256P181 <- df_H256P18 %>% mutate('Histoplasma capsulatum 256_P_18' = percentage)
df_H256P1811 <- df_H256P181 %>% select(-c(n,total,percentage))
##

H316p10 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap316p10.out",col_names = F)
colnames(H316p10) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H316p10 <- H316p10 %>% mutate(length=endq-beginq +1)
write_tsv(H316p10,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H316p10_adjust.tsv")
H316p10 <- dplyr::mutate(H316p10, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H316p10class <- H316p10 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28625580)
write_tsv(H316p10class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H316p10class_Starship.tsv")
H316p10class <- H316p10class %>% mutate('Histoplasma capsulatum 316_p_10' = percentage)
H316p10class1 <- H316p10class %>% select(-c(n,total,percentage))
HcapH316p10 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_316_p_10class.tsv")
totalH316p10 <- HcapH316p10 %>% left_join(H316p10class, by = "class")
dfH316p10 <- as.data.frame(totalH316p10) 
dfH316p10[is.na(dfH316p10)] <- 0
dfH316p10_new <- dfH316p10 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H316p10 <- dfH316p10_new %>% select(c(class,n,total,percentage))
write_tsv(df_H316p10, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H316p10class_TOTAL.tsv")
df_H316p101 <- df_H316p10 %>% mutate('Histoplasma capsulatum 316_p_10' = percentage)
df_H316p1011 <- df_H316p101 %>% select(-c(n,total,percentage))
##

H327P12 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap327P12.out",col_names = F)
colnames(H327P12) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H327P12 <- H327P12 %>% mutate(length=endq-beginq +1)
write_tsv(H327P12,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H327P12_adjust.tsv")
H327P12 <- dplyr::mutate(H327P12, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H327P12class <- H327P12 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28431708)
write_tsv(H327P12class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H327P12class_Starship.tsv")
H327P12class <- H327P12class %>% mutate('Histoplasma capsulatum 327_P_12' = percentage)
H327P12class1 <- H327P12class %>% select(-c(n,total,percentage))
HcapH327P12 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_327_P_12class.tsv")
totalH327P12 <- HcapH327P12 %>% left_join(H327P12class, by = "class")
dfH327P12 <- as.data.frame(totalH327P12) 
dfH327P12[is.na(dfH327P12)] <- 0
dfH327P12_new <- dfH327P12 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H327P12 <- dfH327P12_new %>% select(c(class,n,total,percentage))
write_tsv(df_H327P12, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H327P12class_TOTAL.tsv")
df_H327P121 <- df_H327P12 %>% mutate('Histoplasma capsulatum 327_P_12' = percentage)
df_H327P1211 <- df_H327P121 %>% select(-c(n,total,percentage))
##

H343p18 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap343p18.out",col_names = F)
colnames(H343p18) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H343p18 <- H343p18 %>% mutate(length=endq-beginq +1)
write_tsv(H343p18,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H343p18_adjust.tsv")
H343p18 <- dplyr::mutate(H343p18, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H343p18class <- H343p18 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28732569)
write_tsv(H343p18class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H343p18class_Starship.tsv")
H343p18class <- H343p18class %>% mutate('Histoplasma capsulatum 343_p_18' = percentage)
H343p18class1 <- H343p18class %>% select(-c(n,total,percentage))
HcapH343p18 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_343_p_18class.tsv")
totalH343p18 <- HcapH343p18 %>% left_join(H343p18class, by = "class")
dfH343p18 <- as.data.frame(totalH343p18) 
dfH343p18[is.na(dfH343p18)] <- 0
dfH343p18_new <- dfH343p18 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H343p18 <- dfH343p18_new %>% select(c(class,n,total,percentage))
write_tsv(df_H343p18, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H343p18class_TOTAL.tsv")
df_H343p181 <- df_H343p18 %>% mutate('Histoplasma capsulatum 343_p_18' = percentage)
df_H343p1811 <- df_H343p181 %>% select(-c(n,total,percentage))
##

H388p11 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Histocap388p11.out",col_names = F)
colnames(H388p11) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
H388p11 <- H388p11 %>% mutate(length=endq-beginq +1)
write_tsv(H388p11,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H388p11_adjust.tsv")
H388p11 <- dplyr::mutate(H388p11, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
H388p11class <- H388p11 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28579830)
write_tsv(H388p11class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H388p11class_Starship.tsv")
H388p11class <- H388p11class %>% mutate('Histoplasma capsulatum 388_p_11' = percentage)
H388p11class1 <- H388p11class %>% select(-c(n,total,percentage))
HcapH388p11 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_388_p_11class.tsv")
totalH388p11 <- HcapH388p11 %>% left_join(H388p11class, by = "class")
dfH388p11 <- as.data.frame(totalH388p11) 
dfH388p11[is.na(dfH388p11)] <- 0
dfH388p11_new <- dfH388p11 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_H388p11 <- dfH388p11_new %>% select(c(class,n,total,percentage))
write_tsv(df_H388p11, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/H388p11class_TOTAL.tsv")
df_H388p111 <- df_H388p11 %>% mutate('Histoplasma capsulatum 388_p_11' = percentage)
df_H388p1111 <- df_H388p111 %>% select(-c(n,total,percentage))
##

ES283Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapES283Z.out",col_names = F)
colnames(ES283Z) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
ES283Z <- ES283Z %>% mutate(length=endq-beginq +1)
write_tsv(ES283Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES283Z_adjust.tsv")
ES283Z <- dplyr::mutate(ES283Z, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
ES283Zclass <- ES283Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/33292477)
write_tsv(ES283Zclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES283Zclass_Starship.tsv")
ES283Zclass <- ES283Zclass %>% mutate('Histoplasma capsulatum ES2_83Z' = percentage)
ES283Zclass1 <- ES283Zclass %>% select(-c(n,total,percentage))
HcapES283Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_83Zclass.tsv")
totalES283Z <- HcapES283Z %>% left_join(ES283Zclass, by = "class")
dfES283Z <- as.data.frame(totalES283Z) 
dfES283Z[is.na(dfES283Z)] <- 0
dfES283Z_new <- dfES283Z %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_ES283Z <- dfES283Z_new %>% select(c(class,n,total,percentage))
write_tsv(df_ES283Z, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES283Zclass_TOTAL.tsv")
df_ES283Z1 <- df_ES283Z %>% mutate('Histoplasma capsulatum ES2_83Z' = percentage)
df_ES283Z11 <- df_ES283Z1 %>% select(-c(n,total,percentage))
##

ES285Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapES285Z.out",col_names = F)
colnames(ES285Z) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
ES285Z <- ES285Z %>% mutate(length=endq-beginq +1)
write_tsv(ES285Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES285Z_adjust.tsv")
ES285Z <- dplyr::mutate(ES285Z, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
ES285Zclass <- ES285Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28184349)
write_tsv(ES285Zclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES2_85Zclass_Starship.tsv")
ES285Zclass <- ES285Zclass %>% mutate('Histoplasma capsulatum ES2_85Z' = percentage)
ES285Zclass1 <- ES285Zclass %>% select(-c(n,total,percentage))
HcapES285Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_85Zclass.tsv")
totalES285Z <- HcapES285Z %>% left_join(ES285Zclass, by = "class")
dfES285Z <- as.data.frame(totalES285Z) 
dfES285Z[is.na(dfES285Z)] <- 0
dfES285Z_new <- dfES285Z %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_ES285Z <- dfES285Z_new %>% select(c(class,n,total,percentage))
write_tsv(df_ES285Z, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES285Zclass_TOTAL.tsv")
df_ES285Z1 <- df_ES285Z %>% mutate('Histoplasma capsulatum ES2_85Z' = percentage)
df_ES285Z11 <- df_ES285Z1 %>% select(-c(n,total,percentage))
##

ES286Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapES286Z.out",col_names = F)
colnames(ES286Z) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
ES286Z <- ES286Z %>% mutate(length=endq-beginq +1)
write_tsv(ES286Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES286Z_adjust.tsv")
ES286Z <- dplyr::mutate(ES286Z, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
ES286Zclass <- ES286Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/40716201)
write_tsv(ES286Zclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES2_86Zclass_Starship.tsv")
ES286Zclass <- ES286Zclass %>% mutate('Histoplasma capsulatum ES2_86Z' = percentage)
ES286Zclass1 <- ES286Zclass %>% select(-c(n,total,percentage))
HcapES286Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_86Zclass.tsv")
totalES286Z <- HcapES286Z %>% left_join(ES286Zclass, by = "class")
dfES286Z <- as.data.frame(totalES286Z) 
dfES286Z[is.na(dfES286Z)] <- 0
dfES286Z_new <- dfES286Z %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_ES286Z <- dfES286Z_new %>% select(c(class,n,total,percentage))
write_tsv(df_ES286Z, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES286Zclass_TOTAL.tsv")
df_ES286Z1 <- df_ES286Z %>% mutate('Histoplasma capsulatum ES2_86Z' = percentage)
df_ES286Z11 <- df_ES286Z1 %>% select(-c(n,total,percentage))
##

ES288Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapES288Z.out",col_names = F)
colnames(ES288Z) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
ES288Z <- ES288Z %>% mutate(length=endq-beginq +1)
write_tsv(ES288Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES288Z_adjust.tsv")
ES288Z <- dplyr::mutate(ES288Z, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
ES288Zclass <- ES288Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/25109007)
write_tsv(ES288Zclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES2_88Zclass_Starship.tsv")
ES288Zclass <- ES288Zclass %>% mutate('Histoplasma capsulatum ES2_88Z' = percentage)
ES288Zclass1 <- ES288Zclass %>% select(-c(n,total,percentage))
HcapES288Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_88Zclass.tsv")
totalES288Z <- HcapES288Z %>% left_join(ES288Zclass, by = "class")
dfES288Z <- as.data.frame(totalES288Z) 
dfES288Z[is.na(dfES288Z)] <- 0
dfES288Z_new <- dfES288Z %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_ES288Z <- dfES288Z_new %>% select(c(class,n,total,percentage))
write_tsv(df_ES288Z, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES288Zclass_TOTAL.tsv")
df_ES288Z1 <- df_ES288Z %>% mutate('Histoplasma capsulatum ES2_88Z' = percentage)
df_ES288Z11 <- df_ES288Z1 %>% select(-c(n,total,percentage))
##

ES289Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapES289Z.out",col_names = F)
colnames(ES289Z) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
ES289Z <- ES289Z %>% mutate(length=endq-beginq +1)
write_tsv(ES289Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES289Z_adjust.tsv")
ES289Z <- dplyr::mutate(ES289Z, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
ES289Zclass <- ES289Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/34058413)
write_tsv(ES289Zclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES2_89Zclass_Starship.tsv")
ES289Zclass <- ES289Zclass %>% mutate('Histoplasma capsulatum ES2_89Z' = percentage)
ES289Zclass1 <- ES289Zclass %>% select(-c(n,total,percentage))
HcapES289Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_89Zclass.tsv")
totalES289Z <- HcapES289Z %>% left_join(ES289Zclass, by = "class")
dfES289Z <- as.data.frame(totalES289Z) 
dfES289Z[is.na(dfES289Z)] <- 0
dfES289Z_new <- dfES289Z %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_ES289Z <- dfES289Z_new %>% select(c(class,n,total,percentage))
write_tsv(df_ES289Z, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES289Zclass_TOTAL.tsv")
df_ES289Z1 <- df_ES289Z %>% mutate('Histoplasma capsulatum ES2_89Z' = percentage)
df_ES289Z11 <- df_ES289Z1 %>% select(-c(n,total,percentage))
##

ES290Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapES290Z.out",col_names = F)
colnames(ES290Z) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
ES290Z <- ES290Z %>% mutate(length=endq-beginq +1)
write_tsv(ES290Z,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES290Z_adjust.tsv")
ES290Z <- dplyr::mutate(ES290Z, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
ES290Zclass <- ES290Z %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/36554438)
write_tsv(ES290Zclass, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES2_90Zclass_Starship.tsv")
ES290Zclass <- ES290Zclass %>% mutate('Histoplasma capsulatum ES2_90Z' = percentage)
ES290Zclass1 <- ES290Zclass %>% select(-c(n,total,percentage))
HcapES290Z <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_ES2_90Zclass.tsv")
totalES290Z <- HcapES290Z %>% left_join(ES290Zclass, by = "class")
dfES290Z <- as.data.frame(totalES290Z) 
dfES290Z[is.na(dfES290Z)] <- 0
dfES290Z_new <- dfES290Z %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_ES290Z <- dfES290Z_new %>% select(c(class,n,total,percentage))
write_tsv(df_ES290Z, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/ES290Zclass_TOTAL.tsv")
df_ES290Z1 <- df_ES290Z %>% mutate('Histoplasma capsulatum ES2_90Z' = percentage)
df_ES290Z11 <- df_ES290Z1 %>% select(-c(n,total,percentage))
##

SA15 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapSA15.out",col_names = F)
colnames(SA15) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
SA15 <- SA15 %>% mutate(length=endq-beginq +1)
write_tsv(SA15,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SA15_adjust.tsv")
SA15 <- dplyr::mutate(SA15, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
SA15class <- SA15 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28298332)
write_tsv(SA15class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SA15class_Starship.tsv")
SA15class <- SA15class %>% mutate('Histoplasma capsulatum SA15' = percentage)
SA15class1 <- SA15class %>% select(-c(n,total,percentage))
HcapSA15 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_SA15class.tsv")
totalSA15 <- HcapSA15 %>% left_join(SA15class, by = "class")
dfSA15 <- as.data.frame(totalSA15) 
dfSA15[is.na(dfSA15)] <- 0
dfSA15_new <- dfSA15 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_SA15 <- dfSA15_new %>% select(c(class,n,total,percentage))
write_tsv(df_SA15, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/SA15class_TOTAL.tsv")
df_SA151 <- df_ES290Z %>% mutate('Histoplasma capsulatum SA15' = percentage)
df_SA1511 <- df_ES290Z1 %>% select(-c(n,total,percentage))
##

CB053 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB053.out",col_names = F)
colnames(CB053) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB053 <- CB053 %>% mutate(length=endq-beginq +1)
write_tsv(CB053,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB053_adjust.tsv")
CB053 <- dplyr::mutate(CB053, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB053class <- CB053 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23649534)
write_tsv(CB053class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB053class_Starship.tsv")
CB053class <- CB053class %>% mutate('Histoplasma capsulatum CB053' = percentage)
CB053class1 <- CB053class %>% select(-c(n,total,percentage))
HcapCB053<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB053-2class.tsv")
totalCB053 <- HcapCB053 %>% left_join(CB053class, by = "class")
dfCB053 <- as.data.frame(totalCB053) 
dfCB053[is.na(dfCB053)] <- 0
dfCB053_new <- dfCB053 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB053 <- dfCB053_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB053, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB053class_TOTAL.tsv")
df_CB0531 <- df_CB053 %>% mutate('Histoplasma capsulatum CB053' = percentage)
df_CB05311 <- df_CB0531 %>% select(-c(n,total,percentage))
##

CB055 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB055.out",col_names = F)
colnames(CB055) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB055 <- CB055 %>% mutate(length=endq-beginq +1)
write_tsv(CB055,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB055_adjust.tsv")
CB055 <- dplyr::mutate(CB055, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB055class <- CB055 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/23847241)
write_tsv(CB055class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB055class_Starship.tsv")
CB055class <- CB055class %>% mutate('Histoplasma capsulatum CB055' = percentage)
CB055class1 <- CB055class %>% select(-c(n,total,percentage))
HcapCB055<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB055-2class.tsv")
totalCB055 <- HcapCB055 %>% left_join(CB055class, by = "class")
dfCB055 <- as.data.frame(totalCB055) 
dfCB055[is.na(dfCB055)] <- 0
dfCB055_new <- dfCB055 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB055 <- dfCB055_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB055, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB055class_TOTAL.tsv")
df_CB0551 <- df_CB055 %>% mutate('Histoplasma capsulatum CB055' = percentage)
df_CB05511 <- df_CB0551 %>% select(-c(n,total,percentage))
##

CB062 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB062.out",col_names = F)
colnames(CB062) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB062 <- CB062 %>% mutate(length=endq-beginq +1)
write_tsv(CB062,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB062_adjust.tsv")
CB062 <- dplyr::mutate(CB062, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB062class <- CB062 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28739478)
write_tsv(CB062class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB062class_Starship.tsv")
CB062class <- CB062class %>% mutate('Histoplasma capsulatum CB062' = percentage)
CB062class1 <- CB062class %>% select(-c(n,total,percentage))
HcapCB062<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB062-2class.tsv")
totalCB062 <- HcapCB062 %>% left_join(CB062class, by = "class")
dfCB062 <- as.data.frame(totalCB062) 
dfCB062[is.na(dfCB062)] <- 0
dfCB062_new <- dfCB062 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB062 <- dfCB062_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB062, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB062class_TOTAL.tsv")
df_CB0621 <- df_CB062 %>% mutate('Histoplasma capsulatum CB062' = percentage)
df_CB06211 <- df_CB0621 %>% select(-c(n,total,percentage))
##

CB063 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB063.out",col_names = F)
colnames(CB063) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB063 <- CB063 %>% mutate(length=endq-beginq +1)
write_tsv(CB063,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB063_adjust.tsv")
CB063 <- dplyr::mutate(CB063, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB063class <- CB063 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28655085)
write_tsv(CB063class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB063class_Starship.tsv")
CB063class <- CB063class %>% mutate('Histoplasma capsulatum CB063' = percentage)
CB063class1 <- CB063class %>% select(-c(n,total,percentage))
HcapCB063<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB063-2class.tsv")
totalCB063 <- HcapCB063 %>% left_join(CB063class, by = "class")
dfCB063 <- as.data.frame(totalCB063) 
dfCB063[is.na(dfCB063)] <- 0
dfCB063_new <- dfCB063 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB063 <- dfCB063_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB063, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB063class_TOTAL.tsv")
df_CB0631 <- df_CB063 %>% mutate('Histoplasma capsulatum CB063' = percentage)
df_CB06311 <- df_CB0631 %>% select(-c(n,total,percentage))
##

CB064 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB064.out",col_names = F)
colnames(CB064) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB064 <- CB064 %>% mutate(length=endq-beginq +1)
write_tsv(CB064,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB064_adjust.tsv")
CB064 <- dplyr::mutate(CB064, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB064class <- CB064 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28635935)
write_tsv(CB064class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB064class_Starship.tsv")
CB064class <- CB064class %>% mutate('Histoplasma capsulatum CB064' = percentage)
CB064class1 <- CB064class %>% select(-c(n,total,percentage))
HcapCB064<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB064-2class.tsv")
totalCB064 <- HcapCB064 %>% left_join(CB064class, by = "class")
dfCB064 <- as.data.frame(totalCB064) 
dfCB064[is.na(dfCB064)] <- 0
dfCB064_new <- dfCB064 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB064 <- dfCB064_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB064, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB064class_TOTAL.tsv")
df_CB0641 <- df_CB064 %>% mutate('Histoplasma capsulatum CB064' = percentage)
df_CB06411 <- df_CB0641 %>% select(-c(n,total,percentage))
##

CB065 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB065.out",col_names = F)
colnames(CB065) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB065 <- CB065 %>% mutate(length=endq-beginq +1)
write_tsv(CB065,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB065_adjust.tsv")
CB065 <- dplyr::mutate(CB065, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB065class <- CB065 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/29229412)
write_tsv(CB065class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB065class_Starship.tsv")
CB065class <- CB065class %>% mutate('Histoplasma capsulatum CB065' = percentage)
CB065class1 <- CB065class %>% select(-c(n,total,percentage))
HcapCB065<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB065-2class.tsv")
totalCB065 <- HcapCB065 %>% left_join(CB065class, by = "class")
dfCB065 <- as.data.frame(totalCB065) 
dfCB065[is.na(dfCB065)] <- 0
dfCB065_new <- dfCB065 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB065 <- dfCB065_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB065, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB065class_TOTAL.tsv")
df_CB0651 <- df_CB065 %>% mutate('Histoplasma capsulatum CB065' = percentage)
df_CB06511 <- df_CB0651 %>% select(-c(n,total,percentage))
##

CB066 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB066.out",col_names = F)
colnames(CB066) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB066 <- CB066 %>% mutate(length=endq-beginq +1)
write_tsv(CB066,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB066_adjust.tsv")
CB066 <- dplyr::mutate(CB066, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB066class <- CB066 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28630619)
write_tsv(CB066class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB066class_Starship.tsv")
CB066class <- CB066class %>% mutate('Histoplasma capsulatum CB066' = percentage)
CB066class1 <- CB066class %>% select(-c(n,total,percentage))
HcapCB066<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB066-2class.tsv")
totalCB066 <- HcapCB066 %>% left_join(CB066class, by = "class")
dfCB066 <- as.data.frame(totalCB066) 
dfCB066[is.na(dfCB066)] <- 0
dfCB066_new <- dfCB066 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB066 <- dfCB066_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB066, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB066class_TOTAL.tsv")
df_CB0661 <- df_CB066 %>% mutate('Histoplasma capsulatum CB066' = percentage)
df_CB06611 <- df_CB0661 %>% select(-c(n,total,percentage))
##

CB168 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB168.out",col_names = F)
colnames(CB168) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB168 <- CB168 %>% mutate(length=endq-beginq +1)
write_tsv(CB168,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB168_adjust.tsv")
CB168 <- dplyr::mutate(CB168, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB168class <- CB168 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28670342)
write_tsv(CB168class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB168class_Starship.tsv")
CB168class <- CB168class %>% mutate('Histoplasma capsulatum CB168' = percentage)
CB168class1 <- CB168class %>% select(-c(n,total,percentage))
HcapCB168<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB168-2class.tsv")
totalCB168 <- HcapCB168 %>% left_join(CB168class, by = "class")
dfCB168 <- as.data.frame(totalCB168) 
dfCB168[is.na(dfCB168)] <- 0
dfCB168_new <- dfCB168 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB168 <- dfCB168_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB168, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB168class_TOTAL.tsv")
df_CB1681 <- df_CB168 %>% mutate('Histoplasma capsulatum CB168' = percentage)
df_CB16811 <- df_CB1681 %>% select(-c(n,total,percentage))
##

CB174 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB174.out",col_names = F)
colnames(CB174) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB174 <- CB174 %>% mutate(length=endq-beginq +1)
write_tsv(CB174,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB174_adjust.tsv")
CB174 <- dplyr::mutate(CB174, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB174class <- CB174 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28472442)
write_tsv(CB174class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB174class_Starship.tsv")
CB174class <- CB174class %>% mutate('Histoplasma capsulatum CB174' = percentage)
CB174class1 <- CB174class %>% select(-c(n,total,percentage))
HcapCB174<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB174-2class.tsv")
totalCB174 <- HcapCB174 %>% left_join(CB174class, by = "class")
dfCB174 <- as.data.frame(totalCB174) 
dfCB174[is.na(dfCB174)] <- 0
dfCB174_new <- dfCB174 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB174 <- dfCB174_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB174, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB174class_TOTAL.tsv")
df_CB1741 <- df_CB174 %>% mutate('Histoplasma capsulatum CB174' = percentage)
df_CB17411 <- df_CB1741 %>% select(-c(n,total,percentage))
##

CB180 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB180.out",col_names = F)
colnames(CB180) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB180 <- CB180 %>% mutate(length=endq-beginq +1)
write_tsv(CB180,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB180_adjust.tsv")
CB180 <- dplyr::mutate(CB180, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB180class <- CB180 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28608116)
write_tsv(CB180class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB180class_Starship.tsv")
CB180class <- CB180class %>% mutate('Histoplasma capsulatum CB180' = percentage)
CB180class1 <- CB180class %>% select(-c(n,total,percentage))
HcapCB180<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB180-2class.tsv")
totalCB180 <- HcapCB180 %>% left_join(CB180class, by = "class")
dfCB180 <- as.data.frame(totalCB180) 
dfCB180[is.na(dfCB180)] <- 0
dfCB180_new <- dfCB180 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB180 <- dfCB180_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB180, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB180class_TOTAL.tsv")
df_CB1801 <- df_CB180 %>% mutate('Histoplasma capsulatum CB180' = percentage)
df_CB18011 <- df_CB1801 %>% select(-c(n,total,percentage))
##

CB186 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB186.out",col_names = F)
colnames(CB186) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB186 <- CB186 %>% mutate(length=endq-beginq +1)
write_tsv(CB186,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB186_adjust.tsv")
CB186 <- dplyr::mutate(CB186, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB186class <- CB186 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28220179)
write_tsv(CB186class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB186class_Starship.tsv")
CB186class <- CB186class %>% mutate('Histoplasma capsulatum CB186' = percentage)
CB186class1 <- CB186class %>% select(-c(n,total,percentage))
HcapCB186<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB186-2class.tsv")
totalCB186 <- HcapCB186 %>% left_join(CB186class, by = "class")
dfCB186 <- as.data.frame(totalCB186) 
dfCB186[is.na(dfCB186)] <- 0
dfCB186_new <- dfCB186 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB186 <- dfCB186_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB186, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB186class_TOTAL.tsv")
df_CB1861 <- df_CB186 %>% mutate('Histoplasma capsulatum CB186' = percentage)
df_CB18611 <- df_CB1861 %>% select(-c(n,total,percentage))
##

CB192 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapCB192.out",col_names = F)
colnames(CB192) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
CB192 <- CB192 %>% mutate(length=endq-beginq +1)
write_tsv(CB192,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB192_adjust.tsv")
CB192 <- dplyr::mutate(CB192, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
CB192class <- CB192 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28781169)
write_tsv(CB192class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB192class_Starship.tsv")
CB192class <- CB192class %>% mutate('Histoplasma capsulatum CB192' = percentage)
CB192class1 <- CB192class %>% select(-c(n,total,percentage))
HcapCB192<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_CB192-2class.tsv")
totalCB192 <- HcapCB192 %>% left_join(CB192class, by = "class")
dfCB192 <- as.data.frame(totalCB192) 
dfCB192[is.na(dfCB192)] <- 0
dfCB192_new <- dfCB192 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_CB192 <- dfCB192_new %>% select(c(class,n,total,percentage))
write_tsv(df_CB192, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/CB192class_TOTAL.tsv")
df_CB1921 <- df_CB192 %>% mutate('Histoplasma capsulatum CB192' = percentage)
df_CB19211 <- df_CB1921 %>% select(-c(n,total,percentage))
##

S11 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapWGSS11.out",col_names = F)
colnames(S11) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
S11 <- S11 %>% mutate(length=endq-beginq +1)
write_tsv(S11,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S11_adjust.tsv")
S11 <- dplyr::mutate(S11, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
S11class <- S11 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28690043)
write_tsv(S11class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S11class_Starship.tsv")
S11class <- S11class %>% mutate('Histoplasma capsulatum WGSS11' = percentage)
S11class1 <- S11class %>% select(-c(n,total,percentage))
HcapS11<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11class.tsv")
totalS11 <- HcapS11 %>% left_join(S11class, by = "class")
dfS11 <- as.data.frame(totalS11) 
dfS11[is.na(dfS11)] <- 0
dfS11_new <- dfS11 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_S11 <- dfS11_new %>% select(c(class,n,total,percentage))
write_tsv(df_S11, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S11class_TOTAL.tsv")
df_S111 <- df_S11 %>% mutate('Histoplasma capsulatum S11' = percentage)
df_S1111 <- df_S111 %>% select(-c(n,total,percentage))
##

S14 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapWGSS14.out",col_names = F)
colnames(S14) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
S14 <- S14 %>% mutate(length=endq-beginq +1)
write_tsv(S14,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S14_adjust.tsv")
S14 <- dplyr::mutate(S14, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
S14class <- S14 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28326430)
write_tsv(S14class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S14class_Starship.tsv")
S14class <- S14class %>% mutate('Histoplasma capsulatum WGSS14' = percentage)
S14class1 <- S14class %>% select(-c(n,total,percentage))
HcapS14<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S14class.tsv")
totalS14 <- HcapS14 %>% left_join(S14class, by = "class")
dfS14 <- as.data.frame(totalS14) 
dfS14[is.na(dfS14)] <- 0
dfS14_new <- dfS14 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_S14 <- dfS14_new %>% select(c(class,n,total,percentage))
write_tsv(df_S14, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S14class_TOTAL.tsv")
df_S141 <- df_S14 %>% mutate('Histoplasma capsulatum S14' = percentage)
df_S1411 <- df_S141 %>% select(-c(n,total,percentage))
##

S16 <- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/HistocapWGSS16.out",col_names = F)
colnames(S16) <-c("SampleID","prog.","compare","beginq","endq","leftq","strand","dot","Target","class","leftr","div")
S16 <- S16 %>% mutate(length=endq-beginq +1)
write_tsv(S16,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S16_adjust.tsv")
S16 <- dplyr::mutate(S16, ID = row_number())

#the number used in the percentage=total/n is the genome length. Make sure you know the genome length of all your isolates.
S16class <- S16 %>% group_by(class) %>%
  summarise(n=n_distinct(ID),
            total=sum(length),
            percentage=total/28332259)
write_tsv(S16class, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S16class_Starship.tsv")
S16class <- S16class %>% mutate('Histoplasma capsulatum WGSS16' = percentage)
S16class1 <- S16class %>% select(-c(n,total,percentage))
HcapS16<- read_table("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S16class.tsv")
totalS16 <- HcapS16 %>% left_join(S16class, by = "class")
dfS16 <- as.data.frame(totalS16) 
dfS16[is.na(dfS16)] <- 0
dfS16_new <- dfS16 %>% mutate(n= n.x - n.y) %>% mutate(total = total.x - total.y) %>% mutate(percentage = percentage.x - percentage.y)
df_S16 <- dfS16_new %>% select(c(class,n,total,percentage))
write_tsv(df_S16, "/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/S16class_TOTAL.tsv")
df_S161 <- df_S16 %>% mutate('Histoplasma capsulatum S16' = percentage)
df_S1611 <- df_S161 %>% select(-c(n,total,percentage))
####

library(readr)
library(plyr)
library(ggtree)


#Combine results into one large table
Total_table <- join_all(list(df_S1111,df_S1411,df_S1611,df_CB05311,df_CB05511,df_CB06211,df_CB06311,df_CB06411,df_CB06511,df_CB06611,df_CB16811,df_CB17411,df_CB18011,df_CB18611,df_CB19211,df_ES283Z11,df_ES285Z11,df_ES286Z11,df_ES288Z11,df_ES289Z11,df_ES290Z11,df_H104p0611,df_H117p1211,df_H136P0711,df_H122p10B11,df_H144p0811,df_H1517p1711,df_H256P1811,df_H316p1011,df_H327P1211,df_H343p1811,df_H388p1111,df_Hcap19VMG1511,df_HCH14311,df_HISSP1157111,df_HISSPB0582111,df_HISSPCM640811,df_HISSPFGBIK205111,df_HISSPFGBON200111,df_HISSPFGFAN205911,df_HISSPFGPERS203411,df_HISSPFGPIA205211,df_HISSPFGPIE205511,df_Histo485P2011,df_JB_08328511,df_RJ071211,df_SECH_10011,df_SECH_10111,df_SECH_10211,df_SECH_10311,df_SECH_10411,df_SECH_10511,df_SECH_10711,df_SECH_11011,df_SECH_8111,df_SECH_8211,df_SECH_8311,df_SECH_8411,df_SECH_8511,df_SECH_8611,df_SECH_8711,df_SECH_8811,df_SECH_8911,df_SECH_9011,df_SECH_9111,df_SECH_9211,df_SECH_9311,df_SECH_9411,df_SECH_9511,df_SECH_9611,df_SECH_9711,df_SECH_9811,df_SECH_9911),by="class", type="full")

write_tsv(Total_table,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Total_table.tsv")
write_tsv(FULL_Table,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/FULL_Table.tsv")

Complete_TE<- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Cleaned_Star_percentage.tsv", )
Complete_TE_short <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/Cleaned_Star_percentage_SHORT.tsv")
dComp <- as.data.frame(Complete_TE)
#convert NAs to 0
dComp[is.na(dComp)] <- 0



##Building the stacked barplot
#Found this description to follow https://zebrabi.com/guide/how-to-customize-stacked-bar-chart-in-r-ggplot2/#:~:text=Creating%20a%20basic%20stacked%20bar%20chart%20using%20ggplot2%20in%20R&text=For%20instance%2C%20we%20can%20change,bars%20using%20the%20width%20parameter

library(dplyr)
library(tidyr)
library(stringr)
library(reshape2) 
library(reshape)
library(RColorBrewer)
library(ggplot2)


#Used this tutorial to exlpain melt https://cran.r-project.org/web/packages/data.table/vignettes/datatable-reshape.html
#list the name of your species as seen in your table, include spaces
Complete_TE_meltComp<- melt(dComp, id.vars="class", 
                 measure.vars=c("Histoplasma capsulatum WGSS11","Histoplasma capsulatum WGSS14","Histoplasma capsulatum WGSS16","Histoplasma capsulatum CB053",
                                "Histoplasma capsulatum CB055","Histoplasma capsulatum CB062","Histoplasma capsulatum CB063","Histoplasma capsulatum CB064",
                                "Histoplasma capsulatum CB065","Histoplasma capsulatum CB066","Histoplasma capsulatum CB168","Histoplasma capsulatum CB174",
                                "Histoplasma capsulatum CB180","Histoplasma capsulatum CB186","Histoplasma capsulatum CB192","Histoplasma capsulatum ES2_83Z",
                                "Histoplasma capsulatum ES2_85Z","Histoplasma capsulatum ES2_86Z","Histoplasma capsulatum ES2_88Z","Histoplasma capsulatum ES2_89Z",
                                "Histoplasma capsulatum ES2_90Z","Histoplasma capsulatum 104_p_06","Histoplasma capsulatum 117_p_12","Histoplasma capsulatum 136_P_07",
                                "Histoplasma capsulatum 122_p_10B","Histoplasma capsulatum 144_p_08","Histoplasma capsulatum 1517_p_17","Histoplasma capsulatum 256_P_18",
                                "Histoplasma capsulatum 316_p_10","Histoplasma capsulatum 327_P_12","Histoplasma capsulatum 343_p_18","Histoplasma capsulatum 388_p_11",
                                "Histoplasma capsulatum 19VMG-15","Histoplasma capsulatum HCH143","Histoplasma capsulatum HISSP-11571","Histoplasma capsulatum HISSP-B05821",
                                "Histoplasma capsulatum HISSP-CM6408","Histoplasma capsulatum HISSP-FGBIK2051","Histoplasma capsulatum HISSP-FGBON2001","Histoplasma capsulatum HISSP-FGFAN2059",
                                "Histoplasma capsulatum HISSP-FGPERS2034","Histoplasma capsulatum HISSP-FGPIA2052","Histoplasma capsulatum HISSP-FGPIE2055","Histoplasma capsulatum Histo485P20",
                                "Histoplasma capsulatum JB_083285_2-Hc_083285_2","Histoplasma capsulatum 07-12-RJ","Histoplasma capsulatum SECH_100-G186A","Histoplasma capsulatum SECH_101-G184A",
                                "Histoplasma capsulatum SECH_102-505","Histoplasma capsulatum SECH_103-3_11G","Histoplasma capsulatum SECH_104-27_14","Histoplasma capsulatum SECH_105-21_14",
                                "Histoplasma capsulatum SECH_107-duboisii-B","Histoplasma capsulatum SECH_109","Histoplasma capsulatum SECH_110","Histoplasma mississippiense SECH_81-WU24",
                                "Histoplasma ohiense SECH_82-CI_4","Histoplasma mississippiense SECH_83-CI_7","Histoplasma mississippiense SECH_84-CI_19","Histoplasma mississippiense SECH_85-CI_22",
                                "Histoplasma mississippiense SECH_86-CI_24","Histoplasma mississippiense SECH_87-CI_42","Histoplasma mississippiense SECH_88-CI_43","Histoplasma mississippiense SECH_89-UCLA-531",
                                "Histoplasma mississippiense SECH_90-DOWNS","Histoplasma ohiense SECH_91-G217B","Histoplasma ohiense SECH_92-G222B","Histoplasma ohiense SECH_93-CI_6",
                                "Histoplasma ohiense SECH_94-CI_9","Histoplasma ohiense SECH_95-CI_10","Histoplasma ohiense SECH_96-CI_17","Histoplasma ohiense SECH_97-CI_18",
                                "Histoplasma ohiense SECH_98-CI_30","Histoplasma ohiense SECH_99-CI_35","Histoplasma capsulatum HISSP-CM6015","Histoplasma capsulatum HISSP-CM7256",
                                "Histoplasma capsulatum HISSP-FGAMA2041","Histoplasma capsulatum HISSP-FGFAR0189","Histoplasma capsulatum HISSP-FGFIN2028","Histoplasma capsulatum HISSP-FGJOS2044",
                                "Histoplasma capsulatum HISSP-FGPSO2043","Histoplasma capsulatum HISSP-FGTRO0285","Histoplasma capsulatum NACVFR_Histo_HC1070058_2","Histoplasma capsulatum HISSP-1014-Belem3",
                                "Histoplasma capsulatum HISSP-FGFER2036","Histoplasma capsulatum HISSP-FGLIN2055","Histoplasma capsulatum JB_01752-Hc_01752","Histoplasma capsulatum JB_021091-Hc_021091",
                                "Histoplasma capsulatum JB_031837-Hc_031837","Histoplasma capsulatum JB_042430-Hc_042430","Histoplasma capsulatum JB_062632-Hc_062632","Histoplasma capsulatum JB_062775-Hc_062775",
                                "Histoplasma capsulatum JB_073129-Hc_073129","Histoplasma capsulatum HISSP-FGMAR2044","Histoplasma capsulatum SA15","Blastomyces parvus UAMH130"))

write_tsv(Complete_TE_melt,"/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/FULL_Table.tsv")
Complete_TE_meltComp <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/FULL_Table.tsv")

#Saving results to a PDF file
pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Stacked_barplot_COMPLEX.pdf', width=15, height=20)

ggplot(Complete_TE_meltComp, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip(ylim=c()) + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of Repetative Elements per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') 

dev.off()

png(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Stacked_barplot_COMPLEX.png',width=1500, height=1900, res=100)

ggplot(Complete_TE_meltComp, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of TE per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') 

dev.off()


dd <- as.data.frame(Complete_TE_short)
#convert NAs to 0
dd[is.na(dd)] <- 0

Complete_TE_melt_short<- melt(dd, id.vars="class", 
                        measure.vars=c("Histoplasma capsulatum WGSS11","Histoplasma capsulatum WGSS14","Histoplasma capsulatum WGSS16","Histoplasma capsulatum CB053",
                                       "Histoplasma capsulatum CB055","Histoplasma capsulatum CB062","Histoplasma capsulatum CB063","Histoplasma capsulatum CB064",
                                       "Histoplasma capsulatum CB065","Histoplasma capsulatum CB066","Histoplasma capsulatum CB168","Histoplasma capsulatum CB174",
                                       "Histoplasma capsulatum CB180","Histoplasma capsulatum CB186","Histoplasma capsulatum CB192","Histoplasma capsulatum ES2_83Z",
                                       "Histoplasma capsulatum ES2_85Z","Histoplasma capsulatum ES2_86Z","Histoplasma capsulatum ES2_88Z","Histoplasma capsulatum ES2_89Z",
                                       "Histoplasma capsulatum ES2_90Z","Histoplasma capsulatum 104_p_06","Histoplasma capsulatum 117_p_12","Histoplasma capsulatum 136_P_07",
                                       "Histoplasma capsulatum 122_p_10B","Histoplasma capsulatum 144_p_08","Histoplasma capsulatum 1517_p_17","Histoplasma capsulatum 256_P_18",
                                       "Histoplasma capsulatum 316_p_10","Histoplasma capsulatum 327_P_12","Histoplasma capsulatum 343_p_18","Histoplasma capsulatum 388_p_11",
                                       "Histoplasma capsulatum 19VMG-15","Histoplasma capsulatum HCH143","Histoplasma capsulatum HISSP-11571","Histoplasma capsulatum HISSP-B05821",
                                       "Histoplasma capsulatum HISSP-CM6408","Histoplasma capsulatum HISSP-FGBIK2051","Histoplasma capsulatum HISSP-FGBON2001","Histoplasma capsulatum HISSP-FGFAN2059",
                                       "Histoplasma capsulatum HISSP-FGPERS2034","Histoplasma capsulatum HISSP-FGPIA2052","Histoplasma capsulatum HISSP-FGPIE2055","Histoplasma capsulatum Histo485P20",
                                       "Histoplasma capsulatum JB_083285_2-Hc_083285_2","Histoplasma capsulatum 07-12-RJ","Histoplasma capsulatum SECH_100-G186A","Histoplasma capsulatum SECH_101-G184A",
                                       "Histoplasma capsulatum SECH_102-505","Histoplasma capsulatum SECH_103-3_11G","Histoplasma capsulatum SECH_104-27_14","Histoplasma capsulatum SECH_105-21_14",
                                       "Histoplasma capsulatum SECH_107-duboisii-B","Histoplasma capsulatum SECH_109","Histoplasma capsulatum SECH_110","Histoplasma mississippiense SECH_81-WU24",
                                       "Histoplasma ohiense SECH_82-CI_4","Histoplasma mississippiense SECH_83-CI_7","Histoplasma mississippiense SECH_84-CI_19","Histoplasma mississippiense SECH_85-CI_22",
                                       "Histoplasma mississippiense SECH_86-CI_24","Histoplasma mississippiense SECH_87-CI_42","Histoplasma mississippiense SECH_88-CI_43","Histoplasma mississippiense SECH_89-UCLA-531",
                                       "Histoplasma mississippiense SECH_90-DOWNS","Histoplasma ohiense SECH_91-G217B","Histoplasma ohiense SECH_92-G222B","Histoplasma ohiense SECH_93-CI_6",
                                       "Histoplasma ohiense SECH_94-CI_9","Histoplasma ohiense SECH_95-CI_10","Histoplasma ohiense SECH_96-CI_17","Histoplasma ohiense SECH_97-CI_18",
                                       "Histoplasma ohiense SECH_98-CI_30","Histoplasma ohiense SECH_99-CI_35","Histoplasma capsulatum HISSP-CM6015","Histoplasma capsulatum HISSP-CM7256",
                                       "Histoplasma capsulatum HISSP-FGAMA2041","Histoplasma capsulatum HISSP-FGFAR0189","Histoplasma capsulatum HISSP-FGFIN2028","Histoplasma capsulatum HISSP-FGJOS2044",
                                       "Histoplasma capsulatum HISSP-FGPSO2043","Histoplasma capsulatum HISSP-FGTRO0285","Histoplasma capsulatum NACVFR_Histo_HC1070058_2","Histoplasma capsulatum HISSP-1014-Belem3",
                                       "Histoplasma capsulatum HISSP-FGFER2036","Histoplasma capsulatum HISSP-FGLIN2055","Histoplasma capsulatum JB_01752-Hc_01752","Histoplasma capsulatum JB_021091-Hc_021091",
                                       "Histoplasma capsulatum JB_031837-Hc_031837","Histoplasma capsulatum JB_042430-Hc_042430","Histoplasma capsulatum JB_062632-Hc_062632","Histoplasma capsulatum JB_062775-Hc_062775",
                                       "Histoplasma capsulatum JB_073129-Hc_073129","Histoplasma capsulatum HISSP-FGMAR2044","Histoplasma capsulatum SA15","Blastomyces parvus UAMH130"))

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Stacked_barplot_short.pdf', width=15, height=20)

ggplot(Complete_TE_melt_short, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of TE per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired") 

dev.off()

png(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Stacked_barplot_short.png',width=1500, height=1900, res=100)
    
ggplot(Complete_TE_melt_short, (aes(x=variable, y=value, fill=class))) + geom_bar(position="stack", stat="identity") + coord_flip() + 
  labs(title = "Histoplasma Isolates TE", y = "Percentage of TE per Genome", x = "Histoplasma Isolates") +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired") 
    
dev.off()

#+ scale_y_discrete(expand = c(0,0.8))
## Tree Building
library("ggplot2")
library("ggtree")
library("phytools")

tree_AAB<-read.tree("/nas/longleaf/home/taniak/taniak/Tree_Construction_PHYling/Trial5.AA-Blast/combined_top50.partition.ULTRA.contree.nw")
is.ultrametric(tree_AAB)

new_seq <- tree_AAB$tip.label
dd = data.frame(new_seq)
dd$Species = 0
dd$Species[grep("UAMH130", dd$new_seq)] = "Blastomyces parvus UAMH130"
dd$Species[grep("07_12-RJ", dd$new_seq)] = "Histoplasma capsulatum 07-12-RJ"
dd$Species[grep("HISSP-11571-Belem1", dd$new_seq)] = "Histoplasma capsulatum HISSP-11571"
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
dd$Species[grep("NACVFR_Histo_HC1070058", dd$new_seq)] = "Histoplasma capsulatum NACVFR_Histo_HC1070058_2"
dd$Species[grep("SECH_103", dd$new_seq)] = "Histoplasma capsulatum SECH_103-3_11G"
dd$Species[grep("HISSP-FGPERS2034", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGPERS2034"
dd$Species[grep("HISSP-B05821", dd$new_seq)] = "Histoplasma capsulatum HISSP-B05821"
dd$Species[grep("JB_031837-Hc_031837", dd$new_seq)] = "Histoplasma capsulatum JB_031837-Hc_031837"
dd$Species[grep("JB_042430-Hc_042430", dd$new_seq)] = "Histoplasma capsulatum JB_042430-Hc_042430"
dd$Species[grep("JB_062775-Hc_062775", dd$new_seq)] = "Histoplasma capsulatum JB_062775-Hc_062775"
dd$Species[grep("JB_083285_2-Hc_083285_2", dd$new_seq)] = "Histoplasma capsulatum JB_083285_2-Hc_083285_2"
dd$Species[grep("CB053", dd$new_seq)] = "Histoplasma capsulatum CB053" 
dd$Species[grep("CB055", dd$new_seq)] = "Histoplasma capsulatum CB055"
dd$Species[grep("G184A", dd$new_seq)] = "Histoplasma capsulatum SECH_101-G184A"
dd$Species[grep("HISSP-CM6408", dd$new_seq)] = "Histoplasma capsulatum HISSP-CM6408"
dd$Species[grep("HISSP-1014-Belem3", dd$new_seq)] = "Histoplasma capsulatum HISSP-1014-Belem3"
dd$Species[grep("ES2_83Z", dd$new_seq)] = "Histoplasma capsulatum ES2_83Z"
dd$Species[grep("ES2_90Z", dd$new_seq)] = "Histoplasma capsulatum ES2_90Z"
dd$Species[grep("ES2_86Z", dd$new_seq)] = "Histoplasma capsulatum ES2_86Z"
dd$Species[grep("ES2_88Z", dd$new_seq)] = "Histoplasma capsulatum ES2_88Z"
dd$Species[grep("HCH143", dd$new_seq)] = "Histoplasma capsulatum HCH143"
dd$Species[grep("duboisii", dd$new_seq)] = "Histoplasma capsulatum SECH_107-duboisii-B" 
dd$Species[grep("ES2_85Z", dd$new_seq)] = "Histoplasma capsulatum ES2_85Z"
dd$Species[grep("HISSP-CM7256", dd$new_seq)] = "Histoplasma capsulatum HISSP-CM7256"
dd$Species[grep("ES2_89Z", dd$new_seq)] = "Histoplasma capsulatum ES2_89Z"
dd$Species[grep("SA15", dd$new_seq)] = "Histoplasma capsulatum SA15"
dd$Species[grep("SECH_98", dd$new_seq)] = "Histoplasma ohiense SECH_98-CI_30"
dd$Species[grep("SECH_110", dd$new_seq)] = "Histoplasma capsulatum SECH_110"
dd$Species[grep("G217B", dd$new_seq)] = "Histoplasma ohiense SECH_91-G217B"
dd$Species[grep("SECH_92", dd$new_seq)] = "Histoplasma ohiense SECH_92-G222B"
dd$Species[grep("Hc1986", dd$new_seq)] = "Histoplasma ohiense Hc1986"
dd$Species[grep("SECH_82", dd$new_seq)] = "Histoplasma ohiense SECH_82-CI_4"
dd$Species[grep("SECH_95", dd$new_seq)] = "Histoplasma ohiense SECH_95-CI_10"
dd$Species[grep("SECH_94", dd$new_seq)] = "Histoplasma ohiense SECH_94-CI_9"
dd$Species[grep("SECH_96", dd$new_seq)] = "Histoplasma ohiense SECH_96-CI_17"
dd$Species[grep("SECH_93", dd$new_seq)] = "Histoplasma ohiense SECH_93-CI_6"
dd$Species[grep("SECH_97", dd$new_seq)] = "Histoplasma ohiense SECH_97-CI_18"
dd$Species[grep("SECH_99", dd$new_seq)] = "Histoplasma ohiense SECH_99-CI_35"
dd$Species[grep("SECH_104", dd$new_seq)] = "Histoplasma capsulatum SECH_104-27_14" 
dd$Species[grep("CB174", dd$new_seq)] = "Histoplasma capsulatum CB174"
dd$Species[grep("SECH_109", dd$new_seq)] = "Histoplasma capsulatum SECH_109"
dd$Species[grep("CB062", dd$new_seq)] = "Histoplasma capsulatum CB062"
dd$Species[grep("CB065", dd$new_seq)] = "Histoplasma capsulatum CB065"
dd$Species[grep("CB192", dd$new_seq)] = "Histoplasma capsulatum CB192"
dd$Species[grep("CB064", dd$new_seq)] = "Histoplasma capsulatum CB064"
dd$Species[grep("CB180", dd$new_seq)] = "Histoplasma capsulatum CB180"
dd$Species[grep("CB063", dd$new_seq)] = "Histoplasma capsulatum CB063"
dd$Species[grep("CB186", dd$new_seq)] = "Histoplasma capsulatum CB186" 
dd$Species[grep("SECH_81", dd$new_seq)] = "Histoplasma mississippiense SECH_81-WU24"
dd$Species[grep("CB066", dd$new_seq)] = "Histoplasma capsulatum CB066"
dd$Species[grep("SECH_88", dd$new_seq)] = "Histoplasma mississippiense SECH_88-CI_43"
dd$Species[grep("SECH_86", dd$new_seq)] = "Histoplasma mississippiense SECH_86-CI_24"
dd$Species[grep("CB168", dd$new_seq)] = "Histoplasma capsulatum CB168"
dd$Species[grep("SECH_85", dd$new_seq)] = "Histoplasma mississippiense SECH_85-CI_22"
dd$Species[grep("SECH_83", dd$new_seq)] = "Histoplasma mississippiense SECH_83-CI_7" 
dd$Species[grep("SECH_87", dd$new_seq)] = "Histoplasma mississippiense SECH_87-CI_42"
dd$Species[grep("SECH_90", dd$new_seq)] = "Histoplasma mississippiense SECH_90-DOWNS"
dd$Species[grep("SECH_102", dd$new_seq)] = "Histoplasma capsulatum SECH_102-505"
dd$Species[grep("SECH_100", dd$new_seq)] = "Histoplasma capsulatum SECH_100-G186A"
dd$Species[grep("SECH_84", dd$new_seq)] = "Histoplasma mississippiense SECH_84-CI_19"
dd$Species[grep("SECH_89", dd$new_seq)] = "Histoplasma mississippiense SECH_89-UCLA-531"
dd$Species[grep("JB_01752-Hc_01752", dd$new_seq)] = "Histoplasma capsulatum JB_01752-Hc_01752"
dd$Species[grep("104_p_06", dd$new_seq)] = "Histoplasma capsulatum 104_p_06"
dd$Species[grep("1517_p_17", dd$new_seq)] = "Histoplasma capsulatum 1517_p_17"
dd$Species[grep("117_p_12", dd$new_seq)] = "Histoplasma capsulatum 117_p_12"
dd$Species[grep("388_p_11", dd$new_seq)] = "Histoplasma capsulatum 388_p_11" 
dd$Species[grep("144_p_08", dd$new_seq)] = "Histoplasma capsulatum 144_p_08"
dd$Species[grep("327_P_12", dd$new_seq)] = "Histoplasma capsulatum 327_P_12"
dd$Species[grep("316_p_10", dd$new_seq)] = "Histoplasma capsulatum 316_p_10"
dd$Species[grep("WGS_S14", dd$new_seq)] = "Histoplasma capsulatum WGSS14"
dd$Species[grep("Histo-485P20", dd$new_seq)] = "Histoplasma capsulatum Histo485P20"
dd$Species[grep("136_P_07", dd$new_seq)] = "Histoplasma capsulatum 136_P_07"
dd$Species[grep("343_p_18", dd$new_seq)] = "Histoplasma capsulatum 343_p_18"
dd$Species[grep("WGS_S11", dd$new_seq)] = "Histoplasma capsulatum WGSS11"
dd$Species[grep("WGS_S16", dd$new_seq)] = "Histoplasma capsulatum WGSS16"
dd$Species[grep("122_p_10_B", dd$new_seq)] = "Histoplasma capsulatum 122_p_10B"
dd$Species[grep("256_P_18", dd$new_seq)] = "Histoplasma capsulatum 256_P_18" 
dd$Species[grep("JB_021091-Hc_021091", dd$new_seq)] = "Histoplasma capsulatum JB_021091-Hc_021091"
dd$Species[grep("JB_073129-Hc_073129", dd$new_seq)] = "Histoplasma capsulatum JB_073129-Hc_073129"
dd$Species[grep("HISSP-FGFAR0189", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGFAR0189"
dd$Species[grep("HISSP-FGFIN2028", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGFIN2028"
dd$Species[grep("HISSP-FGPSO2043", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGPSO2043"
dd$Species[grep("HISSP-FGMAR2044", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGMAR2044"
dd$Species[grep("HISSP-FGTRO0285", dd$new_seq)] = "Histoplasma capsulatum HISSP-FGTRO0285"
dd$Species[grep("SECH_105", dd$new_seq)] = "Histoplasma capsulatum SECH_105-21_14"
dd$Species[grep("19VMG-15", dd$new_seq)] = "Histoplasma capsulatum 19VMG-15"
dd$Species[grep("SECH_101", dd$new_seq)] = "Histoplasma capsulatum SECH_101-G184A"
tree_AAB$tip.label <- dd$Species


#Sanity check!
colnames(Complete_TE_short)[which(colnames(Complete_TE_short)  %in% tree_AAB$tip.label==F)]
p_AABBB1 <- ggtree(tree_AAB, layout = "rectangular") + geom_tiplab(size = 3, fontface = 3, align=TRUE) +
  theme_tree2() + theme(legend.position = "none") 
print(p_AABBB1)

colnames(Complete_TE_melt_short)
Complete_TE_melt_short2 <- Complete_TE_melt_short[,c(2,1,3)]

#pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Tree_fan_Stacked_barplot_AAB-new.pdf', width=15, height=15)

p_AABBB2 <- p_AABBB1 + geom_facet(panel = "TE Abundance", data = Complete_TE_melt_short2, geom = geom_col, 
                                  aes(x = value, color = class, 
                                      fill = class), orientation = 'y', width = .6) + geom_nodelab(size = 3, na.rm = TRUE, nudge_x = 0.05) +
  theme_tree2(legend.position=c(.95, .25)) + xlim_expand(c(0,0.8), 'TE Abundance') + xlim_expand(c(0, 4), 'Tree') + geom_cladelabel(node=140, label="NAm1", color="red", offset = 1.3) +
  geom_cladelabel(node=160, label="NAm2", color="#808000", offset = 1.25) + geom_cladelabel(node=172, label="India", color="blue", offset = 1.25) + geom_cladelabel(node=113, label="Africa", color="#FFA500", offset = 1.25) +
  geom_cladelabel(node=105, label="LAmA 1", color="purple", offset = 1.28) + geom_cladelabel(node=121, label="LAmA 2", color="purple", offset = 1.28) + geom_cladelabel(node=187, label="LAmA O", color="purple", offset = 1.25) +
  geom_cladelabel(node=191, label="Africa O", color="#FFA500", offset = 1.25) + geom_cladelabel(node=137, label="LAmA 3", color="purple", offset = 1.28) + geom_cladelabel(node=1, label="Root", color="green", offset = 1.3)

#dev.off()

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Tree_Stacked_barplotAAB_Star1.pdf', width=15, height=15)

facet_widths(p_AABBB2, widths = c(1, 0.40))

dev.off()

png(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Tree_Stacked_barplotAAB_Star1.png', width=1700, height=1500, res=100)

facet_widths(p_AABBB2, widths = c(1, 0.40))

dev.off()



#RUN significance tests on tree_AAB

library('phytools')
colnames(Complete_TE_short) = gsub("[.]", "_", colnames(Complete_TE_short))
tree_AAB$tip.label = gsub("[-]", "_", tree_AAB$tip.label)
all(colnames(Complete_TE_short)  %in% tree_AAB$tip.label)
colnames(Complete_TE_short)[which(colnames(Complete_TE_short)  %in% tree_AAB$tip.label==F)]
#Just first column doesn't agree which is class

#Test TE for phylogenetic signal - but phylosig requires named vector for data
row_DNA = which(Complete_TE_short$class=="DNA TE")
data.vec_DNA = as.numeric(Complete_TE_short[row_DNA, ]); names(data.vec_DNA) = colnames(Complete_TE_short)

# remove the first column, which contains no data
data.vec_DNA=data.vec_DNA[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.DNA = phylosig(tree=tree_AAB, x=data.vec_DNA, method="lambda", test=TRUE)
sig.lam.DNA
plot.phylosig(sig.lam.DNA)


#Phylogenetic signal lambda : 0.852516 
#logL(lambda) : 237.552 
#LR(lambda=0) : 4.49699 
#P-value (based on LR test) : 0.0339546 

sig.k.DNA = (phylosig(tree=tree_AAB, x=data.vec_DNA, method="K", test=TRUE))
sig.k.DNA
plot.phylosig(sig.k.DNA)

#Phylogenetic signal K : 0.50602 
#P-value (based on 1000 randomizations) : 0.144

### RNA-TE
row_RNA = which(Complete_TE_short$class=="RNA TE")
data.vec_RNA = as.numeric(Complete_TE_short[row_RNA, ]); names(data.vec_RNA) = colnames(Complete_TE_short)

# remove the first column, which contains no data
data.vec_RNA=data.vec_RNA[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.RNA = phylosig(tree=tree_AAB, x=data.vec_RNA, method="lambda", test=TRUE)
sig.lam.RNA
plot.phylosig(sig.lam.RNA)

#Phylogenetic signal lambda : 0.829077 
#logL(lambda) : 128.61 
#LR(lambda=0) : 59.7364 
#P-value (based on LR test) : 1.0845e-14 

sig.k.RNA = (phylosig(tree=tree_AAB, x=data.vec_RNA, method="K", test=TRUE))
sig.k.RNA
plot.phylosig(sig.k.RNA)

#Phylogenetic signal K : 0.926126 
#P-value (based on 1000 randomizations) : 0.001 

### LC-TE
row_LC = which(Complete_TE_short$class=="Low_complexity")
data.vec_LC = as.numeric(Complete_TE_short[row_LC, ]); names(data.vec_LC) = colnames(Complete_TE_short)

# remove the first column, which contains no data
data.vec_LC=data.vec_LC[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.LC = phylosig(tree=tree_AAB, x=data.vec_LC, method="lambda", test=TRUE)
sig.lam.LC
plot.phylosig(sig.lam.LC)


#Phylogenetic signal lambda : 0.805051 
#logL(lambda) : 542.184 
#LR(lambda=0) : 63.0754 
#P-value (based on LR test) : 1.98939e-15 

sig.k.LC = (phylosig(tree=tree_AAB, x=data.vec_LC, method="K", test=TRUE))
sig.k.LC
plot.phylosig(sig.k.LC)

#Phylogenetic signal K : 0.958456 
#P-value (based on 1000 randomizations) : 0.001 

### MITE
row_MITE = which(Complete_TE_short$class=="MITE")
data.vec_MITE = as.numeric(Complete_TE_short[row_MITE, ]); names(data.vec_MITE) = colnames(Complete_TE_short)

# remove the first column, which contains no data
data.vec_MITE=data.vec_MITE[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.MITE = phylosig(tree=tree_AAB, x=data.vec_MITE, method="lambda", test=TRUE)
sig.lam.MITE
plot.phylosig(sig.lam.MITE)

#Phylogenetic signal lambda : 0.999927 
#logL(lambda) : 453.416 
#LR(lambda=0) : 16.1056 
#P-value (based on LR test) : 5.9906e-05 

sig.k.MITE = (phylosig(tree=tree_AAB, x=data.vec_MITE, method="K", test=TRUE))
sig.k.MITE
plot.phylosig(sig.k.MITE)

#Phylogenetic signal K : 0.59035 
#P-value (based on 1000 randomizations) : 0.171 

### Satellite
row_SAT = which(Complete_TE_short$class=="Satellite")
data.vec_SAT = as.numeric(Complete_TE_short[row_SAT, ]); names(data.vec_SAT) = colnames(Complete_TE_short)

# remove the first column, which contains no data
data.vec_SAT=data.vec_SAT[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.SAT = phylosig(tree=tree_AAB, x=data.vec_SAT, method="lambda", test=TRUE)
sig.lam.SAT
plot.phylosig(sig.lam.SAT)

#Phylogenetic signal lambda : 0.999927 
#logL(lambda) : 1145.24 
#LR(lambda=0) : 15.1406 
#P-value (based on LR test) : 9.97952e-05 

sig.k.SAT = (phylosig(tree=tree_AAB, x=data.vec_SAT, method="K", test=TRUE))
sig.k.SAT
plot.phylosig(sig.k.SAT)

#Phylogenetic signal K : 0.579515 
#P-value (based on 1000 randomizations) : 0.208 

### Simple Repeat
row_SR = which(Complete_TE_short$class=="Simple_repeat")
data.vec_SR = as.numeric(Complete_TE_short[row_SR, ]); names(data.vec_SR) = colnames(Complete_TE_short)

# remove the first column, which contains no data
data.vec_SR=data.vec_SR[-1]

# we'll use the test=TRUE argument so that it will tell us whether the signal is 'significant'
sig.lam.SR = phylosig(tree=tree_AAB, x=data.vec_SR, method="lambda", test=TRUE)
sig.lam.SR
plot.phylosig(sig.lam.SR)

#Phylogenetic signal lambda : 0.887888 
#logL(lambda) : 420.138 
#LR(lambda=0) : 79.9757 
#P-value (based on LR test) : 3.79044e-19 

sig.k.SR = (phylosig(tree=tree_AAB, x=data.vec_SR, method="K", test=TRUE))
sig.k.SR
plot.phylosig(sig.k.SR)

#Phylogenetic signal K : 1.18819 
#P-value (based on 1000 randomizations) : 0.001 

## Starships
row_Star = which(Complete_TE_short$class=="Starship")
data.vec_Star = as.numeric(Complete_TE_short[row_Star, ]); names(data.vec_Star) = colnames(Complete_TE_short)

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


#Comparing Genome Length RE and Starships
restarlength <- read_tsv("/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/re_star_length.tsv")

#RDD <- as.data.frame(restarlength)
restarlength_melt<- melt(restarlength, id.vars="class", 
                              measure.vars=c("Histoplasma capsulatum WGSS11","Histoplasma capsulatum WGSS14","Histoplasma capsulatum WGSS16","Histoplasma capsulatum CB053",
                                             "Histoplasma capsulatum CB055","Histoplasma capsulatum CB062","Histoplasma capsulatum CB063","Histoplasma capsulatum CB064",
                                             "Histoplasma capsulatum CB065","Histoplasma capsulatum CB066","Histoplasma capsulatum CB168","Histoplasma capsulatum CB174",
                                             "Histoplasma capsulatum CB180","Histoplasma capsulatum CB186","Histoplasma capsulatum CB192","Histoplasma capsulatum ES2_83Z",
                                             "Histoplasma capsulatum ES2_85Z","Histoplasma capsulatum ES2_86Z","Histoplasma capsulatum ES2_88Z","Histoplasma capsulatum ES2_89Z",
                                             "Histoplasma capsulatum ES2_90Z","Histoplasma capsulatum 104_p_06","Histoplasma capsulatum 117_p_12","Histoplasma capsulatum 136_P_07",
                                             "Histoplasma capsulatum 122_p_10B","Histoplasma capsulatum 144_p_08","Histoplasma capsulatum 1517_p_17","Histoplasma capsulatum 256_P_18",
                                             "Histoplasma capsulatum 316_p_10","Histoplasma capsulatum 327_P_12","Histoplasma capsulatum 343_p_18","Histoplasma capsulatum 388_p_11",
                                             "Histoplasma capsulatum 19VMG-15","Histoplasma capsulatum HCH143","Histoplasma capsulatum HISSP-11571","Histoplasma capsulatum HISSP-B05821",
                                             "Histoplasma capsulatum HISSP-CM6408","Histoplasma capsulatum HISSP-FGBIK2051","Histoplasma capsulatum HISSP-FGBON2001","Histoplasma capsulatum HISSP-FGFAN2059",
                                             "Histoplasma capsulatum HISSP-FGPERS2034","Histoplasma capsulatum HISSP-FGPIA2052","Histoplasma capsulatum HISSP-FGPIE2055","Histoplasma capsulatum Histo485P20",
                                             "Histoplasma capsulatum JB_083285_2-Hc_083285_2","Histoplasma capsulatum 07-12-RJ","Histoplasma capsulatum SECH_100-G186A","Histoplasma capsulatum SECH_101-G184A",
                                             "Histoplasma capsulatum SECH_102-505","Histoplasma capsulatum SECH_103-3_11G","Histoplasma capsulatum SECH_104-27_14","Histoplasma capsulatum SECH_105-21_14",
                                             "Histoplasma capsulatum SECH_107-duboisii-B","Histoplasma capsulatum SECH_109","Histoplasma capsulatum SECH_110","Histoplasma mississippiense SECH_81-WU24",
                                             "Histoplasma ohiense SECH_82-CI_4","Histoplasma mississippiense SECH_83-CI_7","Histoplasma mississippiense SECH_84-CI_19","Histoplasma mississippiense SECH_85-CI_22",
                                             "Histoplasma mississippiense SECH_86-CI_24","Histoplasma mississippiense SECH_87-CI_42","Histoplasma mississippiense SECH_88-CI_43","Histoplasma mississippiense SECH_89-UCLA-531",
                                             "Histoplasma mississippiense SECH_90-DOWNS","Histoplasma ohiense SECH_91-G217B","Histoplasma ohiense SECH_92-G222B","Histoplasma ohiense SECH_93-CI_6",
                                             "Histoplasma ohiense SECH_94-CI_9","Histoplasma ohiense SECH_95-CI_10","Histoplasma ohiense SECH_96-CI_17","Histoplasma ohiense SECH_97-CI_18",
                                             "Histoplasma ohiense SECH_98-CI_30","Histoplasma ohiense SECH_99-CI_35","Histoplasma capsulatum HISSP-CM6015","Histoplasma capsulatum HISSP-CM7256",
                                             "Histoplasma capsulatum HISSP-FGAMA2041","Histoplasma capsulatum HISSP-FGFAR0189","Histoplasma capsulatum HISSP-FGFIN2028","Histoplasma capsulatum HISSP-FGJOS2044",
                                             "Histoplasma capsulatum HISSP-FGPSO2043","Histoplasma capsulatum HISSP-FGTRO0285","Histoplasma capsulatum NACVFR_Histo_HC1070058_2","Histoplasma capsulatum HISSP-1014-Belem3",
                                             "Histoplasma capsulatum HISSP-FGFER2036","Histoplasma capsulatum HISSP-FGLIN2055","Histoplasma capsulatum JB_01752-Hc_01752","Histoplasma capsulatum JB_021091-Hc_021091",
                                             "Histoplasma capsulatum JB_031837-Hc_031837","Histoplasma capsulatum JB_042430-Hc_042430","Histoplasma capsulatum JB_062632-Hc_062632","Histoplasma capsulatum JB_062775-Hc_062775",
                                             "Histoplasma capsulatum JB_073129-Hc_073129","Histoplasma capsulatum HISSP-FGMAR2044","Histoplasma capsulatum SA15"))


pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Starships/ALL_Attempt2/new_RM/Star_RE_GL.pdf', width=15, height=20)
ggplot(restarlength_melt, (aes(x=variable, y=value, fill=class)))+ geom_bar(position="stack", stat="identity") + 
  coord_flip() + labs(title = "Histoplasma Repetative Elements vs Genome Length", y = "Genome Assembly (Mb)", x = "Histoplasma Isolates") + 
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = "Paired")
dev.off()

pdf(file='/nas/longleaf/home/taniak/taniak/TE_Analysis/Summ_GL_UNKNOWN_RE.pdf.pdf', width=15, height=20)
Stack_GL_TE_UNKNOWN + geom_facet(panel = "TE Abundance", data = Stack_GL_TE_UNKNOWN, geom = geom_col, aes(x = value, color = class, fill = class), orientation = 'y', width = .6) 
dev.off()