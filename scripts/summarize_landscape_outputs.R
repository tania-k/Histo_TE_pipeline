library(tidyverse)
library(readr)
library(dplyr)
library(ggsci)
library(gridExtra)
library(ggplot2)
library(hrbrthemes)
setwd("/nas/longleaf/home/taniak/taniak/TE_Analysis")

args <- commandArgs()
cat(args, sep = "\n")
args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

divsumtbl <- args[1]
genome_size <- strtoi(args[2])
divsumtbl

outcsv=gsub(pattern = "\\.tbl$", ".kimuraTE_summary.csv", divsumtbl)
outtsv=gsub(pattern = "\\.tbl$", ".kimutaTE_summary.tsv", divsumtbl)
species = gsub(pattern ="\\.tbl$","",basename(divsumtbl))

Org1 <- read_table(divsumtbl,skip=3,col_names = F)

colnames(KimuraDistance) <-c("score","div.","del.","ins.","query","beginq","endq","leftq","strand","family","class","beginr","endr","leftr","ID","overlap")
KimuraDistance <- KimuraDistance %>% mutate(length=endq-beginq +1) 
KimuraDistance1 <- subset(KimuraDistance,class != "Low_complexity")
KimuraDistance2 <- subset(KimuraDistance1,class != "Simple_repeat")
write_csv(outcsv, path="RM_Summarize/table/divsumtbl.tsv")
#KimuraDistanceclass<- KimuraDistance2 %>% group_by(class) %>%
#  summarise(n=n_distinct(ID),
#            total=sum(length),
#            percentage=total/genome_size)
#write_tsv(outtsv, path="RM_Summarize/table/divsumtbl_sum.tsv)
