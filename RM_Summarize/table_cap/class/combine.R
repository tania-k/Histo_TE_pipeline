library(dplyr)
library(readr)

#df <- list.files(path="/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/class", full.names = TRUE) %>% 
#  lapply(read_tsv) %>% bind_rows 


library(data.table)

 fileNames = dir("/nas/longleaf/home/taniak/taniak/TE_Analysis/RM_Summarize/table_cap/class", pattern = "*.tsv", full.names = TRUE)
 listOfTables = lapply(fileNames, fread)

 # convert all tables to a big df
 df = rbindlist(listOfTables, fill=TRUE)

 # filter only the rows you want:
 newDf = df[n %chin% c("DNA/DTC", "DNA/DTA", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron", "LTR/Copia", "LTR/Gypsy", "LTR/TRIM", "LTR/unknown", "Low_complexity", "MITE/DTH", "MITE/DTC", "MITE/DTM", "MITE/DTT", "Simple_repeat", "TIR/PiggyBac", "TIR/Tc1_Mariner", "LINE/unknown", "DNAnona/Helitron", "DNAauto/Helitron", "LTR/Copia_1", "LTR/Copia_2", "LTR/Copia_4", "LTR/Gypsy_1", "LTR/Gypsy_2", "LTR/Gypsy_3", "TIR/PiggyBac"), ]
