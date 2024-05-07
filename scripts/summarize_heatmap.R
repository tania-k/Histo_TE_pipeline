#heatmap
library(RColorBrewer)
library(gplots)
library(pheatmap)
DNATE <- read_tsv("tables/TE composition - 7.7.23.tsv",col_names = T)
as.data.frame(DNATE)
rownames(DNATE)<-DNATE$class
DF<- DNATE %>% select(-class)
heatmap(DF)
DF<-sapply(DF,as.numeric)
dev.off()
par(mfrow=c(2,2))
dpi = 100
fig.width = 5
fig.height = 5
p1<- heatmap.2(DF, scale = "column", col = bluered(100), 
          trace = "none", density.info = "none",Colv = NA, Rowv = NA,labRow = x,margins = c(15, 15))
ggsave("heatmap.pdf",p1)
pheatmap::pheatmap(DF,cutree_rows = 4)

x <- DNATE$class
y <- c("Mutator","CMC","hAT","TcMar","PIF_Harbinger","PiggyBac","Merlin")
z<-DNATE[2:18,2:8]
z <- sapply(z, as.numeric)
stopifnot(all(is.numeric(z)))
heatmap(z,scale = "row",Colv = NA, Rowv = NA,)


data <- expand.grid(X=x, Y=y)
data$Z <- z
data %>% ggplot(aes(X, Y, fill= Z)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdPu") +
  theme_ipsum()


#7.19 heatmap
library(RColorBrewer)
library(pheatmap)
library(NatParksPalettes)
library(wesanderson)
setwd("~/bigdata/TE_composition-EDTA/")
DNATE <- read.table("~/bigdata/TE_composition-EDTA/tables/TE composition - Sheet12.csv",sep=",",row.names = 1,header=T)
as.data.frame(DNATE)
DNATE<-sapply(DNATE,as.numeric)

p1<-pheatmap(DNATE,display_numbers = T,cluster_rows = F,cluster_cols = F,main = "% of genome occupied by different DNA transposons",fontsize_number = 12)


ggsave("plots/DNAheatmap1.pdf",p1)
pheatmap::pheatmap(DF,cutree_rows = 4)

x <- DNATE$class
y <- c("Mutator","CMC","hAT","TcMar","PIF_Harbinger","PiggyBac","Merlin")
z<-DNATE[2:18,2:8]
z <- sapply(z, as.numeric)
stopifnot(all(is.numeric(z)))
heatmap(z,scale = "row",Colv = NA, Rowv = NA,)


data <- expand.grid(X=x, Y=y)
data$Z <- z
data %>% ggplot(aes(X, Y, fill= Z)) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdPu") +
  theme_ipsum()
