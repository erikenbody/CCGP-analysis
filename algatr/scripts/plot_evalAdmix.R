library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
qfile <- args[1]
in_file <- args[2]
out_file <- args[3]
K <- as.numeric(args[4])
vis <- args[5]
prefix <- args[6]

source(vis)

r <- as.matrix(read.table(in_file))

q<-read.table(qfile,stringsAsFactors=T)

q2 <- q %>% mutate(pop = gsub("V","",names(q)[max.col(q)]))


palette(c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F",
          "#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#1B9E77","#999999"))

# order according to population and plot the ADMIXTURE reults
ord<-orderInds(pop = as.vector(q2$pop), q = q)

#make barplot
#plotAdmix(q,ord=ord,pop=q2$pop)
png(out_file, width = 8, height = 6, units = "in", res = 350)

plotCorRes(
  cor_mat = r, 
  pop = as.vector(q2$pop), 
  ord = ord, 
  title = paste("evalAdmix K=",K," -", prefix), 
  max_z = 0.1, 
  min_z = -0.1
)
dev.off()
