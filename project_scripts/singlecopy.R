library(data.table)
library(plyr)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/single_copies/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/single_copies/")}

OGs <- fread("OGs.txt", header = F)
ids <- fread("identifiers.txt",header = F)
reference <- fread("alldata.txt")
combo <- data.table(OGs, ids)
colnames(combo) <- c("OG","query_name")
mynames <- c("amphipod", "bluekingcrab", "coconutcrab", "isopod", "redkingcrab", "lobster", "marbledcrayfish", "swimmingcrab", "whiteshrimp")
combo$organism <- rep(mynames,40)
mytable <- merge(combo,reference,by=c("query_name","organism"), all = T)
mytable <- mytable[mytable$OG!=0,]
mytable <- mytable[,c(3,1,2,5,4)]
mytable$OG <- gsub("\\..*","",mytable$OG)
colnames(mytable) <- c("OG","fasta_header", "organism", "taxon", "KOG")
art <- mytable[mytable$taxon=="Arthropoda",]
