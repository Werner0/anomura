library(mygene)
library(data.table)
library(plyr)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/mygene/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/mygene/")}

expansiongenes <- fread("expansion.csv")
expansiongenesnames <- expansiongenes$Gene

mygenes <- c("KI")
mygenes <- expansiongenesnames
res <- queryMany(mygenes, scopes='symbol', fields=c('entrez', 'go'),species='zebrafish')
res <- as.data.table(res)
restable <- data.table()
for (i in mygenes) {
  CCterm <- toString(res$go.CC[[which(i == mygenes)]]$term)
  CCGO <- toString(res$go.CC[[which(i == mygenes)]]$id)
  MFterm <- toString(res$go.MF[[which(i == mygenes)]]$term)
  MFGO <- toString(res$go.MF[[which(i == mygenes)]]$id)
  BPterm <- toString(res$go.BP[[which(i == mygenes)]]$term)
  BPGO <- toString(res$go.BP[[which(i == mygenes)]]$id)
  currenttable <- data.table(gene = i,CC_term = CCterm, CCgo = CCGO, MFterm = MFterm, MFgo = MFGO, BPterm = BPterm, BPgo = BPGO)
  restable <- rbind(restable,currenttable)
}

#PArt 2
gosCC <- c()
eviCC <- c()
termCC <- c()
queryCC <- c()
scoreCC <- c()
gosMF <- c()
eviMF <- c()
termMF <- c()
queryMF <- c()
scoreMF <- c()
gosBP <- c()
eviBP <- c()
termBP <- c()
queryBP <- c()
scoreBP <- c()
for (i in mygenes) {
  gosCC <- c(gosCC,res$go.CC[[which(i == mygenes)]]$id)
  eviCC <- c(eviCC,res$go.CC[[which(i == mygenes)]]$evidence)
  termCC <- c(termCC,res$go.CC[[which(i == mygenes)]]$term)
  queryCC <- c(queryCC,rep(res$query[which(i == mygenes)],length(res$go.CC[[which(i == mygenes)]]$id)))
  scoreCC <- c(scoreCC,rep(res$X_score[which(i == mygenes)],length(res$go.CC[[which(i == mygenes)]]$id)))
}
bigtableCC <- data.table(SYMBOL = queryCC, SCORE = scoreCC, 
                       GO = gosCC, EVIDENCE = eviCC, DESC = termCC, TYPE = "CC")
for (i in mygenes) {
  gosMF <- c(gosMF,res$go.MF[[which(i == mygenes)]]$id)
  eviMF <- c(eviMF,res$go.MF[[which(i == mygenes)]]$evidence)
  termMF <- c(termMF,res$go.MF[[which(i == mygenes)]]$term)
  queryMF <- c(queryMF,rep(res$query[which(i == mygenes)],length(res$go.MF[[which(i == mygenes)]]$id)))
  scoreMF <- c(scoreMF,rep(res$X_score[which(i == mygenes)],length(res$go.MF[[which(i == mygenes)]]$id)))
}
bigtableMF <- data.table(SYMBOL = queryMF, SCORE = scoreMF, 
                         GO = gosMF, EVIDENCE = eviMF, DESC = termMF, TYPE = "MF")
for (i in mygenes) {
  gosBP <- c(gosBP,res$go.BP[[which(i == mygenes)]]$id)
  eviBP <- c(eviBP,res$go.BP[[which(i == mygenes)]]$evidence)
  termBP <- c(termBP,res$go.BP[[which(i == mygenes)]]$term)
  queryBP <- c(queryBP,rep(res$query[which(i == mygenes)],length(res$go.BP[[which(i == mygenes)]]$id)))
  scoreBP <- c(scoreBP,rep(res$X_score[which(i == mygenes)],length(res$go.BP[[which(i == mygenes)]]$id)))
}
bigtableBP <- data.table(SYMBOL = queryBP, SCORE = scoreBP, 
                         GO = gosBP, EVIDENCE = eviBP, DESC = termBP, TYPE = "BP")


#ONE BY ONE
# FLYBASE <- rbind(bigtableBP,bigtableCC,bigtableMF)
# FLYBASE[,"DB"] <- "FLYBASE"
# 
# ENTREZ_fruitfly <- rbind(bigtableBP,bigtableCC,bigtableMF)
# ENTREZ_fruitfly[,"DB"] <- "ENTREZ_fruitfly"
# 
# ENTREZ_human <- rbind(bigtableBP,bigtableCC,bigtableMF)
# ENTREZ_human[,"DB"] <- "ENTREZ_human"
# 
# ENTREZ_zebrafish <- rbind(bigtableBP,bigtableCC,bigtableMF)
# ENTREZ_zebrafish[,"DB"] <- "ENTREZ_zebrafish"
# 
# ENTREZ_nematode <- rbind(bigtableBP,bigtableCC,bigtableMF)
# ENTREZ_nematode[,"DB"] <- "ENTREZ_nematode"
# 
# ALLDB <- rbind(FLYBASE,ENTREZ_fruitfly,ENTREZ_human,ENTREZ_zebrafish,ENTREZ_nematode)
#fwrite(ALLDB, file = "FUNCTIONALANNO.csv", sep = ",")
