#SNIPPETS (USE BEFORE SCRIPT)
#New gene combination in second section of this script
#awk -F\# '$1!="" { print $0 ;} ' whiteshrimp.emapper.annotations > whiteshrimp.emapper_nocomments.annotations

library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/eggnog_perRead/Whiteshrimp/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/eggnog_perRead/Whiteshrimp/")}

library(plyr)
library(tibble)

whiteshrimp <- fread("whiteshrimp.emapper_nocomments.annotations")
whiteshrimpHP <- fread("HPReadsPerGene.out.tab", skip = 4)
whiteshrimpHP2 <- fread("HP2ReadsPerGene.out.tab", skip = 4)
whiteshrimpES <- fread("ESReadsPerGene.out.tab", skip = 4)
whiteshrimpGILL <- fread("GReadsPerGene.out.tab", skip = 4)
whiteshrimpGILL2 <- fread("G2ReadsPerGene.out.tab", skip = 4)
whiteshrimpMS <- fread("MReadsPerGene.out.tab", skip = 4)
whiteshrimpMS2 <- fread("M2ReadsPerGene.out.tab", skip = 4)
whiteshrimpO <- fread("OReadsPerGene.out.tab", skip = 4)

columnnames <- c("query_name", #1
                       "seed_eggNOG_ortholog", #2
                       "seed_ortholog_evalue", #3
                       "seed_ortholog_score",# 4
                       "best_tax_level", #5
                       "Preferred_name", #6
                       "GOs", #7	
                       "EC", #8
                       "KEGG_ko", #9
                       "KEGG_Pathway", #10
                       "KEGG_Module", #11
                       "KEGG_Reaction", #12
                       "KEGG_rclass", #13
                       "BRITE",	#14
                       "KEGG_TC", #15
                       "CAZy", #16
                       "BiGG_Reaction", #17
                       "taxonomic scope", #18
                       "eggNOG OGs", #19
                       "best eggNOG OG", #20
                       "COG_Functional", #21
                       "eggNOG free text desc.") #22

colnames(whiteshrimp) <- columnnames

#whiteshrimpidentifiers <- whiteshrimp[whiteshrimp$Preferred_name!="",c("query_name", "Preferred_name")]
whiteshrimpidentifiers <- whiteshrimp[whiteshrimp$Preferred_name!="",c("query_name", "Preferred_name")]
whiteshrimpidentifiers$query_name <- gsub(".*\\.", "", whiteshrimpidentifiers$query_name)
whiteshrimpidentifiers$Preferred_name <- toupper(whiteshrimpidentifiers$Preferred_name)
whiteshrimpidentifiers <- unique(whiteshrimpidentifiers)

#HP
colnames(whiteshrimpHP) <- c("query_name","left","right","whiteshrimp")
whiteshrimpHP <- whiteshrimpHP[whiteshrimpHP$left>10&whiteshrimpHP$right>10,c("query_name","whiteshrimp")]
whiteshrimpHP_mappings <- merge(whiteshrimpidentifiers,whiteshrimpHP,by="query_name")
whiteshrimpHP_mappingstotal <- aggregate(whiteshrimpHP_mappings$whiteshrimp~whiteshrimpHP_mappings$Preferred_name,data=whiteshrimpHP_mappings,FUN=sum)
whiteshrimpHP_mappingstotal <- data.table(whiteshrimpHP_mappingstotal[order(whiteshrimpHP_mappingstotal$`whiteshrimpHP_mappings$whiteshrimp`, decreasing = T),])
colnames(whiteshrimpHP_mappingstotal) <- c("gene","HP")

#HP2
colnames(whiteshrimpHP2) <- c("query_name","left","right","whiteshrimp")
whiteshrimpHP2 <- whiteshrimpHP2[whiteshrimpHP2$left>10&whiteshrimpHP2$right>10,c("query_name","whiteshrimp")]
whiteshrimpHP2_mappings <- merge(whiteshrimpidentifiers,whiteshrimpHP2,by="query_name")
whiteshrimpHP2_mappingstotal <- aggregate(whiteshrimpHP2_mappings$whiteshrimp~whiteshrimpHP2_mappings$Preferred_name,data=whiteshrimpHP2_mappings,FUN=sum)
whiteshrimpHP2_mappingstotal <- data.table(whiteshrimpHP2_mappingstotal[order(whiteshrimpHP2_mappingstotal$`whiteshrimpHP2_mappings$whiteshrimp`, decreasing = T),])
colnames(whiteshrimpHP2_mappingstotal) <- c("gene","HP2")

#ES
colnames(whiteshrimpES) <- c("query_name","left","right","whiteshrimp")
whiteshrimpES <- whiteshrimpES[whiteshrimpES$left>10&whiteshrimpES$right>10,c("query_name","whiteshrimp")]
whiteshrimpES_mappings <- merge(whiteshrimpidentifiers,whiteshrimpES,by="query_name")
whiteshrimpES_mappingstotal <- aggregate(whiteshrimpES_mappings$whiteshrimp~whiteshrimpES_mappings$Preferred_name,data=whiteshrimpES_mappings,FUN=sum)
whiteshrimpES_mappingstotal <- data.table(whiteshrimpES_mappingstotal[order(whiteshrimpES_mappingstotal$`whiteshrimpES_mappings$whiteshrimp`, decreasing = T),])
colnames(whiteshrimpES_mappingstotal) <- c("gene","ES")

#GILL
colnames(whiteshrimpGILL) <- c("query_name","left","right","whiteshrimp")
whiteshrimpGILL <- whiteshrimpGILL[whiteshrimpGILL$left>10&whiteshrimpGILL$right>10,c("query_name","whiteshrimp")]
whiteshrimpGILL_mappings <- merge(whiteshrimpidentifiers,whiteshrimpGILL,by="query_name")
whiteshrimpGILL_mappingstotal <- aggregate(whiteshrimpGILL_mappings$whiteshrimp~whiteshrimpGILL_mappings$Preferred_name,data=whiteshrimpGILL_mappings,FUN=sum)
whiteshrimpGILL_mappingstotal <- data.table(whiteshrimpGILL_mappingstotal[order(whiteshrimpGILL_mappingstotal$`whiteshrimpGILL_mappings$whiteshrimp`, decreasing = T),])
colnames(whiteshrimpGILL_mappingstotal) <- c("gene","GILL")

#GILL2
colnames(whiteshrimpGILL2) <- c("query_name","left","right","whiteshrimp")
whiteshrimpGILL2 <- whiteshrimpGILL2[whiteshrimpGILL2$left>10&whiteshrimpGILL2$right>10,c("query_name","whiteshrimp")]
whiteshrimpGILL2_mappings <- merge(whiteshrimpidentifiers,whiteshrimpGILL2,by="query_name")
whiteshrimpGILL2_mappingstotal <- aggregate(whiteshrimpGILL2_mappings$whiteshrimp~whiteshrimpGILL2_mappings$Preferred_name,data=whiteshrimpGILL2_mappings,FUN=sum)
whiteshrimpGILL2_mappingstotal <- data.table(whiteshrimpGILL2_mappingstotal[order(whiteshrimpGILL2_mappingstotal$`whiteshrimpGILL2_mappings$whiteshrimp`, decreasing = T),])
colnames(whiteshrimpGILL2_mappingstotal) <- c("gene","GILL2")

#MS
colnames(whiteshrimpMS) <- c("query_name","left","right","whiteshrimp")
whiteshrimpMS <- whiteshrimpMS[whiteshrimpMS$left>10&whiteshrimpMS$right>10,c("query_name","whiteshrimp")]
whiteshrimpMS_mappings <- merge(whiteshrimpidentifiers,whiteshrimpMS,by="query_name")
whiteshrimpMS_mappingstotal <- aggregate(whiteshrimpMS_mappings$whiteshrimp~whiteshrimpMS_mappings$Preferred_name,data=whiteshrimpMS_mappings,FUN=sum)
whiteshrimpMS_mappingstotal <- data.table(whiteshrimpMS_mappingstotal[order(whiteshrimpMS_mappingstotal$`whiteshrimpMS_mappings$whiteshrimp`, decreasing = T),])
colnames(whiteshrimpMS_mappingstotal) <- c("gene","MS")

#MS2
colnames(whiteshrimpMS2) <- c("query_name","left","right","whiteshrimp")
whiteshrimpMS2 <- whiteshrimpMS2[whiteshrimpMS2$left>10&whiteshrimpMS2$right>10,c("query_name","whiteshrimp")]
whiteshrimpMS2_mappings <- merge(whiteshrimpidentifiers,whiteshrimpMS2,by="query_name")
whiteshrimpMS2_mappingstotal <- aggregate(whiteshrimpMS2_mappings$whiteshrimp~whiteshrimpMS2_mappings$Preferred_name,data=whiteshrimpMS2_mappings,FUN=sum)
whiteshrimpMS2_mappingstotal <- data.table(whiteshrimpMS2_mappingstotal[order(whiteshrimpMS2_mappingstotal$`whiteshrimpMS2_mappings$whiteshrimp`, decreasing = T),])
colnames(whiteshrimpMS2_mappingstotal) <- c("gene","MS2")

#O
colnames(whiteshrimpO) <- c("query_name","left","right","whiteshrimp")
whiteshrimpO <- whiteshrimpO[whiteshrimpO$left>10&whiteshrimpO$right>10,c("query_name","whiteshrimp")]
whiteshrimpO_mappings <- merge(whiteshrimpidentifiers,whiteshrimpO,by="query_name")
whiteshrimpO_mappingstotal <- aggregate(whiteshrimpO_mappings$whiteshrimp~whiteshrimpO_mappings$Preferred_name,data=whiteshrimpO_mappings,FUN=sum)
whiteshrimpO_mappingstotal <- data.table(whiteshrimpO_mappingstotal[order(whiteshrimpO_mappingstotal$`whiteshrimpO_mappings$whiteshrimp`, decreasing = T),])
colnames(whiteshrimpO_mappingstotal) <- c("gene","O")

comparison <- Reduce(function(...) merge(..., all = TRUE, by = "gene"), 
                              list(whiteshrimpES_mappingstotal,
                                   whiteshrimpMS_mappingstotal,
                                   whiteshrimpMS2_mappingstotal,
                                   whiteshrimpGILL_mappingstotal,
                                   whiteshrimpGILL2_mappingstotal,
                                   whiteshrimpHP_mappingstotal,
                                   whiteshrimpHP2_mappingstotal,
                                   whiteshrimpO_mappingstotal))

comparison[is.na(comparison)] <- 0

head(comparison[order(-comparison$MS),], n=5)

colSums(whiteshrimpES[,c("whiteshrimp")]) #Known & unknown gene symbols