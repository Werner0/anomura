#SNIPPETS (USE BEFORE SCRIPT)
#New gene combination in second section of this script
#awk -F\# '$1!="" { print $0 ;} ' lobster.emapper.annotations > lobster.emapper_nocomments.annotations

library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/eggnog_perRead/Lobster/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/eggnog_perRead/Lobster/")}

library(plyr)
library(tibble)

lobster <- fread("lobster.emapper_nocomments.annotations")
lobsterHP <- fread("HPReadsPerGene.out.tab", skip = 4)
lobsterES <- fread("ESReadsPerGene.out.tab", skip = 4)
lobsterGILL <- fread("GReadsPerGene.out.tab", skip = 4)
lobsterOV <- fread("OReadsPerGene.out.tab", skip = 4)
lobsterMS <- fread("MReadsPerGene.out.tab", skip = 4)

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

colnames(lobster) <- columnnames

#lobsteridentifiers <- lobster[lobster$Preferred_name!="",c("query_name", "Preferred_name")]
lobsteridentifiers <- lobster[lobster$Preferred_name=="",c("query_name", "Preferred_name")]
lobsteridentifiers$query_name <- gsub("\\..*", "", lobsteridentifiers$query_name)
lobsteridentifiers$Preferred_name <- toupper(lobsteridentifiers$Preferred_name)
lobsteridentifiers <- unique(lobsteridentifiers)

#HP
colnames(lobsterHP) <- c("query_name","left","right","lobster")
lobsterHP <- lobsterHP[lobsterHP$left>10&lobsterHP$right>10,c("query_name","lobster")]
lobsterHP_mappings <- merge(lobsteridentifiers,lobsterHP,by="query_name")
lobsterHP_mappingstotal <- aggregate(lobsterHP_mappings$lobster~lobsterHP_mappings$Preferred_name,data=lobsterHP_mappings,FUN=sum)
lobsterHP_mappingstotal <- data.table(lobsterHP_mappingstotal[order(lobsterHP_mappingstotal$`lobsterHP_mappings$lobster`, decreasing = T),])
colnames(lobsterHP_mappingstotal) <- c("gene","HP")

#ES
colnames(lobsterES) <- c("query_name","left","right","lobster")
lobsterES <- lobsterES[lobsterES$left>10&lobsterES$right>10,c("query_name","lobster")]
lobsterES_mappings <- merge(lobsteridentifiers,lobsterES,by="query_name")
lobsterES_mappingstotal <- aggregate(lobsterES_mappings$lobster~lobsterES_mappings$Preferred_name,data=lobsterES_mappings,FUN=sum)
lobsterES_mappingstotal <- data.table(lobsterES_mappingstotal[order(lobsterES_mappingstotal$`lobsterES_mappings$lobster`, decreasing = T),])
colnames(lobsterES_mappingstotal) <- c("gene","ES")

#GILL
colnames(lobsterGILL) <- c("query_name","left","right","lobster")
lobsterGILL <- lobsterGILL[lobsterGILL$left>10&lobsterGILL$right>10,c("query_name","lobster")]
lobsterGILL_mappings <- merge(lobsteridentifiers,lobsterGILL,by="query_name")
lobsterGILL_mappingstotal <- aggregate(lobsterGILL_mappings$lobster~lobsterGILL_mappings$Preferred_name,data=lobsterGILL_mappings,FUN=sum)
lobsterGILL_mappingstotal <- data.table(lobsterGILL_mappingstotal[order(lobsterGILL_mappingstotal$`lobsterGILL_mappings$lobster`, decreasing = T),])
colnames(lobsterGILL_mappingstotal) <- c("gene","GILL")

#OV
colnames(lobsterOV) <- c("query_name","left","right","lobster")
lobsterOV <- lobsterOV[lobsterOV$left>10&lobsterOV$right>10,c("query_name","lobster")]
lobsterOV_mappings <- merge(lobsteridentifiers,lobsterOV,by="query_name")
lobsterOV_mappingstotal <- aggregate(lobsterOV_mappings$lobster~lobsterOV_mappings$Preferred_name,data=lobsterOV_mappings,FUN=sum)
lobsterOV_mappingstotal <- data.table(lobsterOV_mappingstotal[order(lobsterOV_mappingstotal$`lobsterOV_mappings$lobster`, decreasing = T),])
colnames(lobsterOV_mappingstotal) <- c("gene","OV")

#MS
colnames(lobsterMS) <- c("query_name","left","right","lobster")
lobsterMS <- lobsterMS[lobsterMS$left>10&lobsterMS$right>10,c("query_name","lobster")]
lobsterMS_mappings <- merge(lobsteridentifiers,lobsterMS,by="query_name")
lobsterMS_mappingstotal <- aggregate(lobsterMS_mappings$lobster~lobsterMS_mappings$Preferred_name,data=lobsterMS_mappings,FUN=sum)
lobsterMS_mappingstotal <- data.table(lobsterMS_mappingstotal[order(lobsterMS_mappingstotal$`lobsterMS_mappings$lobster`, decreasing = T),])
colnames(lobsterMS_mappingstotal) <- c("gene","MS")

comparison <- Reduce(function(...) merge(..., all = TRUE, by = "gene"), 
                              list(lobsterES_mappingstotal,
                                   lobsterMS_mappingstotal,
                                   lobsterGILL_mappingstotal,
                                   lobsterOV_mappingstotal,
                                   lobsterHP_mappingstotal))

comparison[is.na(comparison)] <- 0

head(comparison[order(-comparison$GILL),], n=5)
