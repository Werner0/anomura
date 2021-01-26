#SNIPPETS (USE BEFORE SCRIPT)
#New gene combination in second section of this script
#awk -F\# '$1!="" { print $0 ;} ' kingcrab.emapper.annotations > kingcrab.emapper_nocomments.annotations

library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/eggnog_perRead/Kingcrab//") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/eggnog_perRead/Kingcrab//")}

library(plyr)
library(tibble)

kingcrab <- fread("kingcrab.emapper_nocomments.annotations")
kingcrabHP <- fread("HPReadsPerGene.out.tab", skip = 4)
kingcrabES <- fread("ESReadsPerGene.out.tab", skip = 4)
kingcrabGILL <- fread("GReadsPerGene.out.tab", skip = 4)
kingcrabMS <- fread("MReadsPerGene.out.tab", skip = 4)

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

colnames(kingcrab) <- columnnames

#kingcrabidentifiers <- kingcrab[kingcrab$Preferred_name!="",c("query_name", "Preferred_name")]
kingcrabidentifiers <- kingcrab[kingcrab$Preferred_name=="",c("query_name", "Preferred_name")]
kingcrabidentifiers$query_name <- gsub("\\..*", "", kingcrabidentifiers$query_name)
kingcrabidentifiers$Preferred_name <- toupper(kingcrabidentifiers$Preferred_name)
kingcrabidentifiers <- unique(kingcrabidentifiers)

#HP
colnames(kingcrabHP) <- c("query_name","left","right","kingcrab")
kingcrabHP <- kingcrabHP[kingcrabHP$left>10&kingcrabHP$right>10,c("query_name","kingcrab")]
kingcrabHP_mappings <- merge(kingcrabidentifiers,kingcrabHP,by="query_name")
kingcrabHP_mappingstotal <- aggregate(kingcrabHP_mappings$kingcrab~kingcrabHP_mappings$Preferred_name,data=kingcrabHP_mappings,FUN=sum)
kingcrabHP_mappingstotal <- data.table(kingcrabHP_mappingstotal[order(kingcrabHP_mappingstotal$`kingcrabHP_mappings$kingcrab`, decreasing = T),])
colnames(kingcrabHP_mappingstotal) <- c("gene","HP")

#ES
colnames(kingcrabES) <- c("query_name","left","right","kingcrab")
kingcrabES <- kingcrabES[kingcrabES$left>10&kingcrabES$right>10,c("query_name","kingcrab")]
kingcrabES_mappings <- merge(kingcrabidentifiers,kingcrabES,by="query_name")
kingcrabES_mappingstotal <- aggregate(kingcrabES_mappings$kingcrab~kingcrabES_mappings$Preferred_name,data=kingcrabES_mappings,FUN=sum)
kingcrabES_mappingstotal <- data.table(kingcrabES_mappingstotal[order(kingcrabES_mappingstotal$`kingcrabES_mappings$kingcrab`, decreasing = T),])
colnames(kingcrabES_mappingstotal) <- c("gene","ES")

#GILL
colnames(kingcrabGILL) <- c("query_name","left","right","kingcrab")
kingcrabGILL <- kingcrabGILL[kingcrabGILL$left>10&kingcrabGILL$right>10,c("query_name","kingcrab")]
kingcrabGILL_mappings <- merge(kingcrabidentifiers,kingcrabGILL,by="query_name")
kingcrabGILL_mappingstotal <- aggregate(kingcrabGILL_mappings$kingcrab~kingcrabGILL_mappings$Preferred_name,data=kingcrabGILL_mappings,FUN=sum)
kingcrabGILL_mappingstotal <- data.table(kingcrabGILL_mappingstotal[order(kingcrabGILL_mappingstotal$`kingcrabGILL_mappings$kingcrab`, decreasing = T),])
colnames(kingcrabGILL_mappingstotal) <- c("gene","GILL")

#MS
colnames(kingcrabMS) <- c("query_name","left","right","kingcrab")
kingcrabMS <- kingcrabMS[kingcrabMS$left>10&kingcrabMS$right>10,c("query_name","kingcrab")]
kingcrabMS_mappings <- merge(kingcrabidentifiers,kingcrabMS,by="query_name")
kingcrabMS_mappingstotal <- aggregate(kingcrabMS_mappings$kingcrab~kingcrabMS_mappings$Preferred_name,data=kingcrabMS_mappings,FUN=sum)
kingcrabMS_mappingstotal <- data.table(kingcrabMS_mappingstotal[order(kingcrabMS_mappingstotal$`kingcrabMS_mappings$kingcrab`, decreasing = T),])
colnames(kingcrabMS_mappingstotal) <- c("gene","MS")

comparison <- Reduce(function(...) merge(..., all = TRUE, by = "gene"), 
                              list(kingcrabES_mappingstotal,
                                   kingcrabMS_mappingstotal,
                                   kingcrabGILL_mappingstotal,
                                   kingcrabHP_mappingstotal))

comparison[is.na(comparison)] <- 0

head(comparison[order(-comparison$MS),], n=5)
