#SNIPPETS (USE BEFORE SCRIPT)
#New gene combination in second section of this script
#awk -F\# '$1!="" { print $0 ;} ' coconutcrab.emapper.annotations > coconutcrab.emapper_nocomments.annotations

library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/eggnog_perRead/Coconutcrab//") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/eggnog_perRead/Coconutcrab//")}

library(plyr)
library(tibble)

coconutcrab <- fread("coconutcrab.emapper_nocomments.annotations")
coconutcrabHP <- fread("HPReadsPerGene.out.tab", skip = 4)
coconutcrabES <- fread("ESReadsPerGene.out.tab", skip = 4)
coconutcrabGILL <- fread("GReadsPerGene.out.tab", skip = 4)
coconutcrabMS <- fread("MReadsPerGene.out.tab", skip = 4)

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

colnames(coconutcrab) <- columnnames

#coconutcrabidentifiers <- coconutcrab[coconutcrab$Preferred_name!="",c("query_name", "Preferred_name")]
coconutcrabidentifiers <- coconutcrab[coconutcrab$Preferred_name=="",c("query_name", "Preferred_name")]
coconutcrabidentifiers$query_name <- gsub("\\..*", "", coconutcrabidentifiers$query_name)
coconutcrabidentifiers$Preferred_name <- toupper(coconutcrabidentifiers$Preferred_name)
coconutcrabidentifiers <- unique(coconutcrabidentifiers)

#HP
colnames(coconutcrabHP) <- c("query_name","left","right","coconutcrab")
coconutcrabHP <- coconutcrabHP[coconutcrabHP$left>10&coconutcrabHP$right>10,c("query_name","coconutcrab")]
coconutcrabHP_mappings <- merge(coconutcrabidentifiers,coconutcrabHP,by="query_name")
coconutcrabHP_mappingstotal <- aggregate(coconutcrabHP_mappings$coconutcrab~coconutcrabHP_mappings$Preferred_name,data=coconutcrabHP_mappings,FUN=sum)
coconutcrabHP_mappingstotal <- data.table(coconutcrabHP_mappingstotal[order(coconutcrabHP_mappingstotal$`coconutcrabHP_mappings$coconutcrab`, decreasing = T),])
colnames(coconutcrabHP_mappingstotal) <- c("gene","HP")

#ES
colnames(coconutcrabES) <- c("query_name","left","right","coconutcrab")
coconutcrabES <- coconutcrabES[coconutcrabES$left>10&coconutcrabES$right>10,c("query_name","coconutcrab")]
coconutcrabES_mappings <- merge(coconutcrabidentifiers,coconutcrabES,by="query_name")
coconutcrabES_mappingstotal <- aggregate(coconutcrabES_mappings$coconutcrab~coconutcrabES_mappings$Preferred_name,data=coconutcrabES_mappings,FUN=sum)
coconutcrabES_mappingstotal <- data.table(coconutcrabES_mappingstotal[order(coconutcrabES_mappingstotal$`coconutcrabES_mappings$coconutcrab`, decreasing = T),])
colnames(coconutcrabES_mappingstotal) <- c("gene","ES")

#GILL
colnames(coconutcrabGILL) <- c("query_name","left","right","coconutcrab")
coconutcrabGILL <- coconutcrabGILL[coconutcrabGILL$left>10&coconutcrabGILL$right>10,c("query_name","coconutcrab")]
coconutcrabGILL_mappings <- merge(coconutcrabidentifiers,coconutcrabGILL,by="query_name")
coconutcrabGILL_mappingstotal <- aggregate(coconutcrabGILL_mappings$coconutcrab~coconutcrabGILL_mappings$Preferred_name,data=coconutcrabGILL_mappings,FUN=sum)
coconutcrabGILL_mappingstotal <- data.table(coconutcrabGILL_mappingstotal[order(coconutcrabGILL_mappingstotal$`coconutcrabGILL_mappings$coconutcrab`, decreasing = T),])
colnames(coconutcrabGILL_mappingstotal) <- c("gene","GILL")

#MS
colnames(coconutcrabMS) <- c("query_name","left","right","coconutcrab")
coconutcrabMS <- coconutcrabMS[coconutcrabMS$left>10&coconutcrabMS$right>10,c("query_name","coconutcrab")]
coconutcrabMS_mappings <- merge(coconutcrabidentifiers,coconutcrabMS,by="query_name")
coconutcrabMS_mappingstotal <- aggregate(coconutcrabMS_mappings$coconutcrab~coconutcrabMS_mappings$Preferred_name,data=coconutcrabMS_mappings,FUN=sum)
coconutcrabMS_mappingstotal <- data.table(coconutcrabMS_mappingstotal[order(coconutcrabMS_mappingstotal$`coconutcrabMS_mappings$coconutcrab`, decreasing = T),])
colnames(coconutcrabMS_mappingstotal) <- c("gene","MS")

comparison <- Reduce(function(...) merge(..., all = TRUE, by = "gene"), 
                              list(coconutcrabES_mappingstotal,
                                   coconutcrabMS_mappingstotal,
                                   coconutcrabGILL_mappingstotal,
                                   coconutcrabHP_mappingstotal))

comparison[is.na(comparison)] <- 0

head(comparison[order(-comparison$MS),], n=5)