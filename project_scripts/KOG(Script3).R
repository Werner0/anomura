#SNIPPETS (USE BEFORE SCRIPT)
#New gene combination in second section of this script
#awk -F\# '$1!="" { print $0 ;} ' lobster.emapper.annotations > lobster.emapper_nocomments.annotations

library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/eggnog") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/eggnog")}

library(plyr)
library(tibble)

lobster <- fread("lobster_noHP.emapper_nocomments.annotations")
lobster_ab_initio <- fread("lobster_noHP_abinitio_UTR.emapper_nocomments.annotations")
kingcrab <- fread("kingcrab.emapper_nocomments.annotations")
kingcrab_ab_initio <- fread("kingcrab_abinitio.emapper_nocomments.annotations")
coconutcrab <- fread("coconutcrab.emapper_nocomments.annotations")
coconutcrab_ab_initio <- fread("coconutcrab_abinitio.emapper_nocomments.annotations")
bluekingcrab <- fread("bluekingcrab.emapper_nocomments.annotations")
whiteshrimp <- fread("whiteshrimp.emapper_nocomments.annotations")
#mitten <- fread("mitten.emapper_nocomments.annotations")
marbled <- fread("marbled.emapper_nocomments.annotations")
isopod <- fread("isopod.emapper_nocomments.annotations")
amphipod <- fread("amphipod.emapper_nocomments.annotations")
swimmingcrab <- fread("swimmingcrab_new.emapper_nocomments.annotations")
#singleorthologs <- fread("singles.emapper_nocomments.annotations")


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
lobster$organism <- "lobster"
colnames(lobster_ab_initio) <- columnnames
colnames(kingcrab) <- columnnames
kingcrab$organism <- "redkingcrab"
colnames(kingcrab_ab_initio) <- columnnames
colnames(coconutcrab) <- columnnames
coconutcrab$organism <- "coconutcrab"
colnames(coconutcrab_ab_initio) <- columnnames
colnames(bluekingcrab) <- columnnames
bluekingcrab$organism <- "bluekingcrab"
colnames(whiteshrimp) <- columnnames
whiteshrimp$organism <- "whiteshrimp"
colnames(marbled) <- columnnames
marbled$organism <- "marbledcrayfish"
#colnames(mitten) <- columnnames
colnames(swimmingcrab) <- columnnames
swimmingcrab$organism <- "swimmingcrab"
colnames(isopod) <- columnnames
isopod$organism <- "isopod"
colnames(amphipod) <- columnnames
amphipod$organism <- "amphipod"

alldata <- rbind(lobster,kingcrab,coconutcrab,bluekingcrab,marbled,whiteshrimp,amphipod,isopod,swimmingcrab)
alldata <- alldata[,c("query_name","Preferred_name","COG_Functional","taxonomic scope","organism")]
#fwrite(alldata, file = "alldata.txt", sep = "\t")

#singleorthologs <- singleorthologs[,c("query_name","best_tax_level","Preferred_name","taxonomic scope","COG_Functional","eggNOG free text desc.")]
#single_arthropoda <- singleorthologs[singleorthologs$`taxonomic scope`=="Arthropoda",]

taxonomiccomparison <- Reduce(function(...) merge(..., all = TRUE, by = "x"), 
                              list(count(lobster$`taxonomic scope`),
                                   count(kingcrab$`taxonomic scope`),
                                   count(coconutcrab$`taxonomic scope`),
                                   count(lobster_ab_initio$`taxonomic scope`),
                                   count(kingcrab_ab_initio$`taxonomic scope`),
                                   count(coconutcrab_ab_initio$`taxonomic scope`),
                                   count(bluekingcrab$`taxonomic scope`),
                                   count(whiteshrimp$`taxonomic scope`),
                                   count(marbled$`taxonomic scope`),
                                   #count(mitten$`taxonomic scope`),
                                   count(swimmingcrab$`taxonomic scope`),
                                   count(amphipod$`taxonomic scope`),
                                   count(isopod$`taxonomic scope`)))

KOGbreakdown <- Reduce(function(...) merge(..., all = TRUE, by = "x"), 
                              list(count(lobster$COG_Functional),
                                   count(kingcrab$COG_Functional),
                                   count(coconutcrab$COG_Functional),
                                   count(lobster_ab_initio$COG_Functional),
                                   count(kingcrab_ab_initio$COG_Functional),
                                   count(coconutcrab_ab_initio$COG_Functional),
                                   count(bluekingcrab$COG_Functional),
                                   count(whiteshrimp$COG_Functional),
                                   count(marbled$COG_Functional),
                                   #count(mitten$COG_Functional),
                                   count(swimmingcrab$COG_Functional),
                                   count(amphipod$COG_Functional),
                                   count(isopod$COG_Functional)))

genecomparison <- Reduce(function(...) merge(..., all = TRUE, by = "x"), 
                              list(count(toupper(lobster$Preferred_name)),
                                   count(toupper(kingcrab$Preferred_name)),
                                   count(toupper(coconutcrab$Preferred_name)),
                                   count(toupper(lobster_ab_initio$Preferred_name)),
                                   count(toupper(kingcrab_ab_initio$Preferred_name)),
                                   count(toupper(coconutcrab_ab_initio$Preferred_name)),
                                   count(toupper(bluekingcrab$Preferred_name)),
                                   count(toupper(whiteshrimp$Preferred_name)),
                                   count(toupper(marbled$Preferred_name)),
                                   #count(toupper(mitten$Preferred_name)),
                                   count(toupper(swimmingcrab$Preferred_name)),
                                   count(toupper(amphipod$Preferred_name)),
                                   count(toupper(isopod$Preferred_name))))

genecomparison <- Reduce(function(...) merge(..., all = TRUE, by = "x"), 
                         list(count(toupper(lobster$Preferred_name)),
                              count(toupper(kingcrab$Preferred_name)),
                              count(toupper(coconutcrab$Preferred_name)),
                              count(toupper(lobster_ab_initio$Preferred_name)),
                              count(toupper(kingcrab_ab_initio$Preferred_name)),
                              count(toupper(coconutcrab_ab_initio$Preferred_name)),
                              count(toupper(bluekingcrab$Preferred_name)),
                              count(toupper(whiteshrimp$Preferred_name)),
                              count(toupper(marbled$Preferred_name)),
                              #count(toupper(mitten$Preferred_name)),
                              count(toupper(swimmingcrab$Preferred_name)),
                              count(toupper(amphipod$Preferred_name)),
                              count(toupper(isopod$Preferred_name))))

renamecolumns <- c("Taxon",
                   "Lobster_RNA",
                   "Kingcrab_RNA",
                   "Coconutcrab_RNA",
                   "Lobster_ab",
                   "Kingcrab_ab",
                   "Coconutcrab_ab",
                   "Bluekingcrab",
                   "Whiteshrimp",
                   "Marbledcrayfish",
                   #"Mittencrab",
                   "Swimmingcrab",
                   "Amphipod",
                   "Isopod")

taxonomiccomparison[is.na(taxonomiccomparison)] <- 0
genecomparison[is.na(genecomparison)] <- 0
KOGbreakdown[is.na(KOGbreakdown)] <- 0
colnames(taxonomiccomparison) <- renamecolumns
colnames(genecomparison) <- renamecolumns
colnames(genecomparison)[1] <- "Gene"
colnames(KOGbreakdown) <- renamecolumns
colnames(KOGbreakdown)[1] <- "KOG_category"
gc <- genecomparison[,c(1:4,8:13)]
#GOI selection criteria
gc <- gc[!(gc$Lobster_RNA==0&gc$Kingcrab_RNA==0&gc$Coconutcrab_RNA==0&gc$Bluekingcrab==0&gc$Whiteshrimp==0&gc$Marbledcrayfish==0&gc$Swimmingcrab==0&gc$Amphipod==0&gc$Isopod==0),]


gc <- gc[!(gc$Gene==""),] #####Include this line for cafe ("gc") file written below


#For CAFE
colnames(gc)[1] <- "Desc"
gc <- add_column(gc, "Family ID" = "test", .after = "Desc")
gc$`Family ID` <- gsub("test", NA, gc$`Family ID`)
gc[, "max"] <- apply(gc[, 3:11], 1, max)
gc[, "min"] <- apply(gc[, 3:11], 1, min)
gc[, "sum"] <- apply(gc[, 3:11], 1, sum)
gc$diff <- gc$max - gc$min
gc <- gc[gc$sum>1&gc$diff>99,c(1:11)]
gc$`Family ID` <- 1:nrow(gc)
colnames(gc) <- c("Desc", "Family ID", "lobster","kingcrab","coconutcrab","bluekingcrab","whiteshrimp","marbledcrayfish","swimmingcrab","amphipod","isopod")
#gc <- gc[gc$amphipod>0&gc$coconutcrab>0&gc$isopod>0
#                      &gc$kingcrab>0&gc$bluekingcrab>0&gc$lobster>0&gc$marbledcrayfish>0
#                      &gc$swimmingcrab>0&gc$whiteshrimp>0,]
#fwrite(gc, file = "outliercoding.csv", sep = ",")

#fwrite(taxonomiccomparison, file = "Final_Eggnog_taxa.csv", sep = ",")
#fwrite(KOGbreakdown, file = "Final_KOG_breakdown.csv", sep = ",")
#fwrite(gc, file = "Final_Eggnog_homologs_with_no_symbols.csv", sep = ",")
