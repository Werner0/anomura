#SNIPPETS (USE BEFORE SCRIPT)

library(data.table)
library(tibble)
library(UpSetR)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/mTP_anno/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/mTP_anno/")}

library(plyr)

lobster <- fread("lobster.emapper_nocomments.annotations")
kingcrab <- fread("kingcrab.emapper_nocomments.annotations")
coconutcrab <- fread("coconutcrab.emapper_nocomments.annotations")
swimmingcrab <- fread("swimmingcrab_nocomments_new.annotations")
marbledcrayfish <- fread("marbledcrayfish_nocomments_new.annotations")
bluekingcrab <- fread("bluekingcrab_nocomments_new.annotations")
whiteshrimp <- fread("whiteshrimp_nocomments_new.annotations")

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
colnames(kingcrab) <- columnnames
colnames(coconutcrab) <- columnnames
colnames(swimmingcrab) <- columnnames
colnames(bluekingcrab) <- columnnames
colnames(marbledcrayfish) <- columnnames
colnames(whiteshrimp) <- columnnames

#singleorthologs <- singleorthologs[,c("query_name","best_tax_level","Preferred_name","taxonomic scope","COG_Functional","eggNOG free text desc.")]
#single_arthropoda <- singleorthologs[singleorthologs$`taxonomic scope`=="Arthropoda",]

taxonomiccomparison <- Reduce(function(...) merge(..., all = TRUE, by = "x"), 
                              list(count(lobster$`taxonomic scope`),
                                   count(kingcrab$`taxonomic scope`),
                                   count(coconutcrab$`taxonomic scope`),
                                   count(swimmingcrab$`taxonomic scope`),
                                   count(marbledcrayfish$`taxonomic scope`),
                                   count(bluekingcrab$`taxonomic scope`),
                                   count(whiteshrimp$`taxonomic scope`)))

KOGbreakdown <- Reduce(function(...) merge(..., all = TRUE, by = "x"), 
                              list(count(lobster$COG_Functional),
                                   count(kingcrab$COG_Functional),
                                   count(coconutcrab$COG_Functional),
                                   count(swimmingcrab$COG_Functional),
                                   count(marbledcrayfish$COG_Functional),
                                   count(bluekingcrab$COG_Functional),
                                   count(whiteshrimp$COG_Functional)))

genecomparison <- Reduce(function(...) merge(..., all = TRUE, by = "x"), 
                              list(count(toupper(lobster$Preferred_name)),
                                   count(toupper(kingcrab$Preferred_name)),
                                   count(toupper(coconutcrab$Preferred_name)),
                                   count(toupper(swimmingcrab$Preferred_name)),
                                   count(toupper(marbledcrayfish$Preferred_name)),
                                   count(toupper(bluekingcrab$Preferred_name)),
                                   count(toupper(whiteshrimp$Preferred_name))))

# genecomparison <- Reduce(function(...) merge(..., all = TRUE, by = "x"), 
#                           list(count(toupper(lobster$Preferred_name)),
#                                count(toupper(kingcrab$Preferred_name)),
#                                count(toupper(coconutcrab$Preferred_name))))

renamecolumns <- c("Taxon",
                   "Lobster",
                   "Kingcrab",
                   "Coconutcrab",
                   "Swimmingcrab",
                   "Marbledcrayfish",
                   "Bluekingcrab",
                   "Whiteshrimp")

taxonomiccomparison[is.na(taxonomiccomparison)] <- 0
genecomparison[is.na(genecomparison)] <- 0
KOGbreakdown[is.na(KOGbreakdown)] <- 0
colnames(taxonomiccomparison) <- renamecolumns
colnames(genecomparison) <- renamecolumns
colnames(genecomparison)[1] <- "Gene"
colnames(KOGbreakdown) <- renamecolumns
colnames(KOGbreakdown)[1] <- "KOG_category"
gc <- genecomparison[,c(1:8)]
#GOI selection criteria
gc <- gc[!(gc$Lobster==0&gc$Kingcrab==0&gc$Coconutcrab==0&gc$Swimmingcrab==0&gc$Marbledcrayfish==0&gc$Bluekingcrab==0&gc$Whiteshrimp==0),]
gc <- gc[!(gc$Gene==""),] #####Excluded this line for ("gc") file written below


#CAFE 
gc[, "max"] <- apply(gc[, 2:8], 1, max)
gc[, "min"] <- apply(gc[, 2:8], 1, min)
gc[, "sum"] <- apply(gc[, 2:8], 1, sum)
gc$diff <- gc$max - gc$min
gc <- gc[gc$diff<10000&gc$sum>0,c(1:8)]

#Intersection plot
tt <- as.data.frame(gc[,c(2:8)])
tt[tt > 1] <- 1
tt <- cbind(gc[,1],tt)
colnames(tt) <- c("gene","Panulirus ornatus","Paralithodes camtschaticus","Birgus latro","Portunus trituberculatus","Procambarus virginalis","Paralithodes platypus","Litopenaeus vannamei")
tt[,"sum"] <- tt$`Panulirus ornatus`+tt$`Paralithodes camtschaticus`+tt$`Birgus latro`+tt$`Portunus trituberculatus`+tt$`Procambarus virginalis`+tt$`Paralithodes platypus`+tt$`Litopenaeus vannamei`
tt <- tt[tt$sum>1,]
upset(tt, sets = c("Litopenaeus vannamei","Procambarus virginalis","Panulirus ornatus","Portunus trituberculatus","Paralithodes platypus","Paralithodes camtschaticus","Birgus latro"), 
      order.by = c("freq","degree"),
      mainbar.y.label = "Number of mTPs in intersection",
      sets.x.label = "Total mTPs in intersections",nintersects = NA,
      show.numbers = FALSE,
      point.size = 2,
      keep.order = TRUE,
      main.bar.color=c("#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A","#FFC20A",
                       "#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC","#0C7BDC",
                       "#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A","#E1BE6A",
                       "#40B0A6","#40B0A6","#40B0A6","#40B0A6","#40B0A6","#40B0A6","#40B0A6","#40B0A6","#40B0A6","#40B0A6",
                       "#FEFE62","#FEFE62","#FEFE62"))
#fwrite(taxonomiccomparison, file = "Eggnog_taxa.csv", sep = ",")
#fwrite(KOGbreakdown, file = "KOG_breakdown.csv", sep = ",")
#fwrite(gc, file = "Eggnog_homologs.csv", sep = ",")