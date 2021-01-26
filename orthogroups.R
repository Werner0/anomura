#SNIPPETS (USE BEFORE SCRIPT)

library(data.table)
library(tibble)
library(plyr)
library(gridExtra)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/Orthogroups/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/Orthogroups/")}

#perspecies <- fread("Statistics_PerSpecies.tsv") #Cut file manually before loading
#colSums(perspecies[,2:10])

orthonumbers <- fread("Orthogroups.GeneCount.tsv")
#orthonumbers <- orthonumbers[apply(orthonumbers!=0, 1, all),]
#orthonumbers[, "sum"] <- orthonumbers$kingcrab.aa+orthonumbers$bluekingcrab-(orthonumbers$lobster.aa)-(orthonumbers$coconutcrab.aa)
#orthonumbers <- orthonumbers[apply(orthonumbers[,c(2:10)],1,function(x) all(x>0)),]
orthonumbers <- orthonumbers[apply(orthonumbers[,c(2:10)],1,function(x) all(x==1)),] 
#orthonumbers <- orthonumbers[orthonumbers$amphipod.aa<3&orthonumbers$isopod.aa<3,] 
orthonumbers <- orthonumbers[order(-orthonumbers$coconutcrab.aa),]

orthos <- fread("Orthogroups.tsv", sep="\t")

#singleorthos <- orthonumbers$Orthogroup
alltables_reduced <- fread("alltables_reduced.txt")
#singletable <- data.table()
#for (i in singleorthos) {

#query <- strsplit(gsub(" ","",orthos$bluekingcrab[orthos$Orthogroup=="OG000001"]),",")
#query <- query[[1]]
#subject <- alltables_reduced[alltables_reduced$query_name %in% query,]
#OGdata <- count(subject[subject$organism=="bluekingcrab",c("Preferred_name","eggNOG free text desc.")])
# OGplaceholder <- data.table(Preferred_name = NA,eggNOG.free.text.desc.=NA)
# if (nrow(OGdata)<1) {singletable <- rbind(singletable, OGplaceholder, fill = T)}
# singletable <- rbind(singletable, OGdata[,1:2], fill = T) }

#highcopygenes <- count(alltables_reduced[,c("Preferred_name","organism")])
#highcopygenes <- as.data.table(highcopygenes[order(-highcopygenes$freq),]) #Choose a gene from here for use in next line
#highcopyidentifiers <- alltables_reduced[alltables_reduced$Preferred_name %in% "KIF22"&alltables_reduced$organism=="coconutcrab",c("query_name")]
#orthogroupsofinterest <- orthos[orthos$coconutcrab.aa %in% highcopyidentifiers$query_name, c("Orthogroup")]

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
#library("flextable")
tree <- read.tree("SpeciesTree_rooted_node_labels.txt")
nodedata <- as.data.frame(tree$node.label)
nodedata$duplications <- c(758,846,740,1449,227,581,2935,11042)
nodedata$duplicationsl <- c("","","","","","","","")
colnames(nodedata) <- c("newick_label","duplications","duplicationsl")
groupInfo <- split(tree$tip.label, gsub("_\\w+", "", tree$tip.label))
names(groupInfo) <- gsub("\\..*","",names(groupInfo))
names(groupInfo) <- gsub("\\-.*","",names(groupInfo))
names(groupInfo) <- c("A. vulgare", "P. hawaiensis", "Achelata", "P. virginalis", "P. trituberculatus", "Anomura", "Anomura", "Anomura", "L. vannamei")
tree$tip.label <- c("Armadillidium vulgare", "Parhyale hawaiensis", "Panulirus ornatus", "Procambarus virginalis", "Portunus trituberculatus", "Birgus latro", "Paralithodes camtschaticus", "Paralithodes platypus", "Litopenaeus vannamei")
tree <- groupOTU(tree, groupInfo)
#tree <- root(tree, node = 011, edgelabel = TRUE)
g <- ggtree(tree, size=1.5) %<+% nodedata + geom_treescale() + xlim(NA, 1.1)
g <- rotate(g,14) #Lobster
g <- rotate(g,16) #Anomura
g2 <- g + geom_label(aes(label = duplicationsl, fill = duplications), show.legend = FALSE) +
  theme(legend.position = NULL) + scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu")) +
  #geom_cladelabel(node=12, label="Pleocyemata", color="red2", offset=0.6, align=TRUE, angle = 90, offset.text = 0.05, hjust = 0.5, barsize = 1) + 
  #geom_cladelabel(node=17, label="Anomura", color="red2", offset=0.8, align=TRUE, angle = 90, offset.text = 0.05, hjust = 0.5, barsize = 1) + 
  theme_tree2()
#geom_text(aes(subset=(node==508), label = italic('Acetobacter spp.')), parse=TRUE, colour="blue", hjust=-.02) 

#(isopod.aa_6304:0.253617,(amphipod.aa_8566:0.638071,(((lobster.aa_49883:0.202489,marbled.aa_11462:0.19829)N4_227:0.0536186,(swimmingcrab.aa_5873:0.346618,(coconutcrab.aa_85938:0.181323,(kingcrab.aa_78264:0.00540666,bluekingcrab_84676:0.0272019)N7_11042:0.215072)N6_2935:0.188058)N5_581:0.0669401)N3_1449:0.0654401,whiteshrimp.aa_4397:0.293433)N2_740:0.234702)N1_846:0.253617)N0_758;
