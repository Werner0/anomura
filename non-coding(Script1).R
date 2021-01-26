#SNIPPETS (RUN BEFORE SCRIPT)
#awk '$1 !~ /#/ {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' mrum-genome.tblout > organism.awk_parsed.mrum-genome.tblout

library(data.table)
library(tibble)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/rfam") {
setwd("/Users/wernerveldsman/Desktop/Rsources/rfam")}

noncoding_results <- fread("results_noncoding/Base_clade_results.txt")
coding_results <- fread("results_coding/Base_clade_results.txt")
results_merge <- merge(coding_results,noncoding_results,by="#Taxon_ID")
#fwrite(results_merge, file="family_analysis.csv", sep = ",")

coconutcrab <- fread("coconutcrab.awk_parsed.mrum-genome.tblout")
kingcrab <- fread("kingcrab.awk_parsed.mrum-genome.tblout")
bluekingcrab <- fread("bluekingcrab.awk_parsed.mrum-genome.tblout")
lobster <- fread("lobster.awk_parsed.mrum-genome.tblout")
mitten <- fread("mittencrab.awk_parsed.mrum-genome_2.tblout")
mitten_published <- fread("mittencrab_NCBI.awk_parsed.mrum-genome.tblout")
whiteshrimp <- fread("whiteshrimp.awk_parsed.mrum-genome.tblout")
isopod <- fread("isopod.awk_parsed.mrum-genome.tblout")
amphipod <- fread("amphipod.awk_parsed.mrum-genome.tblout")
swimmingcrab <- fread("swimmingcrab_new.awk_parsed.mrum-genome.tblout")
marbledcrayfish <- fread("marbled_crayfish.awk_parsed.mrum-genome.tblout")

coconutcrab <- coconutcrab[,c(1,2,3,9,10,11,12,17)]
lobster <- lobster[,c(1,2,3,9,10,11,12,17)]
kingcrab <- kingcrab[,c(1,2,3,9,10,11,12,17)]
bluekingcrab <- bluekingcrab[,c(1,2,3,9,10,11,12,17)]
mitten <- mitten[,c(1,2,3,9,10,11,12,17)]
mitten_published <- mitten_published[,c(1,2,3,9,10,11,12,17)]
whiteshrimp <- whiteshrimp[,c(1,2,3,9,10,11,12,17)]
isopod <- isopod[,c(1,2,3,9,10,11,12,17)]
amphipod <- amphipod[,c(1,2,3,9,10,11,12,17)]
swimmingcrab <- swimmingcrab[,c(1,2,3,9,10,11,12,17)]
marbledcrayfish <- marbledcrayfish[,c(1,2,3,9,10,11,12,17)]

colnames(coconutcrab) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(kingcrab) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(bluekingcrab) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(lobster) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(mitten) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(mitten_published) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(whiteshrimp) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(isopod) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(amphipod) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(swimmingcrab) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")
colnames(marbledcrayfish) <- c("name", "accession", "query","start","stop","strand","truncated","evalue")


coconutcrab <- cbind("organism" = "coconutcrab", coconutcrab)
kingcrab <- cbind("organism" = "kingcrab", kingcrab)
bluekingcrab <- cbind("organism" = "bluekingcrab", bluekingcrab)
lobster <- cbind("organism" = "lobster", lobster)
mitten <- cbind("organism" = "mittencrab", mitten)
mitten_published <- cbind("organism" = "mittencrab_published", mitten_published)
whiteshrimp <- cbind("organism" = "whiteshrimp", whiteshrimp)
isopod <- cbind("organism" = "isopod", isopod)
amphipod <- cbind("organism" = "amphipod", amphipod)
swimmingcrab <- cbind("organism" = "swimmingcrab", swimmingcrab)
marbledcrayfish <- cbind("organism" = "marbledcrayfish", marbledcrayfish)

ncTable <- rbind(coconutcrab,kingcrab,bluekingcrab, lobster,whiteshrimp,isopod,amphipod,swimmingcrab,marbledcrayfish)
#fwrite(ncTable, file = "noncodingtable.txt", sep = "\t")
ncTable <- ncTable[ncTable$evalue<0.00000000000001,] #Power
ncTable$sample.id <- do.call(paste, c(ncTable, sep=":"))
ncTable$startm <- pmin(ncTable$start,ncTable$stop)
ncTable$stopm <- pmax(ncTable$start,ncTable$stop)
ncTable <- unique(ncTable, by=c("startm","stopm"))

ncTable2 <- ncTable[,c(4,5,6,10)]
ncTable2$startm <- pmin(ncTable2$start,ncTable2$stop)
ncTable2$stopm <- pmax(ncTable2$start,ncTable2$stop)
ncTable2 <- ncTable2[,c(1,5,6,4)]

#fwrite(ncTable2[ncTable2$sample.id %like% "coconutcrab",], file = "noncodingtable_coconutcrab.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "kingcrab",], file = "noncodingtable_kingcrab.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "bluekingcrab",], file = "noncodingtable_bluekingcrab.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "lobster",], file = "noncodingtable_lobster.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "mitten",], file = "noncodingtable_mitten.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "whiteshrimp",], file = "noncodingtable_whiteshrimp.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "marbledcrayfish",], file = "noncodingtable_marbledcrayfish.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "isopod",], file = "noncodingtable_isopod.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "amphipod",], file = "noncodingtable_amphipod.txt", sep = "\t", col.names = F)
#fwrite(ncTable2[ncTable2$sample.id %like% "swimmingcrab",], file = "noncodingtable_swimmingcrab.txt", sep = "\t", col.names = F)



library(plyr)
longform <- data.table(count(ncTable, c("organism","name")))
wideform <- dcast(longform, name ~ organism, value.var = "freq")
wideform[is.na(wideform),] <- 0
#dt4cafe <- wideform[wideform$amphipod>0&wideform$coconutcrab>0&wideform$isopod>0
#                    &wideform$kingcrab>0&wideform$bluekingcrab>0&wideform$lobster>0&wideform$marbledcrayfish>0
#                    &wideform$mittencrab>0&wideform$swimmingcrab>0&wideform$whiteshrimp>0,]
dt4cafe <- wideform
dt4cafe <- dt4cafe[dt4cafe$name!="tRNA",]
colnames(dt4cafe)[1] <- "Desc"
dt4cafe <- add_column(dt4cafe, "Family ID" = "test", .after = "Desc")
dt4cafe$`Family ID` <- gsub("test", NA, dt4cafe$`Family ID`)
dt4cafe[, "max"] <- apply(dt4cafe[, 3:11], 1, max)
dt4cafe[, "min"] <- apply(dt4cafe[, 3:11], 1, min)
dt4cafe[, "sum"] <- apply(dt4cafe[, 3:11], 1, sum)
dt4cafe$diff <- dt4cafe$max - dt4cafe$min
dt4cafe <- dt4cafe[dt4cafe$sum>1&dt4cafe$sum>99,c(1:11)]
dt4cafe$`Family ID` <- 1:nrow(dt4cafe)

#names(noncoding_results)[1] <- "Family ID"
#test <- merge(dt4cafe,noncoding_results,by="Family ID")

#wideform[is.na(wideform)] <- 0
#fwrite(wideform, file = "noncodingtable_grouped_deduplicated_fifthpower.txt", sep = "\t")
#fwrite(wideform, file = "uncanny_mitten_similarities.csv", sep = ",")
#fwrite(dt4cafe, file = "outliernoncoding.csv", sep = ",")

###SNIPPETS (RUN AFTER SCRIPT)

#Extract rfam sequences for three decapods (after running above r script)
#CMD: bedtools getfasta -nameOnly -fi coconutcrab_masked.fasta -bed noncodingtable_coconutcrab.txt -fo coconutcrab_rfam.fasta 
#CMD: gt sequniq -o coconutcrab_rfam_deduped.fasta coconutcrab_rfam.fasta 

#Extracting sequences with a substring
#CMD: cat coconutcrab_rfam_deduped.fasta | grep tRNA -A1 | grep -v "\--" > coconutcrab_rfam_deduped_trna.fasta

#Run trnascan
#trnascan-1.4 -o kingcrab_rfam_trnascan.txt kingcrab_rfam_deduped_trna.fasta &> trnascan.out

#Printing only first predicted tRNAscan-SE amino acid
#CMD: awk '/sequence name/{print; nr[NR+10]; next}; NR in nr' coconutcrab_rfam_trnascan.txt > rfam_trnascan_first_predicted.txt

#Group and count amino acids
#sort rfam_trnascan_first_predicted.txt | uniq -c | grep "tRNA predict" > kingcrab_aminoacids.txt

#SUMMARY
# bedtools getfasta -nameOnly -fi swimmingcrabgenome.fna -bed noncodingtable_swimming_crab.txt -fo swimming_crab_rfam.fasta
# gt sequniq -o swimming_crab_rfam_deduped.fasta swimming_crab_rfam.fasta 
# cat swimming_crab_rfam_deduped.fasta | grep tRNA -A1 | grep -v "\--" > swimming_crab_rfam_deduped_trna.fasta
# trnascan-1.4 -o swimming_crab_rfam_trnascan.txt swimming_crab_rfam_deduped_trna.fasta &> trnascan.out
# awk '/sequence name/{print; nr[NR+10]; next}; NR in nr' swimming_crab_rfam_trnascan.txt > swimming_crab_rfam_trnascan_first_predicted.txt
# sort swimming_crab_rfam_trnascan_first_predicted.txt | uniq -c | grep "tRNA predict" > swimming_crab_aminoacids.txt