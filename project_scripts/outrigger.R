library(data.table)
library(stringr)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/outrigger/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/outrigger/")}

coconutcrab <- fread("coconutcrab_outrigger.csv", header = T)
coconutcrab_case8 <- coconutcrab[coconutcrab$notes%like%"Case 8",]
coconutcrab_reduced <- coconutcrab_case8[coconutcrab_case8$psi>0.05,c("sample_id","splice_type", "psi")]
coconutcrab_summary <- dcast(coconutcrab_reduced, sample_id ~ splice_type)
coconutcrab_summary$"reads_known" <- c("9366627","5702143","7324992","14632835")
coconutcrab_summary$"reads_unknown" <- c("11517306","8388998","15901278","9728718")
coconutcrab_summary$"se/mxe" <- coconutcrab_summary$se/coconutcrab_summary$mxe
coconutcrab_summary$organism <- "coconutcrab"

kingcrab <- fread("kingcrab_outrigger_summary.csv", header = T)
kingcrab_case8 <- kingcrab[kingcrab$notes%like%"Case 8",]
kingcrab_reduced <- kingcrab_case8[kingcrab_case8$psi>0.05,c("sample_id","splice_type", "psi")]
kingcrab_summary <- dcast(kingcrab_reduced, sample_id ~ splice_type)
kingcrab_summary$"reads_known" <- c("637164","593730","1186256","1317524") 
kingcrab_summary$"reads_unknown" <- c("2179077","2601830","4241297","2765245") 
kingcrab_summary$"se/mxe" <- kingcrab_summary$se/kingcrab_summary$mxe
kingcrab_summary$organism <- "kingcrab"

lobster <- fread("lobster_outrigger_summary.csv", header = T)
lobster_case8 <- lobster[lobster$notes%like%"Case 8",]
lobster_reduced <- lobster_case8[lobster_case8$psi>0.05,c("sample_id","splice_type", "psi")]
lobster_summary <- dcast(lobster_reduced, sample_id ~ splice_type)
lobster_summary$"reads_known" <- c("1368587","64863","602814","3533192","7453861") #ES GILL HP MS OV
lobster_summary$"reads_unknown" <- c("4706513","1111914","4718801","8527274","9249683") #ES GILL HP MS OV
lobster_summary$"se/mxe" <- lobster_summary$se/lobster_summary$mxe
lobster_summary$organism <- "lobster"

whiteshrimp <- fread("whteshrimp_many_outrigger_summary.csv", header = T)
whiteshrimp_case8 <- whiteshrimp[whiteshrimp$notes%like%"Case 8",]
whiteshrimp_reduced <- whiteshrimp_case8[whiteshrimp_case8$psi>0.05,c("sample_id","splice_type", "psi")]
whiteshrimp_summary <- dcast(whiteshrimp_reduced, sample_id ~ splice_type)
whiteshrimp_summary$"reads_known" <- c("2827909","7615008","3234481","1102506","4321546","6929419","5641271","7472469") 
whiteshrimp_summary$"reads_unknown" <- c("2929776","4814582","4980251","1639083","5301309","4548243","3897890","5834585")
whiteshrimp_summary$"se/mxe" <- whiteshrimp_summary$se/whiteshrimp_summary$mxe
whiteshrimp_summary$organism <- "whiteshrimp"

# ES      MS     MS2    GILL   GILL2      HP     HP2       O 
# 2827909 6929419 5641271 7615008 3234481 1102506 4321546 7472469 

all_summaries <- rbind(kingcrab_summary,lobster_summary,coconutcrab_summary, whiteshrimp_summary)
#fwrite(all_summaries, file = "AS_with_spit.csv", sep = ",")

#BED files
# library(stringr)
coco_1 <- coconutcrab_case8[coconutcrab_case8$psi>0.05,c("event_id","sample_id","splice_type")]
coco_1$event_id <- gsub("\\+|\\-\\@.*","",coco_1$event_id)
coco_1$event_id <- gsub(".*Seq","Seq",coco_1$event_id)
coco_1$event_id <- substr(coco_1$event_id,1,nchar(coco_1$event_id)-1)
coco_one <- as.data.table(str_split_fixed(coco_1$event_id, "-", 2))
coco_two <- as.data.table(str_split_fixed(coco_one$V1, ":", 2))
coco_bed <- data.table(coco_two$V1,coco_two$V2,coco_one$V2,coco_1$sample_id,coco_1$splice_type)
colnames(coco_bed) <- c("contig", "start", "stop", "tissue", "splice_type")
# 
# king_1 <- kingcrab_case8[kingcrab_case8$psi>0.05,c("event_id","sample_id","splice_type")]
# king_1$event_id <- gsub("\\+|\\-\\@.*","",king_1$event_id)
# king_1$event_id <- gsub(".*Seq","Seq",king_1$event_id)
# king_1$event_id <- substr(king_1$event_id,1,nchar(king_1$event_id)-1)
# king_one <- as.data.table(str_split_fixed(king_1$event_id, "-", 2))
# king_two <- as.data.table(str_split_fixed(king_one$V1, ":", 2))
# king_bed <- data.table(king_two$V1,king_two$V2,king_one$V2,king_1$sample_id,king_1$splice_type)
# colnames(king_bed) <- c("contig", "start", "stop", "tissue", "splice_type")
# 
# lobster_1 <- lobster_case8[lobster_case8$psi>0.05,c("event_id","sample_id","splice_type")]
# lobster_1$event_id <- gsub("\\+|\\-\\@.*","",lobster_1$event_id)
# lobster_1$event_id <- gsub(".*Seq","Seq",lobster_1$event_id)
# lobster_1$event_id <- substr(lobster_1$event_id,1,nchar(lobster_1$event_id)-1)
# lobster_one <- as.data.table(str_split_fixed(lobster_1$event_id, "-", 2))
# lobster_two <- as.data.table(str_split_fixed(lobster_one$V1, ":", 2))
# lobster_bed <- data.table(lobster_two$V1,lobster_two$V2,lobster_one$V2,lobster_1$sample_id,lobster_1$splice_type)
# colnames(lobster_bed) <- c("contig", "start", "stop", "tissue", "splice_type")

whiteshrimp_1 <- whiteshrimp_case8[whiteshrimp_case8$psi>0.05,c("event_id","sample_id","splice_type")]
whiteshrimp_1$event_id <- gsub("\\+|\\-\\@.*","",whiteshrimp_1$event_id)
whiteshrimp_1$event_id <- gsub(".*LVANscaffold","LVANscaffold",whiteshrimp_1$event_id)
whiteshrimp_1$event_id <- substr(whiteshrimp_1$event_id,1,nchar(whiteshrimp_1$event_id)-1)
whiteshrimp_one <- as.data.table(str_split_fixed(whiteshrimp_1$event_id, "-", 2))
whiteshrimp_two <- as.data.table(str_split_fixed(whiteshrimp_one$V1, ":", 2))
whiteshrimp_bed <- data.table(whiteshrimp_two$V1,whiteshrimp_two$V2,whiteshrimp_one$V2,whiteshrimp_1$sample_id,whiteshrimp_1$splice_type)
colnames(whiteshrimp_bed) <- c("contig", "start", "stop", "tissue", "splice_type")

#fwrite(coco_bed, file = "coconutcrab_AS.bed", sep = "\t")
#fwrite(king_bed, file = "kingcrab_AS.bed", sep = "\t")
#fwrite(lobster_bed, file = "lobster_AS.bed", sep = "\t")
#fwrite(whiteshrimp_bed, file = "whiteshrimp_AS.bed", sep = "\t")
#Remember to sort in bash