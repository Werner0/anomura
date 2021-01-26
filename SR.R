library(data.table)
library(Biostrings)
library(plyr)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/SR/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/SR/")}

coconutcrab <- readAAStringSet("coconutcrab.aa")
coconutcrabAA <- alphabetFrequency(coconutcrab)
coconutcrabAA <- coconutcrabAA[,c("S","R")]
coconutcrabAA <- as.data.table(cbind(names(coconutcrab),width(coconutcrab),coconutcrabAA))
colnames(coconutcrabAA) <- c("id", "width", "serine", "arginine")
coconutcrabAA[,"SR"] <- as.numeric(coconutcrabAA$serine)+as.numeric(coconutcrabAA$arginine)
coconutcrabAA <- coconutcrabAA[coconutcrabAA$SR>0,]
coconutcrabAA[,"SR.W"] <- as.numeric(coconutcrabAA$width)/coconutcrabAA$SR
coconutcrab_e <- fread("coconutcrab.emapper_nocomments.annotations")
coconutcrab_e <- coconutcrab_e[,c(1,6)]
colnames(coconutcrab_e) <- c("id", "symbol")
coconutcrab_m <- merge(coconutcrabAA,coconutcrab_e,by="id", all = TRUE)
#coconutcrab_m <- coconutcrab_m[coconutcrab_m$symbol!="",]
coconutcrab_m <- coconutcrab_m[order(coconutcrab_m$SR.W),]
coconutcrab_m$organism <- "coconutcrab"

kingcrab <- readAAStringSet("kingcrab.aa")
kingcrabAA <- alphabetFrequency(kingcrab)
kingcrabAA <- kingcrabAA[,c("S","R")]
kingcrabAA <- as.data.table(cbind(names(kingcrab),width(kingcrab),kingcrabAA))
colnames(kingcrabAA) <- c("id", "width", "serine", "arginine")
kingcrabAA[,"SR"] <- as.numeric(kingcrabAA$serine)+as.numeric(kingcrabAA$arginine)
kingcrabAA <- kingcrabAA[kingcrabAA$SR>0,]
kingcrabAA[,"SR.W"] <- as.numeric(kingcrabAA$width)/kingcrabAA$SR
kingcrab_e <- fread("kingcrab.emapper_nocomments.annotations")
kingcrab_e <- kingcrab_e[,c(1,6)]
colnames(kingcrab_e) <- c("id", "symbol")
kingcrab_m <- merge(kingcrabAA,kingcrab_e,by="id", all = TRUE)
#kingcrab_m <- kingcrab_m[kingcrab_m$symbol!="",]
kingcrab_m <- kingcrab_m[order(kingcrab_m$SR.W),]
kingcrab_m$organism <- "kingcrab"

lobster <- readAAStringSet("lobster.aa")
lobsterAA <- alphabetFrequency(lobster)
lobsterAA <- lobsterAA[,c("S","R")]
lobsterAA <- as.data.table(cbind(names(lobster),width(lobster),lobsterAA))
colnames(lobsterAA) <- c("id", "width", "serine", "arginine")
lobsterAA[,"SR"] <- as.numeric(lobsterAA$serine)+as.numeric(lobsterAA$arginine)
lobsterAA <- lobsterAA[lobsterAA$SR>0,]
lobsterAA[,"SR.W"] <- as.numeric(lobsterAA$width)/lobsterAA$SR
lobster_e <- fread("lobster.emapper_nocomments.annotations")
lobster_e <- lobster_e[,c(1,6)]
colnames(lobster_e) <- c("id", "symbol")
lobster_m <- merge(lobsterAA,lobster_e,by="id", all = TRUE)
#lobster_m <- lobster_m[lobster_m$symbol!="",]
lobster_m <- lobster_m[order(lobster_m$SR.W),]
lobster_m$organism <- "lobster"

bluekingcrab <- readAAStringSet("bluekingcrab.fa")
bluekingcrabAA <- alphabetFrequency(bluekingcrab)
bluekingcrabAA <- bluekingcrabAA[,c("S","R")]
bluekingcrabAA <- as.data.table(cbind(names(bluekingcrab),width(bluekingcrab),bluekingcrabAA))
colnames(bluekingcrabAA) <- c("id", "width", "serine", "arginine")
bluekingcrabAA[,"SR"] <- as.numeric(bluekingcrabAA$serine)+as.numeric(bluekingcrabAA$arginine)
bluekingcrabAA <- bluekingcrabAA[bluekingcrabAA$SR>0,]
bluekingcrabAA[,"SR.W"] <- as.numeric(bluekingcrabAA$width)/bluekingcrabAA$SR
bluekingcrab_e <- fread("bluekingcrab.emapper_nocomments.annotations")
bluekingcrab_e <- bluekingcrab_e[,c(1,6)]
colnames(bluekingcrab_e) <- c("id", "symbol")
bluekingcrab_m <- merge(bluekingcrabAA,bluekingcrab_e,by="id", all = TRUE)
#bluekingcrab_m <- bluekingcrab_m[bluekingcrab_m$symbol!="",]
bluekingcrab_m <- bluekingcrab_m[order(bluekingcrab_m$SR.W),]
bluekingcrab_m$organism <- "bluekingcrab"

isopod <- readAAStringSet("isopod.aa")
names(isopod) <- gsub("[[:space:]].*", "",names(isopod))
isopodAA <- alphabetFrequency(isopod)
isopodAA <- isopodAA[,c("S","R")]
isopodAA <- as.data.table(cbind(names(isopod),width(isopod),isopodAA))
colnames(isopodAA) <- c("id", "width", "serine", "arginine")
isopodAA[,"SR"] <- as.numeric(isopodAA$serine)+as.numeric(isopodAA$arginine)
isopodAA <- isopodAA[isopodAA$SR>0,]
isopodAA[,"SR.W"] <- as.numeric(isopodAA$width)/isopodAA$SR
isopod_e <- fread("isopod.emapper_nocomments.annotations")
isopod_e <- isopod_e[,c(1,6)]
colnames(isopod_e) <- c("id", "symbol")
isopod_m <- merge(isopodAA,isopod_e,by="id", all = TRUE)
#isopod_m <- isopod_m[isopod_m$symbol!="",]
isopod_m <- isopod_m[order(isopod_m$SR.W),]
isopod_m$organism <- "isopod"

amphipod <- readAAStringSet("amphipod.aa")
amphipodAA <- alphabetFrequency(amphipod)
amphipodAA <- amphipodAA[,c("S","R")]
amphipodAA <- as.data.table(cbind(names(amphipod),width(amphipod),amphipodAA))
colnames(amphipodAA) <- c("id", "width", "serine", "arginine")
amphipodAA[,"SR"] <- as.numeric(amphipodAA$serine)+as.numeric(amphipodAA$arginine)
amphipodAA <- amphipodAA[amphipodAA$SR>0,]
amphipodAA[,"SR.W"] <- as.numeric(amphipodAA$width)/amphipodAA$SR
amphipod_e <- fread("amphipod.emapper_nocomments.annotations")
amphipod_e <- amphipod_e[,c(1,6)]
colnames(amphipod_e) <- c("id", "symbol")
amphipod_e$id <- gsub("dottod", ".",amphipod_e$id)
amphipod_m <- merge(amphipodAA,amphipod_e,by="id", all = TRUE)
#amphipod_m <- amphipod_m[amphipod_m$symbol!="",]
amphipod_m <- amphipod_m[order(amphipod_m$SR.W),]
amphipod_m$organism <- "amphipod"

whiteshrimp <- readAAStringSet("lva.proteins.faa")
names(whiteshrimp) <- gsub("[[:space:]].*", "",names(whiteshrimp))
whiteshrimpAA <- alphabetFrequency(whiteshrimp)
whiteshrimpAA <- whiteshrimpAA[,c("S","R")]
whiteshrimpAA <- as.data.table(cbind(names(whiteshrimp),width(whiteshrimp),whiteshrimpAA))
colnames(whiteshrimpAA) <- c("id", "width", "serine", "arginine")
whiteshrimpAA[,"SR"] <- as.numeric(whiteshrimpAA$serine)+as.numeric(whiteshrimpAA$arginine)
whiteshrimpAA <- whiteshrimpAA[whiteshrimpAA$SR>0,]
whiteshrimpAA[,"SR.W"] <- as.numeric(whiteshrimpAA$width)/whiteshrimpAA$SR
whiteshrimp_e <- fread("whiteshrimp.emapper_nocomments.annotations")
whiteshrimp_e <- whiteshrimp_e[,c(1,6)]
colnames(whiteshrimp_e) <- c("id", "symbol")
whiteshrimp_m <- merge(whiteshrimpAA,whiteshrimp_e,by="id", all = TRUE)
#whiteshrimp_m <- whiteshrimp_m[whiteshrimp_m$symbol!="",]
whiteshrimp_m <- whiteshrimp_m[order(whiteshrimp_m$SR.W),]
whiteshrimp_m$organism <- "whiteshrimp"

swimmingcrab <- readAAStringSet("Portunus_trituberculatus_gene_protein.fasta")
names(swimmingcrab) <- gsub(" EVM.*", "",names(swimmingcrab))
swimmingcrabAA <- alphabetFrequency(swimmingcrab)
swimmingcrabAA <- swimmingcrabAA[,c("S","R")]
swimmingcrabAA <- as.data.table(cbind(names(swimmingcrab),width(swimmingcrab),swimmingcrabAA))
colnames(swimmingcrabAA) <- c("id", "width", "serine", "arginine")
swimmingcrabAA[,"SR"] <- as.numeric(swimmingcrabAA$serine)+as.numeric(swimmingcrabAA$arginine)
swimmingcrabAA <- swimmingcrabAA[swimmingcrabAA$SR>0,]
swimmingcrabAA[,"SR.W"] <- as.numeric(swimmingcrabAA$width)/swimmingcrabAA$SR
swimmingcrab_e <- fread("swimmingcrab_new.emapper_nocomments.annotations")
swimmingcrab_e <- swimmingcrab_e[,c(1,6)]
colnames(swimmingcrab_e) <- c("id", "symbol")
swimmingcrab_m <- merge(swimmingcrabAA,swimmingcrab_e,by="id", all = TRUE)
#swimmingcrab_m <- swimmingcrab_m[swimmingcrab_m$symbol!="",]
swimmingcrab_m <- swimmingcrab_m[order(swimmingcrab_m$SR.W),]
swimmingcrab_m$organism <- "swimmingcrab"

marbledcrayfish <- readAAStringSet("Pvir04.predicted_proteins.fasta")
names(marbledcrayfish) <- gsub("[[:space:]].*", "",names(marbledcrayfish))
marbledcrayfishAA <- alphabetFrequency(marbledcrayfish)
marbledcrayfishAA <- marbledcrayfishAA[,c("S","R")]
marbledcrayfishAA <- as.data.table(cbind(names(marbledcrayfish),width(marbledcrayfish),marbledcrayfishAA))
colnames(marbledcrayfishAA) <- c("id", "width", "serine", "arginine")
marbledcrayfishAA[,"SR"] <- as.numeric(marbledcrayfishAA$serine)+as.numeric(marbledcrayfishAA$arginine)
marbledcrayfishAA <- marbledcrayfishAA[marbledcrayfishAA$SR>0,]
marbledcrayfishAA[,"SR.W"] <- as.numeric(marbledcrayfishAA$width)/marbledcrayfishAA$SR
marbledcrayfish_e <- fread("marbled.emapper_nocomments.annotations")
marbledcrayfish_e <- marbledcrayfish_e[,c(1,6)]
colnames(marbledcrayfish_e) <- c("id", "symbol")
marbledcrayfish_m <- merge(marbledcrayfishAA,marbledcrayfish_e,by="id", all = TRUE)
#marbledcrayfish_m <- marbledcrayfish_m[marbledcrayfish_m$symbol!="",]
marbledcrayfish_m <- marbledcrayfish_m[order(marbledcrayfish_m$SR.W),]
marbledcrayfish_m$organism <- "marbledcrayfish"

alltables <- rbind(lobster_m,kingcrab_m,coconutcrab_m, bluekingcrab_m, isopod_m, amphipod_m, whiteshrimp_m, swimmingcrab_m, marbledcrayfish_m)
alltables <- alltables[order(alltables$SR.W),]
t <- alltables[alltables$SR.W<4,]
count(t[,c("organism")])
x <- count(t[,c("organism","symbol")])
y <- as.data.table(x[order(-x$freq),])
y[y$symbol %in% c("B52")]

#lookup <- fread("../bed/kingcrab2.bed")
#z <- lookup[lookup$V4 %like% "g48705.t2",]
#z <- z[order(z$V4),]
