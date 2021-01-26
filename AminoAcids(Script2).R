library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/rfam_parsed") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/rfam_parsed")}

coconutcrab <- fread("coconutcrab_aminoacids.txt")
coconutcrab <- coconutcrab[,c(7,1,10)]
coconutcrab[, "organism" := "coconutcrab"]
colnames(coconutcrab) <- c("aminoacid", "amount", "anticodon", "organism")

kingcrab <- fread("kingcrab_aminoacids.txt")
kingcrab <- kingcrab[,c(7,1,10)]
kingcrab[, "organism" := "redkingcrab"]
colnames(kingcrab) <- c("aminoacid", "amount", "anticodon", "organism")

bluekingcrab <- fread("bluekingcrab_aminoacids.txt")
bluekingcrab <- bluekingcrab[,c(7,1,10)]
bluekingcrab[, "organism" := "bluekingcrab"]
colnames(bluekingcrab) <- c("aminoacid", "amount", "anticodon", "organism")

lobster <- fread("lobster_aminoacids.txt")
lobster <- lobster[,c(7,1,10)]
lobster[, "organism" := "lobster"]
colnames(lobster) <- c("aminoacid", "amount", "anticodon", "organism")

mitten <- fread("mitten_aminoacids.txt")
mitten <- mitten[,c(7,1,10)]
mitten[, "organism" := "mitten"]
colnames(mitten) <- c("aminoacid", "amount", "anticodon", "organism")

swimming <- fread("swimming_crab_new_aminoacids.txt")
swimming <- swimming[,c(7,1,10)]
swimming[, "organism" := "swimming"]
colnames(swimming) <- c("aminoacid", "amount", "anticodon", "organism")

marbled <- fread("marbled_aminoacids.txt")
marbled <- marbled[,c(7,1,10)]
marbled[, "organism" := "marbled"]
colnames(marbled) <- c("aminoacid", "amount", "anticodon", "organism")

isopod <- fread("isopod_aminoacids.txt")
isopod <- isopod[,c(7,1,10)]
isopod[, "organism" := "isopod"]
colnames(isopod) <- c("aminoacid", "amount", "anticodon", "organism")

amphipod <- fread("amphipod_aminoacids.txt")
amphipod <- amphipod[,c(7,1,10)]
amphipod[, "organism" := "amphipod"]
colnames(amphipod) <- c("aminoacid", "amount", "anticodon", "organism")

whiteshrimp <- fread("whiteshrimp_aminoacids.txt")
whiteshrimp <- whiteshrimp[,c(7,1,10)]
whiteshrimp[, "organism" := "whiteshrimp"]
colnames(whiteshrimp) <- c("aminoacid", "amount", "anticodon", "organism")

aminoacids <- rbind(coconutcrab,lobster,kingcrab,bluekingcrab,isopod,amphipod,swimming,marbled,whiteshrimp)
wideformdetail <- dcast(aminoacids, aminoacid+anticodon ~ organism, value.var = "amount")
wideformdetail[is.na(wideformdetail)] <- 0
wideform <- dcast(aminoacids[,c(1,2,4)], aminoacid ~ organism, value.var = "amount", fun.aggregate = sum)

#fwrite(wideform, file ="aminoacids.csv", sep = ",")
#fwrite(wideformdetail, file ="aminoacids.csv", sep = ",")

#For CAFE
library(tibble)
colnames(wideformdetail)[2] <- "Desc"
wideformdetail <- add_column(wideformdetail, "Family ID" = "test", .after = "Desc")
wideformdetail$`Family ID` <- gsub("test", NA, wideformdetail$`Family ID`)
wideformdetail[, "max"] <- apply(wideformdetail[, 4:12], 1, max)
wideformdetail[, "min"] <- apply(wideformdetail[, 4:12], 1, min)
wideformdetail[, "sum"] <- apply(wideformdetail[, 4:12], 1, sum)
wideformdetail$diff <- wideformdetail$max - wideformdetail$min
wideformreference <- wideformdetail
fwrite(wideformreference, file = "wideformreference.txt", sep = ",")
wideformdetail <- wideformdetail[wideformdetail$diff>99&wideformdetail$sum>0,c(1:12)]
wideformdetail$`Family ID` <- 1:nrow(wideformdetail)
#fwrite(wideformdetail, file = "outlieraminostandalone.csv", sep = ",")
