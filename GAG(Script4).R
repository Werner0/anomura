#PRE-SNIPPETS: awk -F '[[:space:]][[:space:]]+' '{print $1,";",$3}' lobster_hints_utr.stats > lobster_hints_utr.stats.parsed

library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/GAG/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/GAG/")}

coconutcrab_hints_utr <- fread("coconutcrab_hints_utr.stats.parsed", skip = 2, header = F)
coconutcrab_hints_utr[,"prediction":="coconutcrab_hints_utr"]
coconutcrab_hints_noutr <- fread("coconutcrab_hints_noutr.stats.parsed", skip = 2, header = F)
coconutcrab_hints_noutr[,"prediction":="coconutcrab_hints_noutr"]
coconutcrab_abinitio_utr <- fread("coconutcrab_abinitio_utr.stats.parsed", skip = 2, header = F)
coconutcrab_abinitio_utr[,"prediction":="coconutcrab_abinitio_utr"]
coconutcrab_abinitio_noutr <- fread("coconutcrab_abinitio_noutr.stats.parsed", skip = 2, header = F)
coconutcrab_abinitio_noutr[,"prediction":="coconutcrab_abinitio_noutr"]

lobster_hints_utr <- fread("lobster_no_hp.stats.parsed", skip = 2, header = F)
lobster_hints_utr[,"prediction":="lobster_hints_utr"]
lobster_hints_noutr <- fread("lobster_noHP_hints.stats.parsed", skip = 2, header = F)
lobster_hints_noutr[,"prediction":="lobster_hints_noutr"]
lobster_abinitio_utr <- fread("lobster_noHP_abinitio_utr.stats.parsed", skip = 2, header = F)
lobster_abinitio_utr[,"prediction":="lobster_abinitio_utr"]
lobster_abinitio_noutr <- fread("lobster_noHP_abinitio.stats.parsed", skip = 2, header = F)
lobster_abinitio_noutr[,"prediction":="lobster_abinitio_noutr"]
#lobster_no_hp <- fread("lobster_no_hp.stats.parsed", skip = 2, header = F)
#lobster_no_hp[,"prediction":="lobster_no_hp"]

kingcrab_hints_utr <- fread("kingcrab_hints_utr.stats.parsed", skip = 2, header = F)
kingcrab_hints_utr[,"prediction":="kingcrab_hints_utr"]
kingcrab_hints_noutr <- fread("kingcrab_hints_noutr.stats.parsed", skip = 2, header = F)
kingcrab_hints_noutr[,"prediction":="kingcrab_hints_noutr"]
kingcrab_abinitio_utr <- fread("kingcrab_abinitio_utr.stats.parsed", skip = 2, header = F)
kingcrab_abinitio_utr[,"prediction":="kingcrab_abinitio_utr"]
kingcrab_abinitio_noutr <- fread("kingcrab_abinitio_noutr.stats.parsed", skip = 2, header = F)
kingcrab_abinitio_noutr[,"prediction":="kingcrab_abinitio_noutr"]

prediction_table <- rbind(coconutcrab_abinitio_noutr,
                          kingcrab_abinitio_noutr,
                          lobster_abinitio_noutr,
                          coconutcrab_abinitio_utr,
                          kingcrab_abinitio_utr,
                          lobster_abinitio_utr,
                          coconutcrab_hints_noutr,
                          kingcrab_hints_noutr,
                          lobster_hints_noutr,
                          coconutcrab_hints_utr,
                          kingcrab_hints_utr,
                          lobster_hints_utr)
                          #lobster_no_hp)
prediction_table <- transform(prediction_table, V2 = as.character(V2))
colnames(prediction_table) <- c("unit","number","prediction")

library(plyr)
wideform <- dcast(prediction_table, unit ~ prediction, value.var = "number", fun.aggregate = toString)
#fwrite(wideform, file = "gene_predictions_noHP.csv", sep = ",")
