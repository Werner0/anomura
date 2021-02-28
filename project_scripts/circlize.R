library(data.table)
library(circlize)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/circlize/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/circlize/")}

#mynames <- c("P. hawaiensis", "P. platypus", "B. latro", "A. vulgare","P. camtschaticus", "P. ornatus", "P. virginalis", "P. trituberculatus", "L. vannamei")
mynames <- c("A","B","C","D","E","F","G","H","I")

onetoone <- fread("OrthologuesStats_one-to-one.tsv")
onetoone <- onetoone[,2:10]
onetoone <- as.matrix(onetoone)
rownames(onetoone) <- mynames
colnames(onetoone) <- mynames
onetoone <- onetoone[c("C","E","B","H","F","G","I","A","D"),c("C","E","B","H","F","G","I","A","D")]
rownames(onetoone) <- mynames
colnames(onetoone) <- mynames
col_fun <- colorRamp2(range(onetoone), c("#FFC20A", "#0C7BDC"), transparency = 0.5)

onetomany <- fread("OrthologuesStats_one-to-many.tsv")
onetomany <- onetomany[,2:10]
onetomany <- as.matrix(onetomany)
rownames(onetomany) <- mynames
colnames(onetomany) <- mynames
onetomany <- onetomany[c("C","E","B","H","F","G","I","A","D"),c("C","E","B","H","F","G","I","A","D")]
rownames(onetomany) <- mynames
colnames(onetomany) <- mynames

manytoone <- fread("OrthologuesStats_many-to-one.tsv")
manytoone <- manytoone[,2:10]
manytoone <- as.matrix(manytoone)
rownames(manytoone) <- mynames
colnames(manytoone) <- mynames
manytoone <- manytoone[c("C","E","B","H","F","G","I","A","D"),c("C","E","B","H","F","G","I","A","D")]
rownames(manytoone) <- mynames
colnames(manytoone) <- mynames

manytomany <- fread("OrthologuesStats_many-to-many.tsv")
manytomany <- manytomany[,2:10]
manytomany <- as.matrix(manytomany)
rownames(manytomany) <- mynames
colnames(manytomany) <- mynames
manytomany <- manytomany[c("C","E","B","H","F","G","I","A","D"),c("C","E","B","H","F","G","I","A","D")]
rownames(manytomany) <- mynames
colnames(manytomany) <- mynames

par(mar = c(1, 1, 1, 1), mfrow = c(2, 2))
chordDiagram(onetoone, symmetric = TRUE, grid.col = "brown", col = col_fun, annotationTrack = c("grid"),annotationTrackHeight = mm_h(2))
title("one-to-one")
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim)+0.1, mean(ylim)+0.1, si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = FALSE, col = "white", cex = 0.7)
}
chordDiagram(onetomany, transparency = 0.5, directional = 1, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow",
             grid.col = "brown", col = col_fun, annotationTrack = c("grid"),annotationTrackHeight = mm_h(2))
title("one-to-many")
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim)+0.1, mean(ylim)+0.1, si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = FALSE, col = "white", cex = 0.7)
}
chordDiagram(manytoone, transparency = 0.5, directional = 1, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow",
             grid.col = "brown", col = col_fun, annotationTrack = c("grid"),annotationTrackHeight = mm_h(2))
title("many-to-one")
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim)+0.1, mean(ylim)+0.1, si, sector.index = si, track.index = 1,
              facing = "bending.inside", niceFacing = FALSE, col = "white", cex = 0.7)
}
chordDiagram(manytomany, transparency = 0.5, directional = 1, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow",
             grid.col = "brown", col = col_fun, annotationTrack = c("grid"),annotationTrackHeight = mm_h(2))
title("many-to-many")

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim)+0.1, mean(ylim)+0.1, si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = FALSE, col = "white", cex = 0.7)
}