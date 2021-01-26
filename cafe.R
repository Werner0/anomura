library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/CAFE/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/CAFE/")}

HOXpvalues <- fread("results_HOX/Base_family_results.txt")
colnames(HOXpvalues)[1] <- "Family ID"
HOXinput <- fread("HOX4cafe.txt")
HOXsignificant <- merge(HOXpvalues[HOXpvalues$`Significant at 0.05`=="y"],HOXinput,by="Family ID")
HOXsignificant$group <- "HOX" 
HOXsignificant$type <- "final" 

HOXinitialpvalues <- fread("results_HOX_initial/Base_family_results.txt")
colnames(HOXinitialpvalues)[1] <- "Family ID"
HOXinitialinput <- fread("initialHOX4cafewithoutCAD.txt")
HOXinitialsignificant <- merge(HOXinitialpvalues[HOXinitialpvalues$`Significant at 0.05`=="y"],HOXinitialinput,by="Family ID")
HOXinitialsignificant$group <- "HOX" 
HOXinitialsignificant$type <- "initial" 

Aminoinitialpvalues <- fread("results_aminoacids_initial/Base_family_results.txt")
colnames(Aminoinitialpvalues)[1] <- "Family ID"
Aminoinitialinput <- fread("initialaminoacids4cafe.txt")
Aminoinitialsignificant <- merge(Aminoinitialpvalues[Aminoinitialpvalues$`Significant at 0.05`=="y"],Aminoinitialinput,by="Family ID")
Aminoinitialsignificant$group <- "Amino" 
Aminoinitialsignificant$type <- "initial" 

Aminopvalues <- fread("results_aminoacids/Base_family_results.txt")
colnames(Aminopvalues)[1] <- "Family ID"
Aminoinput <- fread("aminoacids4cafe.txt")
Aminosignificant <- merge(Aminopvalues[Aminopvalues$`Significant at 0.05`=="y"],Aminoinput,by="Family ID")
Aminosignificant$group <- "Amino" 
Aminosignificant$type <- "Final" 

noncodingpvalues <- fread("results_noncoding_initial/Base_family_results.txt")
colnames(noncodingpvalues)[1] <- "Family ID"
noncodinginput <- fread("initialnoncoding4cafe.txt")
noncodingsignificant <- merge(noncodingpvalues[noncodingpvalues$`Significant at 0.05`=="y"],noncodinginput,by="Family ID")
noncodingsignificant$group <- "Noncoding" 
noncodingsignificant$type <- "Initial" 

codingpvalues <- fread("results_coding_initial/Base_family_results.txt")
colnames(codingpvalues)[1] <- "Family ID"
codinginput <- fread("initialcoding4cafe.txt")
codingsignificant <- merge(codingpvalues[codingpvalues$`Significant at 0.05`=="y"],codinginput,by="Family ID")
codingsignificant$group <- "Coding" 
codingsignificant$type <- "Initial" 

combo <- rbind(HOXsignificant, Aminoinitialsignificant,codingsignificant, noncodingsignificant)
fwrite(combo, "ExpansionAnalysis.csv", sep= ",")
