#SNIPPETS (USE BEFORE SCRIPT)
#awk -F\# '$1!="" { print $0 ;} ' lobster.emapper.annotations > lobster.emapper_nocomments.annotations

library(tibble)
library(data.table)

if (!getwd()=="/Users/wernerveldsman/Desktop/Rsources/Hox/") {
  setwd("/Users/wernerveldsman/Desktop/Rsources/Hox/")}

library(plyr)

lobster <- fread("lobster_noHP.emapper_nocomments.annotations")
redkingcrab <- fread("kingcrab.emapper_nocomments.annotations")
coconutcrab <- fread("coconutcrab.emapper_nocomments.annotations")
bluekingcrab <- fread("bluekingcrab.emapper_nocomments.annotations")
whiteshrimp <- fread("whiteshrimp.emapper_nocomments.annotations")
marbled <- fread("marbled.emapper_nocomments.annotations")
isopod <- fread("isopod.emapper_nocomments.annotations")
amphipod <- fread("amphipod.emapper_nocomments.annotations")
swimmingcrab <- fread("swimmingcrab_new.emapper_nocomments.annotations")


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
lobster[,("organism") := "lobster"]
colnames(redkingcrab) <- columnnames
redkingcrab[,("organism") := "redkingcrab"]
colnames(coconutcrab) <- columnnames
coconutcrab[,("organism") := "coconutcrab"]
colnames(whiteshrimp) <- columnnames
whiteshrimp[,("organism") := "whiteshrimp"]
colnames(marbled) <- columnnames
marbled[,("organism") := "marbled"]
colnames(bluekingcrab) <- columnnames
bluekingcrab[,("organism") := "bluekingcrab"]
colnames(swimmingcrab) <- columnnames
swimmingcrab[,("organism") := "swimmingcrab"]
colnames(isopod) <- columnnames
isopod[,("organism") := "isopod"]
colnames(amphipod) <- columnnames
amphipod[,("organism") := "amphipod"]

alltables <- rbind(lobster,redkingcrab,bluekingcrab,coconutcrab,whiteshrimp,marbled,swimmingcrab,isopod,amphipod)
alltables_reduced <- alltables[,c(1,5:6,18,21:23)]
alltables_reduced$Preferred_name <- toupper(alltables_reduced$Preferred_name)
#fwrite(alltables_reduced, "alltables_reduced.txt",sep ="\t")
test <- fread("ALLhoxGenes_2.txt")

hoxnames <- c("abd-A", "abdA", "iab-2", "iab-4", "Hab", "iab-3", "Abd-B", "AbdB", "iab-7", "Mcp", "iab-6", "Fab-7","Antp" , "ANT-C", "Hu", "Scx", "DmAntp", "Antennapaedia", "bcd", "btn", "cad","38E.19", "cd", "Dfd", "DmDfd", "eve","E(eve)", "l(2)46Ce", "ftz", "Dm-Ftz", "ind", "lab", "F90-2", "F24", "pb", "l(3)04498", "proboscipedia", "ro", "rough", "Scr", "Msc", "Multiple sex comb", "Ubx","bx", "pbx", "Cbx", "bxd", "bithorax", "zen", "z1", "zerknult", "zen2", "z2", "Complementation group F")
hoxnames2 <- c("ANT-C", "Hu", "Scx", "DmAntp", "Antennapaedia","DmDfd","F90-2", "F24","l(3)04498", "proposcipedia","Msc", "Multiple sex comb","abdA", "iab-2", "iab-4", "Hab", "iab-3","AbdB", "iab-7", "Mcp", "iab-6", "Fab-7","bx", "pbx", "Cbx", "bxd", "bithorax")
hoxnames3 <- c("CDX1","CDX2","CDX4","EN1","EN2","EVX1","EVX2","GBX1","GBX2","GSX1","GSX2","HOXA1","HOXB1","HOXD1","HOXA2","HOXB2","HOXA3","HOXB3","HOXD3","HOXA4","HOXB4","HOXC4","HOXD4","HOXA5","HOXB5","HOXC5","HOXA6","HOXB6","HOXC6","HOXA7","HOXB7","HOXB8","HOXC8","HOXD8","HOXA9","HOXB9","HOXC9","HOXD9","HOXA10","HOXC10","HOXD10","HOXA11","HOXC11","HOXD11","HOXC12","HOXD12","HOXA13","HOXB13","HOXC13","HOXD13","MNX1","MEOX1","MEOX2","PDX1","BARHL1","BARHL2","BARX1","BARX2","BSX","DBX1","DBX2","DLX1","DLX2","DLX3","DLX4","DLX5","DLX6","EMX1","EMX2","HHEX","HLX","LBX1","LBX2","MSX1","MSX2","NANOG","NKX1.1","NKX1.2","NKX2.1","NKX2.4","NKX2.2","NKX2.8","NKX3.1","NKX3.2","NKX2.3","NKX2.5","NKX2.6","HMX1","HMX2","HMX3","NKX6.1","NKX6.2","NKX6.3","NOTO","TLX1","TLX2","TLX3","VAX1","VAX2","VENTX","ALX1","ALX3","ALX4","ARGFX","ARX","DMBX1","DPRX","DRGX","DUXAI","DUXAII","DUXBI","DUXBII","ESX1","GSC","GSC2","HESX1","HOPX","ISX","MIXL","NOBOX","OTP","OTX1","OTX2","CRX","PAX3","PAX7","PAX4","PAX6","PHOX2A","PHOX2B","PITX1","PITX2","PITX3","PROP1","PRRX1","PRRX2","RAX","RAX2","RHOXF1","RHOXF2","RHOXF2B","SEBOX","SHOX","SHOX2","TPRX1","TPRXL","UNCX","VSX1","VSX2","LEUTX","ISL1","ISL2","LHX1","LHX5","LHX2","LHX9","LHX3","LHX4","LHX6","LHX8","LMX1A","LMX1B","HDX","POU1F1","POU2F1","POU2F2","POU2F3","POU3F1","POU3F2","POU3F3","POU3F4","POU4F1","POU4F2","POU4F3","POU5F1","POU5F2","POU6F1","POU6F2","HMBOX1","HNF1A","HNF1B","SIX1","SIX2","SIX3","SIX6","SIX4","SIX5","IRX1","IRX2","IRX3","IRX4","IRX5","IRX6","MEIS1","MEIS2","MEIS3","MKX","PBX1","PBX2","PBX3","PBX4","PKNOX1","PKNOX2","TGIF1","TGIF2","TGIF2LX","TGIF2LY","CUX1","CUX2","ONECUT1","ONECUT2","ONECUT3","SATB1","SATB2","PROX1","PROX2","ADNP","ADNP2","TSHZ1","TSHZ2","TSHZ3","ZEB1","ZEB2","ZFHX2I","ZFHX2II","ZFHX2III","ZFHX3I","ZFHX3II","ZFHX3III","ZFHX3IV","ZFHX4I","ZFHX4II","ZFHX4III","ZFHX4IV","ZHX1I","ZHX1II","ZHX1III","ZHX1IV","ZHX1V","ZHX2I","ZHX2II","ZHX2III","ZHX2IV","ZHX3I","ZHX3II","ZHX3III","ZHX3IV","ZHX3V","HOMEZI","HOMEZII","HOMEZIII","CERS2","CERS3","CERS4","CERS5","CERS6","AbdA","AbdB","achi","acj6","al","Antp","ap","ara","Awh","BH1","BH2","bap","bsh","btn","C15","cad","caup","CG11617","CG12361","CG32105","CG32532","CG33980","CG4136","CG4328","CG7056","CG9876","ct","Dfd","Dll","Dr","E5","ems","en","eve","exd","exex","ey","eyg","ftz","gsb","gsbn","Gsc","H2.0","HGTX","Hmx","hth","ind","inv","IP09201","IP17602","lab","Lag1","lbe","lbl","Lim1","Lim3","mirr","nub","oc","OdsH","onecut","Optix","otp","pb","pdm2","pdm3","PHDP","Pph13","prd","pros","Ptx1","Rx","Scr","scro","Six4","slou","so","toe","toy","tup","Ubx","unc4","unpg","vis","vnd","vvl","zfh1","zfh2I","zfh2II","zfh2III","Xlox","Prep","POU1","HNF","CART1","manacle","Hox3b","NK4")
hoxnames <- toupper(hoxnames)
hoxnames2 <- toupper(hoxnames2)
hoxnames3 <- toupper(hoxnames3)

fullhoxnames <- c("HsCDX1","HsCDX2","HsCDX4","HsEN1","HsEN2","HsEVX1","HsEVX2","HsGBX1","HsGBX2","HsGSX1","HsGSX2","HsHOXA1","HsHOXB1","HsHOXD1","HsHOXA2","HsHOXB2","HsHOXA3","HsHOXB3","HsHOXD3","HsHOXA4","HsHOXB4","HsHOXC4","HsHOXD4","HsHOXA5","HsHOXB5","HsHOXC5","HsHOXA6","HsHOXB6","HsHOXC6","HsHOXA7","HsHOXB7","HsHOXB8","HsHOXC8","HsHOXD8","HsHOXA9","HsHOXB9","HsHOXC9","HsHOXD9","HsHOXA10","HsHOXC10","HsHOXD10","HsHOXA11","HsHOXC11","HsHOXD11","HsHOXC12","HsHOXD12","HsHOXA13","HsHOXB13","HsHOXC13","HsHOXD13","HsMNX1","HsMEOX1","HsMEOX2","HsPDX1","HsBARHL1","HsBARHL2","HsBARX1","HsBARX2","HsBSX","HsDBX1","HsDBX2","HsDLX1","HsDLX2","HsDLX3","HsDLX4","HsDLX5","HsDLX6","HsEMX1","HsEMX2","HsHHEX","HsHLX","HsLBX1","HsLBX2","HsMSX1","HsMSX2","HsNANOG","HsNKX1.1","HsNKX1.2","HsNKX2.1","HsNKX2.4","HsNKX2.2","HsNKX2.8","HsNKX3.1","HsNKX3.2","HsNKX2.3","HsNKX2.5","HsNKX2.6","HsHMX1","HsHMX2","HsHMX3","HsNKX6.1","HsNKX6.2","HsNKX6.3","HsNOTO","HsTLX1","HsTLX2","HsTLX3","HsVAX1","HsVAX2","HsVENTX","HsALX1","HsALX3","HsALX4","HsARGFX","HsARX","HsDMBX1","HsDPRX","HsDRGX","HsDUXAI","HsDUXAII","HsDUXBI","HsDUXBII","HsESX1","HsGSC","HsGSC2","HsHESX1","HsHOPX","HsISX","HsMIXL","HsNOBOX","HsOTP","HsOTX1","HsOTX2","HsCRX","HsPAX3","HsPAX7","HsPAX4","HsPAX6","HsPHOX2A","HsPHOX2B","HsPITX1","HsPITX2","HsPITX3","HsPROP1","HsPRRX1","HsPRRX2","HsRAX","HsRAX2","HsRHOXF1","HsRHOXF2","HsRHOXF2B","HsSEBOX","HsSHOX","HsSHOX2","HsTPRX1","HsTPRXL","HsUNCX","HsVSX1","HsVSX2","HsLEUTX","HsISL1","HsISL2","HsLHX1","HsLHX5","HsLHX2","HsLHX9","HsLHX3","HsLHX4","HsLHX6","HsLHX8","HsLMX1A","HsLMX1B","HsHDX","HsPOU1F1","HsPOU2F1","HsPOU2F2","HsPOU2F3","HsPOU3F1","HsPOU3F2","HsPOU3F3","HsPOU3F4","HsPOU4F1","HsPOU4F2","HsPOU4F3","HsPOU5F1","HsPOU5F2","HsPOU6F1","HsPOU6F2","HsHMBOX1","HsHNF1A","HsHNF1B","HsSIX1","HsSIX2","HsSIX3","HsSIX6","HsSIX4","HsSIX5","HsIRX1","HsIRX2","HsIRX3","HsIRX4","HsIRX5","HsIRX6","HsMEIS1","HsMEIS2","HsMEIS3","HsMKX","HsPBX1","HsPBX2","HsPBX3","HsPBX4","HsPKNOX1","HsPKNOX2","HsTGIF1","HsTGIF2","HsTGIF2LX","HsTGIF2LY","HsCUX1","HsCUX2","HsONECUT1","HsONECUT2","HsONECUT3","HsSATB1","HsSATB2","HsPROX1","HsPROX2","HsADNP","HsADNP2","HsTSHZ1","HsTSHZ2","HsTSHZ3","HsZEB1","HsZEB2","HsZFHX2I","HsZFHX2II","HsZFHX2III","HsZFHX3I","HsZFHX3II","HsZFHX3III","HsZFHX3IV","HsZFHX4I","HsZFHX4II","HsZFHX4III","HsZFHX4IV","HsZHX1I","HsZHX1II","HsZHX1III","HsZHX1IV","HsZHX1V","HsZHX2I","HsZHX2II","HsZHX2III","HsZHX2IV","HsZHX3I","HsZHX3II","HsZHX3III","HsZHX3IV","HsZHX3V","HsHOMEZI","HsHOMEZII","HsHOMEZIII","HsCERS2","HsCERS3","HsCERS4","HsCERS5","HsCERS6","DmAbdA","DmAbdB","Dmachi","Dmacj6","Dmal","DmAntp","Dmap","Dmara","DmAwh","DmBH1","DmBH2","Dmbap","Dmbsh","Dmbtn","DmC15","Dmcad","Dmcaup","DmCG11617","DmCG12361","DmCG32105","DmCG32532","DmCG33980","DmCG4136","DmCG4328","DmCG7056","DmCG9876","Dmct","DmDfd","DmDll","DmDr","DmE5","Dmems","Dmen","Dmeve","Dmexd","Dmexex","Dmey","Dmeyg","Dmftz","Dmgsb","Dmgsbn","DmGsc","DmH2.0","DmHGTX","DmHmx","Dmhth","Dmind","Dminv","DmIP09201","DmIP17602","Dmlab","DmLag1","Dmlbe","Dmlbl","DmLim1","DmLim3","Dmmirr","Dmnub","Dmoc","DmOdsH","Dmonecut","DmOptix","Dmotp","Dmpb","Dmpdm2","Dmpdm3","DmPHDP","DmPph13","Dmprd","Dmpros","DmPtx1","DmRx","DmScr","Dmscro","DmSix4","Dmslou","Dmso","Dmtoe","Dmtoy","Dmtup","DmUbx","Dmunc4","Dmunpg","Dmvis","Dmvnd","Dmvvl","Dmzfh1","Dmzfh2I","Dmzfh2II","Dmzfh2III","PsXlox","AmPrep","NvPOU1","NvHNF","NvCART1","Hvmanacle","SmHox3b","PdNK4")
fullhoxnames <- substr(fullhoxnames,1,2)
dt <- data.table(hoxnames3,fullhoxnames)
colnames(dt) <- c("Preferred_name","Source")

hoxnames <- toupper(hoxnames)
hoxnames2 <- toupper(hoxnames2)
hoxnames3 <- toupper(hoxnames3)

HOXfound <- alltables_reduced[alltables_reduced$Preferred_name %in% hoxnames|alltables_reduced$Preferred_name %in% hoxnames2|alltables_reduced$Preferred_name %in% hoxnames3,c(3,7)]
HOXfound$amount <- 1
wideform <- dcast(HOXfound, Preferred_name ~ organism, value.var = "amount")
wideform_2 <- merge(wideform,dt,by="Preferred_name", all = T)
wideform_3 <- wideform_2[!is.na(wideform_2$amphipod),]
wideform_3$Source[1:2] <- "Dm" #Tweaking
wideform_4 <- wideform_3[wideform_3$Source=="Dm",c(1:10)]
colnames(wideform_4)[1] <- "Desc"
wideform_4 <- add_column(wideform_4, "Family ID" = "test", .after = "Desc")
wideform_4$`Family ID` <- gsub("test", NA, wideform_4$`Family ID`)
wideform_4$`Family ID` <- 1:nrow(wideform_4)
#wideform_4 <- wideform_4[wideform_4$amphipod>0&wideform_4$coconutcrab>0&wideform_4$isopod>0
#         &wideform_4$kingcrab>0&wideform_4$bluekingcrab>0&wideform_4$lobster>0&wideform_4$marbled>0
#         &wideform_4$swimmingcrab>0&wideform_4$whiteshrimp>0,]

#foundnames <- c("Abdominal-A","Abdominal-B","Antennapedia","Caudal","Deformed","Even skipped","Fushi tarazu","Labial","Sex combs reduced", "Ultrabithorax")
#section <- c("Bithorax", "Bithorax", "Antennapedia", "U", "Antennapedia", "U", "U", "Antennapedia","Antennapedia", "Bithorax")
#wideform$longname <- foundnames
#wideform$complex <- section

#fwrite(alltables_reduced, file="TABLE4HOX_new.csv", sep = ",")

#wideform_4 <- wideform_4[wideform_4$Desc!="CAD",]
#fwrite(wideform_4, file ="initialHOX4cafewithoutCAD.txt", sep="\t")

#fwrite(wideform_4, file ="HOX4cafe.txt", sep="\t")
