rm(list = ls())

LOFresults <- read.delim('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Results/genelistLOF.txt', sep = " ", header = F)
CNVresults <- read.delim('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Results/genelistHomosapiens.txt',sep = " ", header = F)

LOFresults <- LOFresults[!duplicated(LOFresults),]
CNVresults <- CNVresults[!duplicated(CNVresults),]
combinedResults <- intersect(LOFresults, CNVresults)

#Add genes NOT in OMIM
notinOMIM <- read.csv('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Results/notinOMIM.txt', sep = " ", header = F)$V1
notinOMIM2 <- read.csv('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Results/NEWnotinOMIM.txt', sep = " ", header = F)$V1
notinOMIMall <- c(as.character(notinOMIM), as.character(notinOMIM2))
notinOMIMhomosapiens <- read.csv('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Results/notinOMIMhomosapiens.txt', sep = " ", header = F)$V1
notinOMIMall <- notinOMIMall[!duplicated(notinOMIMall)]
notinOMIMhomosapiens <- notinOMIMhomosapiens[!duplicated(notinOMIMhomosapiens)]
combinedNotInOMIM <- intersect(notinOMIMall, notinOMIMhomosapiens)

#Save results
FinalResults <- c(combinedResults, as.character(combinedNotInOMIM))
write(FinalResults, file = "/Users/steliostheophanous/Desktop/FinalResults.txt", sep = " ")

#write(combinedNotInOMIM, file = "/Users/steliostheophanous/Desktop/FINALnotinOMIM.txt", sep = " ")
#LOFresults <- read.delim('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Results/genelistLOF.txt', sep = " ", header = F)
#notinOMIM <- read.csv('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Results/notinOMIM.txt', sep = " ", header = F)$V1


# Method 2 - Same results
AllLoF <- c(as.character(LOFresults$V1), as.character(notinOMIMall))
AllCNV <- c(as.character(CNVresults$V1), as.character(notinOMIMhomosapiens))
AllLoF <- AllLoF[!duplicated(AllLoF)]
AllCNV <- AllCNV[!duplicated(AllCNV)]
combinedResults <- intersect(AllLoF, AllCNV)
