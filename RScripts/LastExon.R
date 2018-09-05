rm(list=ls())
library(data.table)
library(stringr)
library(biomaRt)

#Read ExAC variants
allvariants <- fread("gunzip -c /Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/exac/ExAC.r0.3.1.sites.vep.canonical.table.gz")
allvariants2 <- allvariants[,c(1, 2, 61, 62, 63, 66, 68, 76)] #Only columns I'm interested in.
hcvariants <- subset(allvariants2, LoF == 'HC') #Only high confidence LoF variants
lcvariants <- subset(allvariants2, LoF == 'LC') #Only low confidence LoF variants
allLOFvariants <- rbind(hcvariants, lcvariants) #Combine
allLOFvariants <- allLOFvariants[complete.cases(allLOFvariants$EXON), ] #Remove entries with exon = NA

#Split exon column to 2: exon_num, exon_total
exoncolumns <- str_split_fixed(allLOFvariants$EXON, "/", 2)
exoncolumns <- as.data.frame(exoncolumns)
setnames(exoncolumns, c('exon_num','exon_total'))
exoncolumns <- as.matrix(as.data.frame(lapply(exoncolumns, as.character)))
#exoncolumns <- as.matrix(as.data.frame(lapply(exoncolumns, as.numeric)))
allLOFvariants <- cbind(allLOFvariants, exoncolumns)

#Sort by chromosome and position + CLEAN
allLOFvariants <- allLOFvariants[complete.cases(allLOFvariants$exon_num), ]
allLOFvariants <- allLOFvariants[complete.cases(allLOFvariants$exon_total), ]
allLOFvariants <- allLOFvariants[order(allLOFvariants[,1], allLOFvariants[,2]), ]

#Last Exon?
allLOFvariants$last_exon <- NA
allLOFvariants$last_exon <- as.character(allLOFvariants$exon_num) == as.character(allLOFvariants$exon_total)

#Satisfies criteria?
allLOFvariants$criteria <- (allLOFvariants$Consequence == "splice_acceptor_variant") | (allLOFvariants$Consequence == "splice_donor_variant") | (allLOFvariants$last_exon == T & allLOFvariants$Consequence == "frameshift_variant") | (allLOFvariants$last_exon == T & allLOFvariants$Consequence == "stop_gained")

#All genes minus genes that do not satisfy criteria
allgenes <- allLOFvariants$Gene
allgenes <- allgenes[!duplicated(allgenes)]
allFalse <- subset(allLOFvariants, allLOFvariants$criteria == F)
allFalse <- allFalse$Gene
allFalse <- allFalse[!duplicated(allFalse)]

#Only genes that satisfy criteria
finalgenelist <- setdiff(allgenes, allFalse)

#Double check
allTrue <- subset(allLOFvariants, allLOFvariants$criteria == T)
allTrue <- allTrue$Gene
allTrue <- allTrue[!duplicated(allTrue)]
doublecheck <- intersect(finalgenelist, allTrue) #3006 genes -> all of the genes in the final list satisfy criteria

#=============================================================

#Genes that have no phenotype in OMIM
#Upload OMIM files
genemap <- data.frame(read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/omim/genemap.txt', sep = '\t', header = T, fill = T))
mim2gene <- data.frame(read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/omim/mim2gene.txt', sep = '\t', header = T, fill = T))

#Find genes without phenotypes
#METHOD 1
match <- subset(mim2gene, Ensembl.Gene.ID %in% finalgenelist)
MIMlist <- match$MIM.Number
phenotypes <- subset(genemap, MIM.Number %in% MIMlist)
#Include genes whose 'Phenotypes' start with:
inclBraces <- subset(phenotypes, grepl("^\\{", phenotypes$Phenotypes))
inclBrackets <- subset(phenotypes, grepl("^\\(", phenotypes$Phenotypes))
inclQmark <- subset(phenotypes, grepl("^\\?", phenotypes$Phenotypes))
inclNull <- subset(phenotypes, Phenotypes == "")
#All genes without phenotypes
NoPhenotypes <- rbind(inclBraces, inclBrackets, inclQmark, inclNull)
NoPhenotypes <- NoPhenotypes[order(NoPhenotypes$Sort),]
#Back to Gene_ID
NoPhenotypeList <- as.character(NoPhenotypes$MIM.Number)
NoPhenotypesGeneID <- subset(mim2gene, MIM.Number %in% NoPhenotypeList)
ListNoPhenotypes <- as.character(NoPhenotypesGeneID$Ensembl.Gene.ID)

#Not in OMIM1
notinOMIM <- setdiff(finalgenelist, match$Ensembl.Gene.ID)
notinOMIM <- notinOMIM[!duplicated(notinOMIM)]

#Not in OMIM2
mimlist2 <- MIMlist
phenlist2 <- as.numeric(as.character(phenotypes$MIM.Number))
different <- setdiff(mimlist2, phenlist2)
#Back to Gene_ID
notinOMIM2GeneID <- subset(mim2gene, MIM.Number %in% different)
notinOMIM2 <- as.character(notinOMIM2GeneID$Ensembl.Gene.ID)


#Combine results
finalresults <- c(as.character(ListNoPhenotypes), as.character(notinOMIM), as.character(notinOMIM2))
finalresults <- finalresults[!duplicated(finalresults)]

#Save final list
write(as.character(finalresults), file = "/Users/steliostheophanous/Desktop/finalgenelistLASTEXON.txt", sep = " ")







#========================================================================================
#Satisfies criteria?
splice_acceptors <- subset(allLOFvariants, Consequence == "splice_acceptor_variant")
splice_donors <- subset(allLOFvariants, Consequence == "splice_donor_variant")
last_exon <- subset(allLOFvariants, last_exon == T)
lastexon_stop <- subset(last_exon, Consequence == "stop_gained")
lastexon_frameshift <- subset(last_exon, Consequence == "frameshift_variant")
allTrue <- rbind(splice_acceptors, splice_donors, lastexon_stop, lastexon_frameshift)
allTrue <- allTrue$Gene
allTrue <- allTrue[!duplicated(allTrue)]
#========================================================================================





