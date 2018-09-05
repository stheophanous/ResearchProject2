rm(list = ls())
library(biomaRt)
library(stringr)


#Exac (Only TranscriptID + GeneName can be used for comparison)
exac <- read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/exac/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt', header = T)
noLOF <- subset(exac, n_lof == 0)
#highPLI <- subset(exac, pLI >= 0.90)
#noLOFhighPLI <- subset(noLOF, pLI >= 0.90)

#OMIM
genemap <- data.frame(read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/omim/genemap.txt', sep = '\t', header = T, fill = T))
mim2gene <- data.frame(read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/omim/mim2gene.txt', sep = '\t', header = T, fill = T))

#*****Convert Transcript IDs*****#
#METHOD 1
#Using BiomaRt + Transcript ID
transcriptlist <- as.character(noLOF$transcript)
newtranscriptlist <- gsub("\\..*","",transcriptlist)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
convertedList <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),
                 filters = "ensembl_transcript_id", values = newtranscriptlist,
                 mart = mart)
match <- subset(mim2gene, Approved.Gene.Symbol %in% convertedList$external_gene_name)
match2 <- subset(mim2gene, Ensembl.Gene.ID %in% convertedList$ensembl_gene_id)
#Direct match of gene names
genenamelist <- as.character(noLOF$gene)
match3 <- subset(mim2gene, Approved.Gene.Symbol %in% genenamelist)

#Combine all matches and remove duplicates
allmatches <- rbind(match, match2, match3)
allmatchesNoDupl <- allmatches[!duplicated(allmatches), ]

#Find genes without phenotypes
MIMlist <- allmatchesNoDupl$MIM.Number  #MIMlist (from MIM2gene dataset): 857 genes
phenotypes <- subset(genemap, MIM.Number %in% MIMlist) #Finds 373 matches only??
#Include genes whose 'Phenotypes' start with:
inclBraces <- subset(phenotypes, grepl("^\\{", phenotypes$Phenotypes))
inclBrackets <- subset(phenotypes, grepl("^\\(", phenotypes$Phenotypes))
inclQmark <- subset(phenotypes, grepl("^\\?", phenotypes$Phenotypes))
inclNull <- subset(phenotypes, Phenotypes == "")
#All genes without phenotypes
NoPhenotypes <- rbind(inclBraces, inclBrackets, inclQmark, inclNull)
NoPhenotypes <- NoPhenotypes[order(NoPhenotypes$Sort),]


#METHOD 2
#Direct conversion from Exac to OMIM (matching gene names)
directconversion <- genemap[grep(paste("\\b",genenamelist,"\\b", sep = "", collapse="|"), genemap$Gene.Symbols), ]
inclBraces2 <- subset(directconversion, grepl("^\\{", directconversion$Phenotypes))
inclBrackets2 <- subset(directconversion, grepl("^\\(", directconversion$Phenotypes))
inclQmark2 <- subset(directconversion, grepl("^\\?", directconversion$Phenotypes))
inclNull2 <- subset(directconversion, Phenotypes == "")
NoPhenotypes2 <- rbind(inclBraces2, inclBrackets2, inclQmark2, inclNull2)
NoPhenotypes2 <- NoPhenotypes2[order(NoPhenotypes2$Sort),]

#Combine the 2 methods
combNoPhen <- rbind(NoPhenotypes, NoPhenotypes2)
combNoPhenNoDupl <- combNoPhen[!duplicated(combNoPhen), ]

#Reverse - These are NOT found in the OMIM database (No Phenotype)
reverse <- noLOF[!grepl(paste(allmatchesNoDupl$Approved.Gene.Symbol, collapse="|"), noLOF$gene), ]
transcriptlist2 <- as.character(reverse$transcript)
newtranscriptlist <- gsub("\\..*","",transcriptlist2)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
convertedList2 <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),
                       filters = "ensembl_transcript_id", values = newtranscriptlist,
                       mart = mart)
NoPhenotypeList2 <- as.character(convertedList2$ensembl_gene_id)
write(as.character(NoPhenotypeList2), file = "/Users/steliostheophanous/Desktop/notinOMIM.txt", sep = " ")


#Pathway analysis
NoPhenotypeList <- combNoPhenNoDupl$MIM.Number
back2gene <- subset(mim2gene, MIM.Number %in% NoPhenotypeList)
pathwayanalysis <- c(as.character(back2gene$Ensembl.Gene.ID), NoPhenotypeList2)
write(pathwayanalysis, file = "/Users/steliostheophanous/Desktop/genelistLOF.txt", sep = " ")


#Reverse2 - These are NOT found in the OMIM database (No Phenotype)
mimlist2 <- MIMlist
phenlist2 <- as.numeric(as.character(phenotypes$MIM.Number))
different <- setdiff(mimlist2, phenlist2)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
convertedList3 <- getBM(attributes = c("mim_gene_accession","ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),
                        filters = "mim_gene_accession", values = different,
                        mart = mart)
convertedList3 <- convertedList3[!duplicated(convertedList3), ]
newList3 <- as.character(convertedList3$ensembl_gene_id)
newList3 <- newList3[!duplicated(newList3)]
write(as.character(newList3), file = "/Users/steliostheophanous/Desktop/NEWnotinOMIM.txt", sep = " ")




