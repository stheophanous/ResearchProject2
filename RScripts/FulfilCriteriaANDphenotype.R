rm(list = ls())

library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#ExACandDGV cobined results
ExACandDGV <- read.table('/Users/steliostheophanous/Desktop/Results/Misc/CombinedExACandDGV.txt', header = F, fill = T)
#OMIM
genemap <- data.frame(read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/omim/genemap.txt', sep = '\t', header = T, fill = T))
mim2gene <- data.frame(read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/omim/mim2gene.txt', sep = '\t', header = T, fill = T))
ensemblID <- as.data.frame(org.Hs.egENSEMBL)
OMIMid <- as.data.frame(org.Hs.egOMIM)

results <- subset(ensemblID, ensembl_id %in% ExACandDGV$V1)
finalresults <- subset(OMIMid, gene_id %in% results$gene_id)

OMIMresults <- subset(genemap, MIM.Number %in% finalresults$omim_id)
GenesWithPhenotypes <- OMIMresults[!(is.na(OMIMresults$Phenotypes) | OMIMresults$Phenotypes==""), ]

#Convert to EnsemblIDs
BackToEns <- subset(OMIMid, omim_id %in% GenesWithPhenotypes$MIM.Number)
BackToEns2 <- subset(ensemblID, gene_id %in% BackToEns$gene_id)

write(as.character(BackToEns2$ensembl_id), file = "/Users/steliostheophanous/Desktop/ExACandDGVwithPhenotype.txt", sep = " ")

# 5 genes found to fulfil ExAC + DGV + DECIPHER criteria AND have a phenotype in OMIM:
# ENSG00000122025 - FLT3 - Class III Receptor Tyrosine Kinase (regulates hematopoiesis) - Acute myeloid leukaemia; Acute lymphoblastic leukaemia
# ENSG00000169032 - MAP2K1 - Mitogen Activated Protein Kinase - Cardiofaciocutaneous syndrome 3
# ENSG00000122026 - RPL21 - Component of the 60S ribosomal subunit - Hypotrichosis 12
# ENSG00000124614 - RPS10 - Component of the 40S ribosomal subunit - Diamond-Blackfan anemia 9
# ENSG00000138326 - RPS24 - Component of the 40S ribosomal subunit - Diamond-blackfan anemia 3
