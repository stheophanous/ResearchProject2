rm(list = ls())
library(bedr)
library(GenomicRanges)
library(biomaRt)

#-------------------------------------------------------------------------------------------------------------------------------
check.binary("bedtools")
check.binary("bedops")

#Dataset 1: All regions
DGV <- read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/DGV/GRCh38_hg38_variants_2016-08-31.txt', header = T, fill = T)
sortedDGV <- DGV[order(as.numeric(as.character(DGV[,2])),as.numeric(as.character(DGV[,3])) ),]
sortedDGV$region <- paste("chr", as.character(sortedDGV$chr),":", as.character(sortedDGV$start),"-", as.character(sortedDGV$end), sep = "")
allregions <- as.character(sortedDGV$region)
notneeded1 <- grep("_alt", sortedDGV$region, value=TRUE)
notneeded2 <- grep("_random", sortedDGV$region, value=TRUE)
notneeded3 <- grep("Un_KI", sortedDGV$region, value=TRUE)
notneeded4 <- grep("\\[", sortedDGV$region, value=TRUE)
notneeded5 <- grep("chr:-", sortedDGV$region, value=TRUE)
remove <- c(notneeded1, notneeded2, notneeded3, notneeded4, notneeded5)
allregions <- allregions[!allregions %in% remove]
notneeded6 <- subset(sortedDGV, as.numeric(as.character(sortedDGV$start)) == as.numeric(as.character(sortedDGV$end)))
allregions <- allregions[!allregions %in% notneeded6$region]

# Dataset 2: Loss and Loss+Gain regions
OnlyLoss <- subset(sortedDGV, variantsubtype == "loss")
GainAndLoss <- subset(sortedDGV, variantsubtype == "gain+loss")
Deletion <- subset(sortedDGV, variantsubtype == "deletion")
AllCNVs <- rbind(OnlyLoss, GainAndLoss, Deletion)
AllCNVs$genesize <- as.numeric(as.character(AllCNVs$end)) - as.numeric(as.character(AllCNVs$start)) #Add gene size column
AllCNVs <- subset(AllCNVs, genesize <= 5000000)
sortedAllCNVs <- AllCNVs[order(as.numeric(as.character(AllCNVs[,2])),as.numeric(as.character(AllCNVs[,3])) ),]
myregions <- as.character(sortedAllCNVs$region)
notneeded11 <- grep('_alt', myregions, value=TRUE)
notneeded12 <- c("chr4_GL000008v2_random:43-5510", "chr4_GL000008v2_random:72-2074", "chr4_GL000008v2_random:87-7500", "chr4_GL000008v2_random:2161-3080", "chr4_GL000008v2_random:122221-195558", "chr4_GL000008v2_random:140926-141998", "chr4_GL000008v2_random:141066-154438", "chr4_GL000008v2_random:143364-184071", "chr4_GL000008v2_random:156159-178479", "chr4_GL000008v2_random:2820-5652", "chr4_GL000008v2_random:4649-6851", "chr4_GL000008v2_random:4824-5892", "chr1_KI270706v1_random:414-18422", "chr1_KI270706v1_random:5933-9736", "chr15_KI270727v1_random:17624-19073", "chrUn_KI270742v1:30-10729", "chrUn_KI270742v1:104390-105003")
notneeded13 <- grep("_random", myregions, value=TRUE)
notneeded14 <- grep("Un_KI", myregions, value=TRUE)
myregions <- myregions[!myregions %in% notneeded11]
myregions <- myregions[!myregions %in% notneeded12]
myregions <- myregions[!myregions %in% notneeded13]
myregions <- myregions[!myregions %in% notneeded14]

# Check if valid
is.allregions.valid <- is.valid.region(allregions)
is.myregions.valid <- is.valid.region(myregions)

# Sort them
allregions.sort <- bedr.sort.region(allregions)
myregions.sort <- bedr.sort.region(myregions)
# Check if sorted
#is.allregions.sorted <- is.sorted.region(allregions)
#is.myregions.sorted <- is.sorted.region(myregions)

# Merge overlapping regions
allregions.merge <- bedr.merge.region(allregions.sort)
myregions.merge <- bedr.merge.region(myregions.sort)
write(myregions.merge, file = "/Users/steliostheophanous/Desktop/lossregions.txt", sep = "\t")


# Check if merged
#is.allregions.merged <- is.merged.region(allregions)
#is.myregions.merged <- is.merged.region(myregions)

# Subtract Dataset 2 from Dataset 1
THEregions <- bedr.subtract.region(allregions.merge, myregions.merge, remove.whole.feature = F)
write(THEregions, file = "/Users/steliostheophanous/Desktop/bedtoolsregions.txt", sep = " ")
#------------------------------------------------------------------------------------------------------------------------------------

#---ReadRegions---#
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
region.df <- read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Results/bedtoolsregionsNEW.txt')
region.df$V1 <- gsub("-", ":", region.df$V1)
#region.df$V1 <- gsub("chr", "", region.df$V1)
region.dfnew <- str_split_fixed(region.df$V1, ":", 3)
region.dfnew <- as.data.frame(region.dfnew)
colnames(region.dfnew) <- c("chr", "start", "end")
regiongr <- makeGRangesFromDataFrame(region.dfnew)
genes <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), regiongr)

ensemblID <- as.data.frame(org.Hs.egENSEMBL)
OMIMid <- as.data.frame(org.Hs.egOMIM)
results <- subset(ensemblID, gene_id %in% genes$gene_id)
finalresults <- subset(OMIMid, gene_id %in% results$gene_id)

#OMIM
genemap <- data.frame(read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/omim/genemap.txt', sep = '\t', header = T, fill = T))
mim2gene <- data.frame(read.table('/Users/steliostheophanous/Documents/Bioinformatics/Project2/Stelios/omim/mim2gene.txt', sep = '\t', header = T, fill = T))

#Find genes without phenotypes
MIMlist <- finalresults$omim_id  #MIMlist (from MIM2gene dataset): 857 genes
phenotypes <- subset(genemap, MIM.Number %in% MIMlist) #Finds 373 matches only??
#Include genes whose 'Phenotypes' start with:
inclBraces <- subset(phenotypes, grepl("^\\{", phenotypes$Phenotypes))
inclBrackets <- subset(phenotypes, grepl("^\\(", phenotypes$Phenotypes))
inclQmark <- subset(phenotypes, grepl("^\\?", phenotypes$Phenotypes))
inclNull <- subset(phenotypes, Phenotypes == "")
#All genes without phenotypes
NoPhenotypes <- rbind(inclBraces, inclBrackets, inclQmark, inclNull)
NoPhenotypes <- NoPhenotypes[order(NoPhenotypes$Sort),]

#Pathway analysis
NoPhenotypeList <- NoPhenotypes$MIM.Number
back2gene <- subset(mim2gene, MIM.Number %in% NoPhenotypeList)
pathwayanalysis <- as.character(back2gene$Ensembl.Gene.ID)
write(pathwayanalysis, file = "/Users/steliostheophanous/Desktop/genelistHomosapiens.txt", sep = " ")

#Genes NOT found in OMIM
reverse <- Reduce(intersect, list(finalresults$omim_id, phenotypes$MIM.Number))
notCommon <- setdiff(finalresults$omim_id, reverse)
notCommon2 <- subset(mim2gene, MIM.Number %in% notCommon)
notinOMIM <- notCommon2[!(is.na(notCommon2$Ensembl.Gene.ID) | notCommon2$Ensembl.Gene.ID==""), ]
write(as.character(notinOMIM$Ensembl.Gene.ID), file = "/Users/steliostheophanous/Desktop/notinOMIMhomosapiens.txt", sep = " ")

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
write(as.character(newList3), file = "/Users/steliostheophanous/Desktop/NEWnotinOMIMhs.txt", sep = " ")