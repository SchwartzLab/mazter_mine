#!/usr/bin/env Rscript
# Title: MASTER-MINE: Calculate cleavage efficiencies based on cleavage motif ##
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############
library("optparse")

# Parsing Arguments ############################################################
option_list = list(
    make_option(c("-i", "--countDataFile"), type="character", default = NULL,
                help=".Rdata file with count data (bam2readEnds.R output)", metavar="character"),
    make_option(c("-g", "--geneAnnotation"), type="character", default = NULL,
                help="Gene annotation file in BED-12 format", metavar="character"),
    make_option(c("-f", "--faGenome"), type="character", default = NULL,
                help="FASTA genome file to retrieve gene seqs", metavar="character"),
    make_option(c("-c", "--clvMotif"), type="character", default = "ACA",
                help="Cleavage motif to measure at [default= %default]", metavar="character"),
    make_option(c("-m", "--minCov"), type="integer", default = 15,
                help="Minimum coverage to quantify site [default= %default]", metavar="numeric"),
    make_option(c("-u", "--upSThr"), type="integer", default = NA,
                help="Up-stream distance to closest motif threshold", metavar="numeric"),
    make_option(c("-d", "--doSThr"), type="integer", default = NA,
                help="Down-stream distance to closest motif threshold", metavar="numeric"),
    make_option(c("-n", "--nCores"), type="integer", default = 2,
                help="Number of cores to be used [default= %default]", metavar="numeric"),
    make_option(c("-t", "--tagName"), type="character", default = "DEFAULT",
                help="Tag used for output files naming", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

nCores <- opt$nCores
countFile <- opt$countDataFile
geneAnnotation <- opt$geneAnnotation
faGenome <- opt$faGenome
clvMotif <- opt$clvMotif
upSThr <- opt$upSThr
doSThr <- opt$doSThr
rtThr <- as.numeric(opt$minCov)
nick <- opt$tagName
if(nick == "DEFAULT"){
    nick <- gsub(pattern = "_Aligned.out.Rdata", replacement = "", x = countFile)
}

#Functions #####################################################################
source("https://drive.google.com/uc?export=download&id=1IpatZ2_AK7Bwor6BzSpYyi6UeCfiA5xx")

# Packages #####################################################################
#CRAN packages
CRAN_packs <- c("magrittr", "parallel", "stringr", "RColorBrewer", "ggplot2",
                "grid", "gtools")
dummy <- sapply(CRAN_packs, function(x) installLoad_CRAN(x))

# MAIN program #################################################################
# Gene Annotation
geneAnnot <- readBED(geneAnnotation)
colnames(geneAnnot) <- c("chr", "start", "end", "name", "score", "strand",
                         "startCodon", "stopCodon", "itemRGB", "blockCount",
                         "blockSizes", "blockStarts")
rownames(geneAnnot) <- geneAnnot$name

# Retrieve gene sequences
tmpFName <- tempfile()
system(paste0("bedtools getfasta -s -split -fi ", faGenome,
              " -bed ", geneAnnotation, " > ", tmpFName, ".fseq.fa"))
dummy <- read.delim(paste0(tmpFName, ".fseq.fa"), header = F, stringsAsFactors = F)
txOmeSEQS <- dummy[!as.logical(1:nrow(dummy)%%2),]; rm(dummy)
if(length(txOmeSEQS) != length(geneAnnot$name)){
    stop("Number of gene sequences retrieved not equal to number of genes in annotation file")
}
names(txOmeSEQS) <- geneAnnot$name
txOmeSEQS <- sapply(txOmeSEQS, function(x) str_to_upper(x))
dummy <- (file.remove(paste0(tmpFName, ".fseq.fa")))

# MotifOme
motifOme <- lapply(txOmeSEQS, function(word) { 
    str_locate_all(word, clvMotif)[[1]][,1]
})
names(motifOme) <- geneAnnot[,4]
motifOme <- lapply(motifOme, function(x) x[x>1]) # Filter start of seq
motifOme <- motifOme[sapply(motifOme, function(x) length(x) > 0)] # Filter genes with no motif

# All Sites motif distance 
allSites <- data.frame(coorNames = unlist(sapply(1:length(motifOme), function(x){
    paste(names(motifOme)[x], motifOme[[x]], sep = "_")})),
    upS_motifDist = lapply(names(motifOme), function(iGene){
        if(length(motifOme[[iGene]]) == 1){return(motifOme[[iGene]] - 1)}
        return(motifOme[[iGene]] - c(1, motifOme[[iGene]][1:length(motifOme[[iGene]])-1]))}) %>% unlist,
    doS_motifDist = lapply(names(motifOme), function(iGene){if(length(motifOme[[iGene]]) == 1){
        return(nchar(txOmeSEQS[iGene]) - motifOme[[iGene]])}
        return(c(motifOme[[iGene]][2:length(motifOme[[iGene]])],
                 nchar(txOmeSEQS[iGene])) - motifOme[[iGene]])}) %>% unlist, stringsAsFactors = F)

allSites$upS_motifDist <- as.integer(allSites$upS_motifDist)
rownames(allSites) <- allSites$coorNames

# Raw cleavage efficiencies 
print("Calculating raw cleavage efficiencies")
print(paste("Minimum coverage required to measure cleavage site =", rtThr))

RAW <- rawClvEff_motif(countDataFile = countFile, RTthr = rtThr, motifOme = motifOme)
if(nrow(RAW) == 0){stop("Premature stop due to no measurable sites available")}

# Optimal distance to motif site, for 3'- 5' measurement correlation 
optMOT <- optiMotifDist(rawCE = RAW, distLimit = 100)

# Aggregate 3'& 5' data 
if(is.na(doSThr)){doSThr <- min(optMOT$minD[c(1,3)])}
if(is.na(upSThr)){upSThr <- min(optMOT$minD[c(1,2)])}
print(paste("Up-stream threshold used =", upSThr))
print(paste("Down-stream threshold used =", doSThr))
AVG_CL_EFF <- avgCleavEff(rawClvEff = RAW, doS_thr = doSThr, upS_thr = upSThr)
rownames(AVG_CL_EFF) <- AVG_CL_EFF$coorNames
if(nrow(AVG_CL_EFF) == 0){stop("Premature stop due to absence of sites after aggregation")}

# Preparing output tables 
cl <- makeCluster(nCores, type = "FORK")
genTABS <- parSapply(cl, AVG_CL_EFF$coorNames, relToGenCoors) %>% t
sampSeqs <- parSapply(cl, AVG_CL_EFF$coorNames, function(y){
    extracSeq(y, window =  10, txOmeSEQS)
})
stopCluster(cl)

outTables <- cbind(genTABS, AVG_CL_EFF, sampSeqs)
colnames(outTables)[1:6] <- colnames(geneAnnot)[1:6]
colnames(outTables)[17] <- "seqs"

# Write output table
write.table(x = outTables, file = paste0(nick, "_clvEffTable.txt"), quote = F, 
            sep = "\t", col.names = F)
writeMyTable(X = outTables, filename = paste0(nick, "_clvEffTable.txt"))

# Quality Control ##############################################################
# Measurements too close to cleavage motif
tmp <- outTables
tmp$farDo <- ifelse(tmp$doS_motifDist < doSThr, "do", "no")
tmp$farUp <- ifelse(tmp$upS_motifDist < upSThr, "up", "no")
tmp$tooCloseToMotif <- "None"
tmp$tooCloseToMotif[tmp$farDo == "do"] <- "DownStream"
tmp$tooCloseToMotif[tmp$farUp == "up"] <- "UpStream"
tmp$tooCloseToMotif[tmp$farDo == "do" & tmp$farUp == "up"] <- "Both"
tmp$tooCloseToMotif <- factor(tmp$tooCloseToMotif, levels = c("UpStream", "DownStream", "Both", "None"))
corUnF <- cor(tmp[,c("clvEff_5", "clvEff_3")] , use = "p")[2,1] %>% round(2)
corNo <- cor(subset(tmp, tooCloseToMotif == "None")[,c("clvEff_5", "clvEff_3")] , use = "p")[2,1] %>% round(2)

gg1 <- ggplot(tmp, aes(x = clvEff_5, y = clvEff_3, colour = tooCloseToMotif)) + 
    geom_point(size = 1) + 
    scale_color_manual(values=c(
        adjustcolor("blue", alpha.f = 0.2),
        adjustcolor("red", alpha.f = 0.2),
        adjustcolor("purple", alpha.f = 0.2),
        adjustcolor("black", alpha.f = 0.4)
    )) + theme_classic() + xlab("Raw 5' Cleavage Efficiency") + 
    ylab("Raw 3' Cleavage Efficiency") + 
    ggtitle("Correlation 3' and 5' signal", 
            paste("Unfiltered Correlation =", corUnF, "\nFiltered Correlation =", corNo))

tmp2 <- subset(tmp, doS_motifDist >= doSThr)
gg2 <- ggplot(tmp2, aes(x = clvEff_5, y = avgClvEff)) + geom_bin2d(bins =100) +
    scale_fill_gradientn(trans = "log", colours = rev(brewer.pal(9, "Spectral"))) +
    theme_classic() + ylab("Average Cleavage Efficiency") + xlab("Raw 5' Cleavage Efficiency")

tmp2 <- subset(tmp, upS_motifDist >= upSThr)
gg3 <- ggplot(tmp2, aes(x = clvEff_3, y = avgClvEff)) + geom_bin2d(bins = 100) +
    scale_fill_gradientn(trans = "log", colours = rev(brewer.pal(9, "Spectral"))) +
    theme_classic() + ylab("Average Cleavage Efficiency") + xlab("Raw 3' Cleavage Efficiency")

# Cleavage motif frequency
cl <- makeCluster(nCores, type = "FORK")
nF <- clvNucFreq(cl, countFile = countFile, txOmeSEQS, nchar(clvMotif))
stopCluster(cl)

ggPie_3 <- ggPieNucFreq(nF$nucsFreq_3, 3) + ggtitle("Top 5 N-Nucleotides library composition", "3' read-ends")
ggPie_5 <- ggPieNucFreq(nF$nucsFreq_5, 3) + ggtitle("Top 5 N-Nucleotides library composition", "5' read-ends")

# Write QC-plots
pdf(paste0(nick, "_MASTER-seq_QC.pdf"), width = 12, height = 4)
plot_OptimMotifDist(optMOT, nick)
multiplot(gg1, gg2, gg3, cols = 3)
multiplot(ggPie_5, ggPie_3, cols = 2)
dev.off()

print("MASTER-MINE analysis complete!")
(file.remove("Rplots.pdf"))
