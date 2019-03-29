#!/usr/bin/env Rscript
# Title: Start-End Paired BED genomic to transcriptomic exon annotation already generated ########
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############
library("optparse")

# Parsing Arguments ############################################################
option_list = list(
    make_option(c("-i", "--BAMfile"), type="character", default = NULL,
                help="BAM alignment file", metavar="character"),
    make_option(c("-g", "--geneAnnotation"), type="character", default = NULL,
                help="Gene annotation file in BED-12 format", metavar="character"),
    make_option(c("-l", "--maxInsLen"), type="integer", default = 200,
                help="Maximum mapped insert length [default= %default]", metavar="numeric"),
    make_option(c("-m", "--minInsG"), type="integer", default = 15,
                help="Minimum number of inserts to report gene [default= %default]", metavar="numeric"),
    make_option(c("-n", "--nCores"), type="integer", default = 2,
                help="Number of cores to be used [default= %default]", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

nCores <- opt$nCores
BAMfile <- opt$BAMfile
geneAnnotation <- opt$geneAnnotation
maxInsertLength <- opt$maxInsLen
minIns <- opt$minInsG

#Functions #####################################################################
source("https://drive.google.com/uc?export=download&id=1IpatZ2_AK7Bwor6BzSpYyi6UeCfiA5xx")

# Packages #####################################################################
#CRAN packages
CRAN_packs <- c("magrittr", "parallel", "optparse", "plyr")
sapply(CRAN_packs, function(x) installLoad_CRAN(x))

# MAIN program #################################################################
# Create bed paired end alignment ####
BEDPfile <- tempfile()
comm <- paste0("samtools view -bf 0x2 ", BAMfile," | bedtools bamtobed -i stdin -bedpe -mate1 | cut -f1-3,5-6,9 > ", BEDPfile)
system(command = comm, wait = T)

# Annotation file ####
geneAnnot <- readBED(geneAnnotation)
colnames(geneAnnot) <- c("chr", "start", "end", "name", "score", "strand", 
                         "startCodon", "stopCodon", "itemRGB", "blockCount", 
                         "blockSizes", "blockStarts")
rownames(geneAnnot) <- geneAnnot$name

# Process sample ####
cl <- makeCluster(nCores, type = "FORK")
procSample <- pairEndReadsToGenes(cl, bedPEGzFile = BEDPfile, gAnnot = geneAnnot)
stopCluster(cl)

# Get count data ####
cl <- makeCluster(nCores, type = "FORK")
countData <- getCountData(cluster = cl, 
             dat = procSample, 
             GENES = names(procSample),
             geneAnnot = geneAnnot, 
             minInserts = minIns,
             maxInsLen = maxInsertLength) 
stopCluster(cl)
names(countData) <- names(procSample)
tmpSel <- sapply(countData, function(x) length(x$cov)) > 1
countData <- countData[tmpSel]

# Save Rdata object ############################################################
newFileName <- gsub(pattern = ".bam", replacement =  ".Rdata", x = BAMfile)
save(countData, file = newFileName)
file.remove(BEDPfile) #remove bedPfile
