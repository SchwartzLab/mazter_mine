# Helper Functions for MAZTER-seq ############################################################
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############

# Write table with TAB separated files standards
writeMyTable <- function(X, filename){
  write.table(X, file = filename, sep= "\t", col.names = NA, quote = F)
}

# Extracting the last n characters from a string in R
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# List files in work dir that match pattern pat
listFilePatt <- function(pattern, path = "."){
  files <- list.files(path)[grep(pattern = pattern, x = list.files(path))]
  return(files)
}

# Object size
oSize <- function(x){
  print(object.size(x), units = "auto") 
}

# Turns a character matrix into a dataframe more "smartly"
asDataFrame <- function(x){
  tmpCol <- colnames(x)
  out <- data.frame(lapply(split(x, col(x)), type.convert, as.is = TRUE),
             stringsAsFactors = FALSE)
  colnames(out) <- tmpCol
  return(out)
}

# Crate BED file with range of positions for sequence retrieval with BEDTOOLS
bedRange <- function(bedDF, d){
    data.frame(bedDF[,1], as.integer(bedDF[,3]) - d -1, 
               as.integer(bedDF[,3]) + d, bedDF[,5],
               bedDF[,4], bedDF[,6], stringsAsFactors = F)
}

# Read a BED-6 file and transform base-0 coordinates to base-1
readBED <- function(fileName, idAsRowNames = T){
    tmp <- read.delim(fileName, header = F, stringsAsFactors = F)
    tmp[,2] <- tmp[,2] + 1
    colnames(tmp)[1:6] <- c("chr", "start", "end", "id", "score", "strand")
    tmp$strand <- as.character(tmp$strand)
    tmp$chr <- as.character(tmp$chr)
    if(idAsRowNames){rownames(tmp) <- tmp$id}
    return(tmp)
}

# Function write bed file (Input = base1, output = base-0)
writeBED <- function(bed, filename){
    bed[,2] <- bed[,2] -1
    write.table(bed, filename, sep= "\t", col.names = F, row.names = F, quote = F)
}


# Crate BED file with range of positions for sequence retrieval with BEDTOOLS
bedRange <- function(bedDF, d){
    cbind(as.character(bedDF[,1]), as.integer(bedDF[,3]) - d -1, 
          as.integer(bedDF[,3]) + d, 
          as.character(bedDF[,5]),
          as.character(bedDF[,4]),
          as.character(bedDF[,6]))
}

# Some logical functions
great <- function(x, Than, orEqual = F, value = F){
    if(value){
        if(orEqual){x[x >= Than]}else if(!orEqual){x[x > Than]}
    }else if(!value){
        if(orEqual){x >= Than}else if(!orEqual){x > Than}
    }
}

less <- function(x, Than, orEqual = F, value = F){
    if(value){
        if(orEqual){x[x <= Than]}else if(!orEqual){x[x < Than]}
    }else if(!value){
        if(orEqual){x <= Than}else if(!orEqual){x < Than}
    }
}


equal <- function(x,y, value = F){
    if(value){x[x == y]}else if(!value){x == y}
}

# Function to go up along a vector until finding lower values, with a threshold of fails going up
stepWiseUp <- function(x, failThr, start = 2){
    best <- start
    fails <- 0
    for(i in start:length(x)){
        if(fails > failThr){break()}
        if(x[i] > x[i-1]){
            best <- i
            fails <- 0
        }else{fails <- fails + 1}
    }
    return(best)
}


# Function to change the range of numeric values from 0 to 1 - Feature Scaling
range0_1 <- function(x){
    z <- (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
    return(z)
}

# Function for plotting a 2 axis line plot with baseR
# Based on https://stackoverflow.com/questions/6142944/how-can-i-plot-with-2-different-y-axes#6143251
secAxisPlot <- function(xAxis, yAxis, yAxis2, xlab = "", ylab = "", y2lab = "", main = "", secCol= 2){
    ## add extra space to right margin of plot within frame
    par(mar=c(5, 5, 4, 6) + 0.1)
    ## Plot first set of data and draw its axis
    plot(xAxis, yAxis, pch = 16, axes = FALSE, ylim = range(yAxis), xlab= "", ylab = "", 
         type = "b", col = "black", main = main)
    axis(2, ylim = range(yAxis), col = "black", las=1)  ## las=1 makes horizontal labels
    mtext(ylab, side = 2, line = 4)
    box()
    ## Allow a second plot on the same graph
    par(new = TRUE)
    ## Plot the second plot and put axis scale on right
    plot(xAxis, yAxis2, pch=15, xlab="", ylab = "", ylim = range(yAxis2), 
         axes= FALSE, type="b", col= secCol)
    ## a little farther out (line=4) to make room for labels
    mtext(y2lab, side = 4, col = secCol, line=4) 
    axis(4, ylim = range(yAxis2), col = secCol, col.axis= secCol, las=1)
    ## Draw the xAxis axis
    axis(1, pretty(range(xAxis), 10))
    mtext(xlab, side = 1, col = "black", line = 2.5)  
}

#Merge two data frames by rownames and renames rows
mergebyRowNames <- function(x, y, keepAll = F){
    tmp <- merge(x, y, by = "row.names", all = keepAll)
    rownames(tmp) <- tmp[,1]
    tmp <- tmp[,-1]
    return(tmp)
}

# Merge by rownames using dplyr
fullJoinByRowNames <- function(x, y, keepAll = F){
    x$rowNames <- rownames(x)
    y$rowNames <- rownames(y)
    tmpO <- dplyr::full_join(x, y, by = "rowNames")
    rownames(tmpO) <- tmpO$rowNames
    return(tmpO[,!colnames(tmpO) == "rowNames"])
}
leftJoinByRowNames <- function(x, y, keepAll = F){
    x$rowNames <- rownames(x)
    y$rowNames <- rownames(y)
    tmpO <- dplyr::left_join(x, y, by = "rowNames")
    rownames(tmpO) <- tmpO$rowNames
    return(tmpO[,!colnames(tmpO) == "rowNames"])
}

# Read a BEDPE file and transform base-0 coordinates to base-1
readBEDPE <- function(fileName, idAsRowNames = T){
    tmp <- read.delim(fileName, header = F, stringsAsFactors = F)
    colnames(tmp) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2")
    tmp$start1 <- tmp$start1 + 1
    tmp$start2 <- tmp$start2 + 1
    if(idAsRowNames){rownames(tmp) <- tmp$name}
    return(tmp)
}

# Read BedPE file for getting read counts (strand specific)
readBEDPE3 <- function(fileName){
    tmp <- read.delim(fileName, header = F, stringsAsFactors = F)
    colnames(tmp) <- c("chr", "start_r1", "end_r1", "start_r2", "end_r2", "strand")
    tmp$start_r1 <- tmp$start_r1 + 1
    tmp$start_r2 <- tmp$start_r2 + 1
    return(tmp)
}

# Generate exon coordinates block
exonBlockGen <- function(iGene, geneAnnot){
    iStart <- geneAnnot[iGene,]$start
    iEnd <- geneAnnot[iGene,]$end
    iStrand <- geneAnnot[iGene,]$strand
    tmpA <- geneAnnot[iGene,]$blockStarts %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
    tmpB <- geneAnnot[iGene,]$blockSizes %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
    iBlocks <- sapply(1:length(tmpA), function(i){
        c(tmpA[i],(tmpA[i] + tmpB[i] -1))
    }) %>% t
    iBlocks <- iStart + iBlocks
    if(iEnd %in% as.vector(iBlocks)){
        if(iStrand == "+"){
            sapply(1:nrow(iBlocks), function(k){
                iBlocks[k,1]:iBlocks[k,2]
            }) %>% unlist
        }else if(iStrand == "-"){
            sapply(1:nrow(iBlocks), function(k){
                iBlocks[k,1]:iBlocks[k,2]
            }) %>% unlist %>% rev
        }
    }else{stop(paste("Malformed exon structure at gene", iGene))}
}

# Function paired end BEDfile to transcriptome
pairEndReadsToGenes <- function(parCluster, bedPEGzFile, gAnnot){
    inputFile <- readBEDPE3(gzfile(bedPEGzFile))
    commChrom <- intersect(unique(gAnnot$chr), unique(inputFile$chr))
    genPReads <- parLapply(parCluster, commChrom, function(iChr){
        tmpIF <- subset(inputFile, chr == iChr)
        tmpGA <- subset(gAnnot, chr == iChr)
        tmpO <- lapply(tmpGA$name, function(iGene){
            exBlock <- exonBlockGen(iGene, geneAnnot)
            if(geneAnnot[iGene, "strand"] == "+"){
                tmpDF <- subset(tmpIF, strand == "+" & start_r1 %in% exBlock & end_r2 %in% exBlock)
                tmpDF$start <- match(tmpDF$start_r1, exBlock)
                tmpDF$end <- match(tmpDF$end_r2, exBlock)
                tmpDF <- tmpDF[,c("chr", "start", "end")]
            }else if(geneAnnot[iGene, "strand"] == "-"){
                tmpDF <- subset(tmpIF, strand == "-" & start_r2 %in% exBlock & end_r1 %in% exBlock)
                tmpDF$start <- match(tmpDF$end_r1, exBlock)
                tmpDF$end <- match(tmpDF$start_r2, exBlock)
                tmpDF <- tmpDF[,c("chr", "start", "end")]
            }
            tmpDF <- ddply(tmpDF,.(chr, start, end), nrow)
            if(ncol(tmpDF) == 4){
                colnames(tmpDF)[4] <- "times"
                tmpDF$len <- tmpDF$end - tmpDF$start
            }
            rownames(tmpDF) <- NULL
            return(tmpDF)
        })
        names(tmpO) <- tmpGA$name
        return(tmpO)
    })
    # Flatten list
    tmp <- do.call(c, genPReads)
    genPReads <- tmp
    return(genPReads)
}

# RowMeans by column groups
rowMeansColG <- function(DF, colGroups, na.rm = T){
    cG <- unique(colGroups)
    out <- sapply(cG, function(x){
        rowMeans(DF[,colGroups == x], na.rm = na.rm)
    }) %>% as.data.frame()
    colnames(out) <- cG
    rownames(out) <- rownames(DF)
    return(out)
}

# Extract countData from DATA object
getCountData <- function(cluster, dat, GENES, geneAnnot, minInserts, maxInsLen){
    rawData <- parLapply(cluster, GENES, function(iGene) {
        iStrand <- geneAnnot[iGene, "strand"]
        tmpD <- dat[[iGene]]
        tmpD <- subset(tmpD, tmpD$len <= maxInsLen)
        if(nrow(tmpD) == 0){
            tmpRE <- list(rE5 = NA, rE3 = NA, cov = NA)
        }else if(sum(tmpD$times) < minInserts){
            tmpRE <- list(rE5 = NA, rE3 = NA, cov = NA)
        }else{
            tmpRE_5 <- tmpRE_3 <- tmpRE_c <- rep(0L, length(exonBlockGen(iGene, geneAnnot)))
            for(k in 1:nrow(tmpD)){
                tmpRE_5[tmpD[k,]$start] <- tmpRE_5[tmpD[k,]$start] + tmpD[k,]$times
                tmpRE_3[tmpD[k,]$end] <- tmpRE_3[tmpD[k,]$end] + tmpD[k,]$times
                tmpRE_c[tmpD[k,]$start:tmpD[k,]$end] <- tmpRE_c[tmpD[k,]$start:tmpD[k,]$end] + tmpD[k,]$times
            }
            tmpRE <- list(rE5 = tmpRE_5, rE3 = tmpRE_3, cov = tmpRE_c)
        }
        return(tmpRE)
    })
    return(rawData)
}


# Load/install CRAN packages
installLoad_CRAN <- function(package){
    if (!require(package, character.only = T)) {
        install.packages(package, dependencies = TRUE)
        library(package, character.only = T, quietly = T)
    }
}

# Load/install Bioconductor packages
installLoad_BioC <- function(package){
    if (!require(package, character.only = T)) {
        if(version$minor %>% as.numeric() >= 3.5){
            if (!requireNamespace("BiocManager", quietly = TRUE)){
                installLoad_CRAN("BiocManager")
            }
            BiocManager::install(package)
            library(package, character.only = T, quietly = T)
        }
        if(version$minor %>% as.numeric() < 3.5){
            source("https://bioconductor.org/biocLite.R")
            biocLite(package, ask = F)
            library(package, character.only = T, quietly = T)
        }
    }
}

# Number of inserts in sample
insertsNum <- function(pairEndReads){
    sapply(pairEndReads, function(x){
        sum(x$times)
    })
}

# Get all the relative coordinates for all genes using a BED-like matrix row
getRelPos <- function(BEDrow){
    tmp <- BEDrow
    iChr <- tmp$chr
    iStart <- tmp$start
    iStrand <- tmp$strand
    tmpA <- subset(geneAnnot, chr == iChr & strand == iStrand & start <= iStart & end >= iStart)
    tmpCoor <- sapply(tmpA$name, function(x){
        tmp <- which(exonBlocks[[x]] == iStart)
        if(!as.logical(length(tmp))){
            tmp <- NA
        }
        return(tmp)
    }) %>% unlist
    paste0(tmpA$name,"_", tmpCoor)[!is.na(tmpCoor)]
}

# From relative to genomic coordinate
relToGenCoors <- function(relCoor, genAnnot = geneAnnot){
    splt <- unlist(strsplit(relCoor, split = "_"))
    gene <- splt[1] %>% as.character()
    relPos <- splt[2] %>% as.numeric()
    chr <- genAnnot[gene,]$chr
    exBlock <- exonBlockGen(gene, genAnnot)
    strand <- genAnnot[gene,]$strand
    genCoor <- exBlock[relPos]
    return(c(chr, genCoor, genCoor, relCoor, 0, strand))
}

# Extract window of sequences from gene sequence vector
extracSeq <- function(relCoor, window = 10, fullSeqs = txOmeSEQS){
    tmp <- unlist(strsplit(relCoor, split = "_"))
    iGene <- tmp[1]
    iCoor <- as.numeric(tmp[2])
    substr(fullSeqs[iGene], start = iCoor - window, stop =  iCoor + window)
}

# raw cleavage efficiencies in motif (motifOme object required)
rawClvEff_motif <- function(countDataFile, RTthr = 15, motifOme = motifOme){
    load(countDataFile)
    tmpSel <- sapply(countData, function(x) length(x$cov)) > 1
    countData <- countData[tmpSel]
    genesI <- intersect(names(countData), names(motifOme))
    rawClvEff <- lapply(genesI, function(geneI){
        tmp <- data.frame(coorNames = paste0(geneI,"_", motifOme[[geneI]]),
                          rE5 = (countData[[geneI]]$rE5)[motifOme[[geneI]]],
                          rT5 = (countData[[geneI]]$cov)[motifOme[[geneI]]],
                          rE3 = (countData[[geneI]]$rE3)[motifOme[[geneI]] -1],
                          rT3 = (countData[[geneI]]$cov)[motifOme[[geneI]] -1], stringsAsFactors = F)
        tmp$clvEff_5 <- tmp$rE5 / tmp$rT5
        tmp$clvEff_3 <- tmp$rE3 / tmp$rT3
        tmp$clvEff_5[tmp$rT5 < RTthr] <- NA
        tmp$clvEff_3[tmp$rT3 < RTthr] <- NA
        out <- tmp[(is.na(tmp[,c("clvEff_5", "clvEff_3")]) %>% rowSums()) < 2,]
        return(out)
    })
    rawClvEff <- do.call(rbind, rawClvEff)
    rawClvEff <- merge(allSites, rawClvEff)
    return(rawClvEff)
}

# Optimal distance analysis
optiMotifDist <- function(rawCE, distLimit){
    TMP_both <- sapply(1:distLimit, function(thr){
        TMPdata <- cbind(rawCE$clvEff_5[rawCE$doS_motifDist >= thr & rawCE$upS_motifDist >= thr],
                         rawCE$clvEff_3[rawCE$doS_motifDist >= thr & rawCE$upS_motifDist >= thr]) %>% na.omit
        c(cor(TMPdata[,1], TMPdata[,2]), nrow(TMPdata))
    })
    
    minBoD <- stepWiseUp(TMP_both[1,], failThr = 1, start = 5)
    TMP_up <- sapply(1:distLimit, function(thr){
        TMP_up <- cbind(rawCE$clvEff_5[rawCE$doS_motifDist >= minBoD & rawCE$upS_motifDist >= thr],
                        rawCE$clvEff_3[rawCE$doS_motifDist >= minBoD & rawCE$upS_motifDist >= thr]) %>% na.omit
        c(cor(TMP_up[,1], TMP_up[,2]), nrow(TMP_up))
    })
    TMP_do <- sapply(1:distLimit, function(thr){
        TMP_do <- cbind(rawCE$clvEff_5[rawCE$doS_motifDist >= thr & rawCE$upS_motifDist >= minBoD],
                        rawCE$clvEff_3[rawCE$doS_motifDist >= thr & rawCE$upS_motifDist >= minBoD]) %>% na.omit
        c(cor(TMP_do[,1], TMP_do[,2]), nrow(TMP_do))
    })
    minUpD <- stepWiseUp(TMP_up[1,], failThr = 1, start = 10)
    minDoD <- stepWiseUp(TMP_do[1,], failThr = 1, start = 10)
    tmp <- c(minBoD, minUpD, minDoD); names(tmp) <- c("minBoD", "minUpD", "minDoD")
    return(list(boD = TMP_both, upD = TMP_up, doD = TMP_do, minD = tmp))
}

# Plotting optimal distance analyisis
plot_OptimMotifDist <- function(opMoDi, smpName){
    palette(brewer.pal(9, "Set1")) #RColorBrewer
    par(mfrow=c(1,3))
    distLimit <- ncol(opMoDi$boD)
    secAxisPlot(1:distLimit, opMoDi$boD[2,], opMoDi$boD[1,],
                xlab = "Distance to ACA",
                ylab = "n ACA sites eval",
                y2lab = "Pearson Correlation",
                main = "P. Corr ~ ACA dist.-BothS")
    abline(v = opMoDi$minD["minBoD"], col = "red")
    abline(h = opMoDi$boD[1, opMoDi$minD["minBoD"]], col = "red")
    
    secAxisPlot(1:distLimit, opMoDi$upD[2,], opMoDi$upD[1,],
                xlab = "Distance to ACA",
                ylab = "n ACA sites eval",
                y2lab = "Pearson Correlation",
                main = "P. Corr ~ ACA dist.-UpS")
    title(sub = smpName)
    abline(v = opMoDi$minD["minUpD"], col = "red")
    abline(h = opMoDi$upD[1, opMoDi$minD["minUpD"]], col = "red")
    
    secAxisPlot(1:distLimit, opMoDi$doD[2,], opMoDi$doD[1,],
                xlab = "Distance to ACA",
                ylab = "n ACA sites eval",
                y2lab = "Pearson Correlation",
                main = "P. Corr ~ ACA dist.-DownS")
    abline(v = opMoDi$minD["minDoD"], col = "red")
    abline(h = opMoDi$doD[1, opMoDi$minD["minDoD"]], col = "red")
    par(mfrow=c(1,1))
}

# Calculating average cleavage efficiency using rawClvEff_motif output and 
# downStream and upStream thresholds
avgCleavEff <- function(rawClvEff, doS_thr, upS_thr){
    tmp <- rawClvEff
    tmp$clvEff_5[tmp$doS_motifDist < doS_thr] <- NA
    tmp$clvEff_3[tmp$upS_motifDist < upS_thr] <- NA
    tmp$avgClvEff <- rowMeans(tmp[,c("clvEff_5", "clvEff_3")], na.rm = T)
    tmp$avgClvEff[is.na(tmp$avgClvEff)] <- NA
    out <- rawClvEff
    out$avgClvEff <- tmp$avgClvEff
    return(out)
}

# Cleavage motif frequency (based on first 'N' nucleotides)
clvNucFreq <- function(cluster, countFile, geneSeqs, nNucs){
    load(countFile)
    triNucs <- permutations(n = 4,r = nNucs, v = c("A", "T", "C", "G"),
                            repeats.allowed = T) %>% apply(MARGIN = 1, function(x){paste(x, collapse = "")})
    iGenes <- intersect(names(countData)[sapply(countData, function(x) sum(x$rE3)>0)], names(geneSeqs))
    nucsFreq_5 <- parSapply(cluster, iGenes, function(iG){
        test <- sapply(which(countData[[iG]]$rE5 > 0), function(i){
            c(substring(text = geneSeqs[iG], first = i, last = i+nNucs-1), countData[[iG]]$rE5[i])
        }) %>% t %>% asDataFrame()
        if(is.null(names(test))){
            test <- data.frame(a= NA, b = NA)
        }
        names(test) <- c("nucs", "times")
        sapply(triNucs, function(x){
            sum(subset(test, nucs == x)$times)
        })
    }) %>% rowSums()
    nucsFreq_3 <- parSapply(cluster, iGenes, function(iG){
        test <- sapply(which(countData[[iG]]$rE3 > 0), function(i){
            c(substring(text = geneSeqs[iG], first = i+1, last = i+nNucs), countData[[iG]]$rE3[i])
        }) %>% t %>% asDataFrame()
        if(is.null(names(test))){
            test <- data.frame(a= NA, b = NA)
        }
        names(test) <- c("nucs", "times")
        sapply(triNucs, function(x){
            sum(subset(test, nucs == x)$times)
        })
    }) %>% rowSums()
    list(nucsFreq_5= nucsFreq_5, nucsFreq_3 = nucsFreq_3)
}

# Which positions are top N
whichTopN <- function(x, N){
    names(x) <- 1:length(x)
    as.numeric(names(tail(sort(x), N)))
}
# Which positions are bottom N
whichBottomN <- function(x, N){
    names(x) <- 1:length(x)
    as.numeric(names(head(sort(x), N)))
}
# Omit infinite values
inf.omit <- function(x){
    x[!is.infinite(x)]
}

# Keep rows and column names when using normalize.quantiles()
normalize.quantiles.keepNames <- function(x){
    x[is.na(x)] <- NA
    require(preprocessCore)
    rN <- rownames(x)
    cN <- colnames(x)
    tmpO <- data.frame(normalize.quantiles(as.matrix(x)))
    rownames(tmpO) <- rN; colnames(tmpO) <- cN
    return(tmpO)
}

# Clear temporary variables
clearTmp <- function(){rm(list = grep("tmp", ls(), value = T)); gc(verbose = T)}

# GGPLOT functions #############################################################
# Multiplot for ggplot2 by winston@stdout.org from Cookbook for R
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

# Boxplot shortcut function (ggplot2)
ggBoxplot <- function(matrix, title = "Title", xlab = "x", ylab = "y", outLCol = NA){
    ggplot(data=melt(as.data.frame(matrix)), aes(variable, value)) + 
        geom_boxplot(outlier.colour= outLCol, outlier.size = 1) + xlab(xlab) + ylab(ylab) +
        ggtitle(title) + theme_classic() +  stat_n_text(size = 3, angle = 90) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Pie Chart Nucleotide frequency
ggPieNucFreq <- function(nucFreq, labSize = 5){
    palette(brewer.pal(9, "Set1")) #RColorBrewer
    tmpDF <- data.frame(nucs = names(nucFreq), Percent = nucFreq, stringsAsFactors = F)
    tmpDF <- tmpDF[order(tmpDF[,2], decreasing = T),]
    tmpDF <- data.frame(rbind(tmpDF[1:5,], c("All Others", sum(tmpDF[-1:-5,2]))))
    tmpDF[,2] <- as.numeric(tmpDF[,2]) / sum(as.numeric(tmpDF[,2]))
    tmpDF[,1] <- factor(tmpDF[,1], levels = tmpDF[,1])
    ggPie <- ggplot(tmpDF, aes(x="", y=Percent, fill=nucs)) +
        geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y",start = 0,direction = 1) + 
        geom_text(aes(label = round(Percent,2)), size= labSize, position = position_stack(vjust = 0.5)) +
        theme(axis.text.x =element_blank()) + theme_classic()
    return(ggPie)
}

# GGplot alternative to pairs function (additionally it fits linear models to all pair-wise comparisons)
ggPairs <- function(DF, alpha = 1){
    iCol <- colnames(DF)
    matD <- combinations(n = length(iCol), r = 2, v = 1:length(iCol))  
    ggSC <- lapply(1:nrow(matD), function(x){
        tmpL <- lm(DF[,matD[x,2]] ~ DF[,matD[x,1]])
        if(tmpL$coefficients[1]>=0){
            linModEq = paste("y = x *",tmpL$coefficients[2] %>% signif(2), "+", tmpL$coefficients[1]  %>% signif(2))
        }else if(tmpL$coefficients[1]<0){linModEq = paste("y = x *", signif(tmpL$coefficients[2],2), "-", 
                                                           tmpL$coefficients[1]  %>% signif(2) %>% abs)}
        tmpC <- cor(DF[,matD[x,1]], DF[,matD[x,2]], use = "p") %>% round(4)
        tmpP <- cor.test(DF[,matD[x,1]], DF[,matD[x,2]], use = "p")$p.value %>% signif(4)
        tmpC2 <- cor(DF[,matD[x,1]], DF[,matD[x,2]], use = "p", method = "spearman") %>% round(4)
        tmpP2 <- cor.test(DF[,matD[x,1]], DF[,matD[x,2]], use = "p", method = "spearman")$p.value %>% signif(4)
        ggplot(DF, aes(x= DF[,matD[x,1]], y= DF[,matD[x,2]])) + 
            geom_point(alpha = alpha, shape = 16) + 
            geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1)) +
            geom_abline(intercept = 0, slope = 1, colour = "gray") + 
            theme_classic() + xlab(iCol[matD[x,1]]) + ylab(iCol[matD[x,2]]) + 
            ggtitle(paste("R =", tmpC, "p = ", tmpP, "rho =", tmpC2, "p =", tmpP2), 
                    subtitle = linModEq) +
            coord_cartesian(ylim = range(DF, na.rm = T), xlim = range(DF, na.rm = T))
    })
    ggLabs <- lapply(iCol, function(x){
        df <- data.frame(x = 1, y = 1, text = x)
        ggO <- ggplot(df, aes(x, y)) +
            geom_text(aes(label = text), size = 5) + theme_classic() +
            theme(panel.border = element_rect(colour = 1, fill = NA), axis.line = element_line())+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) + 
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
        return(ggO)
    })  
    ggCf <- lapply(1:nrow(matD), function(x){
        return()
    })  
    lOut <- matrix(NA, ncol = ncol(DF), nrow = ncol(DF))
    for(i in 1:nrow(matD)){lOut[matD[i,2], matD[i,1]] <- i}
    for(i in 1:length(iCol)){lOut[i, i] <- length(ggSC) + i}
    for(i in 1:nrow(matD)){lOut[matD[i,1], matD[i,2]] <- length(ggSC) + length(iCol) + i}
    multiplot(plotlist = c(ggSC, ggLabs), layout = lOut)
}

# Simple Barplot function
ggBarplot <- function(x, ci = NA, title = NULL, subt = NULL, xLab = "Names", yLab = "Values"){
    if(is.null(names(x))){names(x) <- 1:length(x)}
    df <- data.frame(names = names(x), value=x, CI = ci)
    outGG <- ggplot(data=df, aes(x=names, y=value)) + 
        geom_bar(stat="identity") + theme_classic() + 
        geom_errorbar(aes(ymin=value-CI, ymax=value+CI), width=.2, position=position_dodge(.9)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(title, subt) +
        ylab(yLab) + xlab(xLab)
    return(outGG)
}

# Scatterplot with linear model fitted line and calculates correlation
ggScattLinePlot <- function(x, y, title = "", xLab = "", yLab = "", alpha = 1){
    tmpC <- cor(x, y, use = "p") %>% round(4)
    tmpP <- cor.test(x, y, use = "p")$p.value %>% signif(3)
    tmpC2 <- cor(x, y, use = "p", method = "spearman") %>% round(4)
    tmpP2 <- cor.test(x, y, use = "p", method = "spearman")$p.value %>% signif(3)
    tmpDF <- data.frame(var1 = x, var2 = y)
    ggSCLINE <-  ggplot(tmpDF, aes(x = var1, y = var2)) + geom_point(alpha = alpha) + 
        geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1)) +
        ggtitle(title, paste("R =", tmpC, "p = ", tmpP, "\nrho =", tmpC2, "p =", tmpP2)) +
        ylab(yLab) + xlab(xLab) + theme_classic()
    return(ggSCLINE)
}


# Load/install dependencies
CRAN_packs <- c("magrittr", "parallel", "optparse", "plyr", "ggplot2", "grid", "gtools", "reshape2", "EnvStats")
sapply(CRAN_packs, function(x) installLoad_CRAN(x))
