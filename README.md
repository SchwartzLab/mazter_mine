# MAZTER-mine

MAZTER-mine is a computational pipeline to analize MAZTER-seq derived data, 
a methodology that profiles m6A quantitatively across transcriptomes in a 
single-base manner.

This repository holds two programs to run the MAZTER-seq computational pipeline,
a set of helper functions, and an additional folder which includes a tutorial
and test data to run it:

### helperFunctions.R

An R script with helper functions. For an easier handling this functions are 
loaded from the online repository when running bam2ReadEnds.R and mazter_seq.R.
If you want to run these while offline, please change the path in the source() 
function at the beginning of the main programs to your local copy of 
"helperFunctions.R".

### bam2ReadEnds.R

This is an R script that should be used as a command line tool in a UNIX system, 
providing the necessary arguments.

e.g.
Rscript bam2ReadEnds.R -i Sample1.bam -g geneAnnot.bed

Its main input is an **unsorted** alignment file in BAM format and a gene 
annotation in [BED-12](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. 
This processing step may be lengthy, depending on computational power and file size.

```sh
Rscript bam2ReadEnds.R --help

Usage: bam2ReadEnds.R [options]

Options:
        -i CHARACTER, --BAMfile=CHARACTER
                BAM alignment file

        -g CHARACTER, --geneAnnotation=CHARACTER
                Gene annotation file in BED-12 format

        -l NUMERIC, --maxInsLen=NUMERIC
                Maximum mapped insert length [default= 200]

        -m NUMERIC, --minInsG=NUMERIC
                Minimum number of inserts to report gene [default= 15]

        -n NUMERIC, --nCores=NUMERIC
                Number of cores to be used [default= 2]

        -h, --help
                Show this help message and exit

```

### mazter_mine.R

This is a program that should be used as a command line tool in a UNIX providing the necessary arguments.

mazter_mine's main input is the ".Rdata" file output from the bam2ReadEnds.R program, 
and it outputs a QC report and a cleavage efficiency table which can be used for 
downstream statistical analysis.

```sh
Rscript master_mine.R --help
Usage: master_mine.R [options]


Options:
        -i CHARACTER, --countDataFile=CHARACTER
                .Rdata file with count data (bam2readEnds.R output)

        -g CHARACTER, --geneAnnotation=CHARACTER
                Gene annotation file in BED-12 format

        -f CHARACTER, --faGenome=CHARACTER
                FASTA genome file to retrieve gene seqs

        -c CHARACTER, --clvMotif=CHARACTER
                Cleavage motif to measure at [default= ACA]

        -m NUMERIC, --minCov=NUMERIC
                Minimum coverage to quantify site [default= 15]

        -u NUMERIC, --upSThr=NUMERIC
                Up-stream distance to closest motif threshold

        -d NUMERIC, --doSThr=NUMERIC
                Down-stream distance to closest motif threshold

        -n NUMERIC, --nCores=NUMERIC
                Number of cores to be used [default= 2]

        -t CHARACTER, --tagName=CHARACTER
                Tag used for output files naming

        -h, --help
                Show this help message and exit

```

### tutorial

An example MAZTER-mine run with test data.

Dependencies:

* Bedtools (tested using samtools 1.3.1)
* SAMtools (tested using bedtools 2.26.0)
