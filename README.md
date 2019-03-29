# MAZTER-mine

This repository holds three files to run the MAZTER-seq computational pipeline.

### helperFunctions.R

An R script with helper functions

### bam2ReadEnds.R

This is a script that should be used as a command line tool in a UNIX providing the necessary arguments.

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

Dependencies:

* Bedtools (tested using samtools 1.3.1)
* SAMtools (tested using bedtools v2.26.0)

### mazter_mine.R

This is the main function which outputs cleavage efficiency tables and QC reports.

This is a script that should be used as a command line tool in a UNIX providing the necessary arguments.

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
