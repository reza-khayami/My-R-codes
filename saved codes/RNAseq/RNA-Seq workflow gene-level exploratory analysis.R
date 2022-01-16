#doi: 10.12688/f1000research.7035.1     read if you have problems!

# Locating alignment files ---------------------------------------------------------------

library(airway)

# The R function system.file can be used to find out where on your computer the files from a package have been installed.
dir <- system.file("extdata", package="airway", mustWork=TRUE)


# In this directory, we find the eight BAM files (and some other files):
list.files(dir)


# Typically, we have a table with detailed information for each of our samples that links samples to the associated FASTQ and BAM files. For your own project, you might create such a comma-separated value (CSV) file using a text editor or spreadsheet software such as Excel.

csvfile <- file.path(dir,"sample_table.csv") #load datail table
(sampleTable <- read.csv(csvfile,row.names=1))
# the parentheses () around the entire code of the last line above is an R trick to print the output of the function in addition to saving it to sampleTable

# The following tools can be used generate count matrices: summarizeOverlaps, featureCounts, or htseq-count

# We now proceed using the summarizeOverlaps method of counting.

# construct the full paths to the files

filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)

library(Rsamtools)

bamfiles <- BamFileList(filenames, yieldSize=2000000) 
# only process 2 million reads at a time

#We need to check if the chromosome names are matched between bam and gtf (here called "seqnames")
#in the alignment files
seqinfo(bamfiles[1])


# Defining gene models ----------------------------------------------------
# will read the gene model from an Ensembl GTF file
# TxDb object : a database that can be used to generate
# a variety of range-based objects, such as exons, transcripts, and genes

# There are other options for constructing a TxDb. For the known genes track from the UCSC Genome Browser

# one can use the pre-built Transcript DataBase: TxDb.Hsapiens.UCSC.hg19.knownGene. 

# If the annotation file is accessible from AnnotationHub (as is the case for the Ensembl genes), a pre-scanned GTF file can be imported using makeTxDbFromGRanges. 

# Finally, the makeTxDbFromBiomart function can be used to automatically pull a gene

# model from Biomart using biomaRt15.

library("GenomicFeatures")

#We indicate that none of our sequences (chromosomes) are circular using a 0-length character vector

gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))

#The following line produces a GRangesList of all the exons grouped by gene. Each element of the list is a GRanges object of the exons for a gene.
(ebg <- exonsBy(txdb, by="gene"))


# Read counting step ------------------------------------------------------

library(GenomicAlignments) # for single core alingment 
library(BiocParallel) #multiple core alignment windows can't do it!

register(SerialParam())

# The following call creates the SummarizedExperiment 
# object with counts:
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

# The mode argument describes what kind of read overlaps will be counted

# Setting singleEnd to FALSE indicates that the experiment produced paired-end reads, and we want to count a pair of reads (a fragment) only once toward the count for a gene.

# ignore.strand to TRUE
# The fragments argument can be used when singleEnd=FALSE to specify if unpaired reads should be counted (yes if fragments=TRUE).

??GenomicAlignments # see Counting reads with summarizeOverlaps
# vignette


# SummarizedExperiment ----------------------------------------------------

se
assayNames(se)
head(assay(se), 3)
colSums(assay(se))

#The rowRanges, when printed, only shows the first GRanges, and tells us there are 19 more elements:
rowRanges(se)

# The rowRanges also contains metadata about the construction of the gene model in the metadata slot. Here we use a helpful R function, str, to display the metadata compactly:
  str(metadata(rowRanges(se)))
  colData(se)  
  # The colData slot, so far empty, should contain all the metadata. Because we used a column of sampleTable to produce the bamfiles vector, we know the columns of se are in the same order as the rows of sampleTable.
  # We can assign the sampleTable as the colData of the summarized experiment, by converting it into a DataFrame and using the assignment function:
    (colData(se) <- DataFrame(sampleTable))

# The DESeqDataSet, sample information, and the design formula ------------

# The simplest design formula for differential expression would be ~ condition, where condition is a column in colData(dds) that specifies which of two (or more groups) the samples belong to. 
# For the airway experiment, we will specify ~ cell + dex meaning that we want to test for the effect of dexamethasone (dex) controlling for the effect of different cell line (cell). 
# We can see each of the columns just using the $ directly on the SummarizedExperiment or DESeqDataSet: 
  
se$cell
se$dex

# it is prefered in R that the first level of a factor be the reference level (e.g. control, or untreated samples), so we can relevel the dex factor like so:

se$dex <- relevel(se$dex, "untrt")
se$dex
#For running DESeq2 models, you can use R's formula notation to express any fixed-effects experimental design
# Note that DESeq2 uses the same formula notation as, for instance, the lm function of base R.
??results
# If the research aim is to determine for which genes the effect of treatment is different across groups, then interaction terms can be included and tested using a design such as ~ group + treatment +
#   group:treatment. 
# See the manual page for ?results for more examples

# For a full example of using the HTSeq Python package for read counting, please see the pasilla vignette.
# For an example of generating the DESeqDataSet from files produced by htseq-count, please see the DESeq2 vignette.

# In the following sections, we will demonstrate the construction of the DESeqDataSet from two starting points:
# . from a SummarizedExperiment object
# . from a count matrix and a sample information table


# Starting from SummarizedExperiment --------------------------------------


