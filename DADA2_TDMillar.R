

### Trevor D. Millar
### Project Start Mar 17, 2021
### This is a code through of the DADA2 walkthrough 
### https://benjjneb.github.io/dada2/tutorial.html


### Here is how you install the package
#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16")

library(Biostrings) ## Loading this up before dada2 was the trick to getting everything working. 
library(dada2)   ## Load up dada2
library(ggplot2) ## We need ggplot2 for the plotQualityProfile function
library(tidyverse)


path <- "./DADA2_Data/MiSeq_SOP" ## Define path as the relative file path to the data
list.files(path) ## Check the list of files to make sure my file path works


# Forward and reverse fastq filenames have format: 
# SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
## Forward Sequences are sorted list of files in path that have "_R1_001.fastq, we use the full names
print(fnFs)
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
print(fnRs)

sample.names <- sapply(strsplit(basename(fnFs), "_"), '[', 1) ## Delimiter = _ , take 1st item 
print(sample.names)

fnFs[1:2] # This is the first two elements of the fnFs object
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
## Create a list of file paths from the sample names (Paste sample name and ending)
## "filtered" just creates a sub-directory within the WD that is 'path'
## This directory is really created when we run these data through filterAndTrim below
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
## Create a list of file paths from the sample names
## "filtered" just creates a sub-directory within the WD that is 'path'
## This directory is really created when we run these data through filterAndTrim below

names(filtFs) <- sample.names
## Assign the sample names of the filtFs object 
names(filtRs) <- sample.names
## Assign the sample names of the filtFs object 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # multithreading blows up my computer
## The filtered output lands in the filtered directory we created before. 

head(out) ## This object is just a summary of how many reads we fed into filterAndTrim and 
## how many reads we got out. The output that we are most interested in is in the filtered directory


################################################################################
##### LEARN THE ERROR RATES
################################################################################

errF <- learnErrors(filtFs, multithread = FALSE)
### This is using a fancy algorithm to find the errors in the dataset based on 
errR <- learnErrors(filtRs, multithread = FALSE)


plotErrors(errF, nominalQ = TRUE) ## I need to look more into this to understand what is happening here
## The black line is the estimated error rate from the machine learning algoritm in learnErrors
## The red line is the error rates that are expected based on the Q-score alone. (Phred Quality score)

################################################################################
##### SAMPLE INTERFERENCE
################################################################################

## the DADA algorithm uses the dereplicated amplicon seq reads and the errors learned from the machine 
## learning algorithm to remove all sequencing errors in the input file. 

dadaFs <- dada(filtFs, err=errF, multithread = FALSE)
dadaRs <- dada(filtRs, err=errR, multithread = FALSE)

dadaFs[[1]] ## This is the first object created by the dada algorithm. 
## What was before thought to be 1979 unique sequences was found to only be 
## 128 unique sequences once the algorithm worked on the data. (True sequence variant)


################################################################################
##### MERGE PAIRED READS
################################################################################

## We are going to merge the forward and reverse reads here. 
## The reads must have overlap, and we will used the merged reads to create our contigs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
## for the first file, it looks like there are 106 unique pairings, 
## does that mean we have 106 unique contigs?

################################################################################
##### CONSTRUCT SEQUENCE TABLE
################################################################################

## Here we are going to make a sequence table from our merged dataset
## This sequence table has 20 rows (one for each sample), and 293 columns, one for each
## amplicon sequence variant (ASV). The numbers in the table indicate the count of ASVs in each sample. 

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))
getSequences(seqtab)


## This shows you the 293 sequences in our seqtab table. 
##nchar is a function that counts the number of characters in the character parameter
## Table is making a table of the different lengths of the denoised contigs. 


## I want to take a better look at seq tab
df_seqtab <- data.frame(seqtab)
sum(df_seqtab[,1]) ## THis is the sum of the number of seqs in the first file; F3D0
view(seqtab)

################################################################################
##### REMOVE CHIMERAS
################################################################################

## A chimera in this context is a PCR amplicon whose extension is aborted. Say 
## the polymerase falls off mid reaction. The aborted section will be a primer in the
## next round of amplification. If it anneals to the wrong template you end up with half
## of the read comming from one source, and the other half comming from another source. 

## Chimeras need to be identified, as they are not true sequence variants, rather 
## a sequence that is product of errors during PCR amplification. 

## DADA2 has a function that can identify the two parts of the chimera as being present 
## in the more abundant true sequence variants and thus remove the chimeras from the dataset. 

seqtab_nochimeras <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose = TRUE)
## 61 bimeras (chimeras) were identified out of a total of 293 input sequences. 

dim(seqtab_nochimeras)

sum(seqtab_nochimeras) ## total number of reads in the no chimeras dataset
sum(seqtab) ## total number of reads in the dataset, including chimeras

## lets take the quotient of the two to see the percent of real reads. 

sum(seqtab_nochimeras) / sum(seqtab)
## It looks like about 3.5% of our reads were actually chimeras. 


################################################################################
##### TRACK READS THROUGH THE PIPELINE. 
################################################################################

## Each step of our pipeline should be removing some reads for one reason or another. 
## We can track how many reads are input and output of every step of our pipeline up to this point


getN <- function(x) sum(getUniques(x))
## getUniques is going to return a vector named by the sequence and valued by the abundance of that sequence. 
## when we take the sum of those values, we are getting the total number of seqs. 

## cbind is column bind. Bind by column if you will. 

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_nochimeras))


## It's hard to understand what the rows and columns mean in track, lets add the row and column names

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names


################################################################################
##### ASSIGN TAXONOMY
################################################################################

## Just as the name says, we have some nice reads, now lets have DADA2 
## help us assign taxonomy. 
## Reading up on this makes me think its a machine learning algorithm and needs some "training"
## to get the right answer. 
## In order to train the algorithm, we are going to use a data set with known seqs and taxonomy
## The walkthrough suggests we use data from this web page: https://zenodo.org/record/4587955#.YFPnaPtKhhE
## I'm personally using v138.1, but it looks like the walkthrough used another set of data. 


taxa <- assignTaxonomy(seqtab_nochimeras, "./DADA2_Data/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)







