## Overview
This repository contains a few examples of code I have written for projects during my course in computational biology. 
While the code is my own, each example here is either a product of collaboration with others in the course, or exists 
as part of a larger group project. 

### Mutational Landscape Characterisation & Survivial Analysis (genomics_1.R)
My contribution to a group-based project aimed at replicating results from "The mutational landscape of the SCAN-B 
real-world primary breast cancer transcriptome" (Brueffer et al.). It includes code for recreating the paper's waterfall 
plot characterising mutations by mutation type, sample type and subtypes, and relative frequency, in addition to 
replicating the survival analysis. 

At the bottom of the script, there is also code for my individual extension where I characterised the co-mutational 
landscape of the dataset, and performed a subsequent survivial analysis based on co-mutational status of select gene-
pairs. 

For the most part, this code is self-contained. Where data was readily available on GEO, the GEOquery package is used. 
Any other data is obtained from the data exploration website provided by the reference paper (https://oncogenomics.bmc.lu.se/MutationExplorer/). 

### Differential Expression Analysis (genomics_2.R)
This file contains code for comparing the differential expression workflows in edgeR and DESeq2. The data is un-normalised RNA-seq 
counts from https://bowtie-bio.sourceforge.net/recount/countTables/modencodefly_count_table.txt. 


### Gene Annotation, Bash to R Interfacing (genomics_3.R)
This file contains code for fully-automating gene annotation on the genomes of three drosophila species using an outdated 
(course-required) annotation programme Genscan. While this is not the most technical work I've done, 
I include it here becuase it demonstrates my approach to large-scale data processing and is an example of my flexibility in 
interfacing between bash and R script. The end goal of the functions in this script is to get a single .gff file from the scaffolds of a drosophila genome. 
Since this was performed on three species, the procedure had to be easily reproducible and verifiable. The functions 
defined here allow us to run genscan on a directory of scaffolds using 'genscanner()', validate it has been performed on 
all scaffolds using 'genchecker()', and then convert these outputs to a .gff file using 'gffinator()'. The 
function 'filter()', 'extract_scafid()', and 'genpost()' are all called within the three main functions. 

