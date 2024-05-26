## Overview
This repository contains a few examples of code I have written for projects during my course in computational biology. 
While the code is my own, each example here is either a product of collaboration with others in the course, or exists 
as part of a larger project worked on with others. 

### Genomics_1
My contribution to a group-based project aimed at replicating results from "The mutational landscape of the SCAN-B 
real-world primary breast cancer transcriptome" (Brueffer et al.). It includes code for recreating the paper's waterfall 
plot characterising mutations by mutation type, sample type and subtypes, and relative frequency, in addition to 
replicating the survival analysis. 

At the bottom of the script, there is also code for my individual extension where I characterised the co-mutational 
landscape of the dataset, and performed a subsequent survivial analysis based on co-mutational status of select gene-
pairs. 

For the most part, this code is self-contained. Where data was readily available on GEO, the GEOquery package is used. 
Any other data is obtained from the data exploration website provided by the reference paper (https://oncogenomics.bmc.lu.se/MutationExplorer/). 

### 
