Independent-eQTLs
=================

My slides from the August 14th, 2014 lab meeting can be found attached to my page in the [lab wiki](https://medwiki.stanford.edu/display/montgomerylab/Finding+Independent+eQTLs).

## Data ##
Sample data is located in the `sample_data` directory. This data is dummy example data adapted from Matrix eQTL's sample data. This dummy data is only meant for getting a sense of the format and not for running any test analysis. For real data, you can use the attachments on my page in the lab wiki which is data from the DGN cohort: [Current projects -> Finding Independent eQTLs](https://medwiki.stanford.edu/display/montgomerylab/Finding+Independent+eQTLs)

## Plots ##

## Scripts ##

### findeQTLs.R ###
```
source("findeQTLs.R")
```
This script, adapted from [http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/), finds cis/trans-eQTLs and stores them in the R environment as `me`. Edit the script to fit your own preferences. Further help can be found at [http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/).

**Note: You may need to edit the script depending on the format of your data.**

### setUpEnvironment.R ###
```
setUpEnvironment(geneExpressionFile, genotypeFile, genelocFile, snpslocFile, me)
```
This script takes in the input and output data of Matrix eQTL and outputs a list of data frames needed for the other scripts. **It assumes your data has the same format (i.e. contains headers and row names) as the sample data, so modify your sample data or the script accordingly.**

##### Input: ######
- **geneExpressionFile**: The location of the gene expression data. Equivalent to `GE.txt` in the sample data.
- **genotypeFile**: The location of the genotype data. Equivalent to `SNP.txt` in the sample data.
- **genelocFile**: The location of the gene location data. Equivalent to `geneloc.txt` in the sample data.
- **snpslocFile**: The location of the SNP location data. Equivalent to `snpsloc.txt` in the sample data.
- **me**: The object that Matrix eQTL outputs. Equivalent to `me` in `setUpEnvironment.R`.

##### Returns: #####
A list containing these data frames *(you can use these data frames in the scripts after this)*:
- **SNP.t**: The transpose of the data from `SNP.txt`.
- **SNP.train**: A subset of `SNP.t` to act as a training set.
- **SNP.test**: A subset of `SNP.t` to act as a test set.
- **GE**: The data from `GE.txt`.
- **GE.train**: A subset of `GE` to act as a training set.
- **GE.test**: A subset of `GE` to act as a test set.
- **geneloc**: The gene location data from `geneloc.txt`.
- **snpsloc**: The SNP location data from `snpsloc.txt`.
- **eQTLs**: A data frame of the cis-eQTLs that Matrix eQTL outputs.
- **besteQTLs**: A data frame of the cis-eQTLs that Matrix eQTL outputs, with only one eQTL (the lowest p-value) per gene.

### testAll.R ###
```
testAll(GE, SNPs, eQTLs, method = "pvalue", direction = "forward", steps = 100, verbose = TRUE)
```
This script finds independent eQTLs for every gene in your data and returns information about these eQTLs.

##### Input: #####
- **GE**: Gene expression data.
- **SNPs**: Corresponds to SNP.t; the genotype data.
- **eQTLs**: A data frame of the cis-eQTLs that Matrix eQTL outputs. 
- **method**: Can take on the value "pvalue" or "AIC", correponding to the criterion by which the algorithm selects SNPs.
- **direction**: Can take on the value "forward" or "both". "forward" corresponds to a unidirectional method and "both" corresponds to a bidirectional method. 
- **steps**: The number of maximum steps the algorithm can take, so it doesn't go on forever in a loop deleting and adding a SNP. The default is 100.
- **verbose**: If `TRUE`, the script will print its progress. Otherwise, the script does not print anything.

##### Returns: #####
A list containing a list and a data frame:
- **indEQTLs**: A data frame with three columns:
..- genes: The gene.
..- numind: The number of independent eQTLs.
..- rsquared: The r-squared value that the SNPs yield for that gene.
- **snps**: A list of lists. The i-th list corresponds to the list of SNPs for the i-th gene in `indeQTLs`.
