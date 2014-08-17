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

**Note: You may need to edit the script depending on what format your data is in.**

### setUpEnvironment.R ###
```
setUpEnvironment(geneExpressionFile, genotypeFile, genelocFile, snpslocFile, me)
```
This script takes in the input and output data of Matrix eQTL and outputs a list of data frames needed for the other scripts. **It assumes your data has the same format (i.e. contains headers and row names) as the sample data, so modify your sample data or the script accordingly.**

##### Input ######
- **geneExpressionFile**: The location of the gene expression data. Equivalent to `GE.txt` in the sample data.
- **genotypeFile**: The location of the genotype data. Equivalent to `SNP.txt` in the sample data.
- **genelocFile**: The location of the gene location data. Equivalent to `geneloc.txt` in the sample data.
- **snpslocFile**: The location of the SNP location data. Equivalent to `snpsloc.txt` in the sample data.
- **me**: The object that Matrix eQTL outputs. Equivalent to `me` in `setUpEnvironment.R`.

##### Output #####
A list containing these data frames:
- **SNP.t**: The transpose of `SNP.txt`.
- **SNP.train**: A subset of `SNP.t` to act as a training set.
- **SNP.test**: A subset of `SNP.t` to act as a test set.
- **GE**: The data from `GE.txt`.
- **GE.train**: A subset of `GE` to act as a training set.
- **GE.test**: A subset of `GE` to act as a test set.


## TODO ##
- Add more text to slides for context
- Add plots and scripts
