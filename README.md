# LVEXdetector

## Detecting DNA methylation patterns suggestive of variable escape from X-chromosome inactivation
<img alt="LVEXdetector_logo" src="https://github.com/user-attachments/assets/c2d0a063-db75-42bc-bf12-d2cc743df487" width="420" />



The X chromosome is often excluded from studies analyzing associations between traits and DNA methylation. In females, one copy of most genes on X is inactivated through the establishment of X-chromosome inactivation (XCI), implemented by DNA methylation of the gene promoter on the inactive X, and this leads to challenges in analyzing and interpreting DNA methylation data patterns. Particularly for sex-biased diseases and traits, there may be many loci of interest on chromosome X.  To address the need for appropriate methods for analysis of DNA methylation data on the X chromosome, we developed a statistical approach to infer locus-specific escape from XCI that varies depending on phenotype or covariate values.  Performance of this method is illustrated by analysis of data from two sex-biased traits: rheumatoid arthritis which is 3-fold more common in females, and recurrent venous thromboembolism which occurs 2.5 times more in males. Analyses of these two datasets identify new trait-associated loci on X, demonstrate the capabilities of the new method for both bisulfite sequencing data and Illumina EPIC data, suggest at least one locus where variable escape may explain a sex-specific disease association, and rule out variable escape as a potential explanation at other loci. 

## Installation

Install the development version of **LVEXdetector** from GitHub:

```r
# install devtools if not already installed
install.packages("devtools")

# install LVEXdetector
devtools::install_github("qc-zhao/LVEXdetector")

# load the package
library(LVEXdetector)
