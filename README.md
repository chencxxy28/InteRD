# InteRD
Integrated and robust cell-type deconvolution
 
**InteRD** unifies multiple deconvolution schemes to infer cell type proportions from the target bulk RNA-seq data. Three unique features are embraced in this algorithm: first, **InteRD** is able to incorporate extra biological information from external data sources enables (e.g., scRNA_seq and other bulk RNA_seq data from independent studies); second, **InteRD** calibrates the reference-free algorithm by taking into account the proportion estimates from a reference-based approach; third, **InteRD** is robust to incorrect information from any of the provided data sources.

# Installation
You can install the released version of InteRD with:
```
#install devtools if necessary
install.packages("devtools")

#install the InteRD package
devtools::install_github("chencxxy28/InteRD")

#load
library(InteRD)
```

If [_Biobase_](https://bioconductor.org/packages/release/bioc/html/Biobase.html) package is not available, please install it first before installation of **InteRD**
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase")
```

# Vignettes
Please visit [Tutorial](updating).
