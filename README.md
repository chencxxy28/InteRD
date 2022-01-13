# InteRD
Integrated and robust cell-type deconvolution
 
**InteRD** has three primary advantages. First, it is able to effectively integrate deconvolution results from multiple scRNA-seq datasets. Second, **InteRD** calibrates estimates from reference-based deconvolution by taking into account extra biological information as priors. Third, the proposed algorithm is equipped with a data-driven mechanism of self-control designed to be robust to the introduction of inaccurate information in the deconvolution system.

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
