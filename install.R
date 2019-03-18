pkgs_needed = c("devtools","tidyverse","magrittr","FlowRepositoryR",
                "flowCore","openCyto","scales","parallel",
                "RColorBrewer","ggcorrplot","SummarizedExperiment",
                "lme4","lmerTest","knitr")
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs_needed)
devtools::install_github("ChristofSeiler/CytoGLMM")
devtools::install_github("RGLab/ggcyto", ref="trunk")
