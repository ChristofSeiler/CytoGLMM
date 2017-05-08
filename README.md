# CytoGLMM

## R Package for Differential Analysis of Mass Cytometry Experiments

This packge uses generalized linear mixed model to relate ~40 protein marker expressions (the explanatory variables) to a simulation condition (the response variable, e.g. infection of vaccination). Both maximum likelihood and full Bayesian posterior (using Stan) implementations are available. The key point is to account for donor-specific variablity. Estimating donor-specific parameters is straightforward because the number of cells exceed the number of markers, whereas estimating population-level parameters is challenging because the number of markers usually exceed the number of donors.

## Installation

The R package is available from github and can be installed by running the following command in R:

```
install.packages("devtools")
devtools::install_github("ChristofSeiler/CytoGLMM",build_vignettes = TRUE)
```

## Getting Started

Read the vignette for a step-by-step example workflow:

```
vignette("workflow", package = "CytoGLMM")
```
