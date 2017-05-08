# CytoGLMM

## R Package for Differential Analysis of Mass Cytometry Experiments

This packge uses generalized linear mixed model to relate simulation condition (the response variable, e.g. infection of vaccination) to ~40 expression of protein markers (the explanatory variables). Both maximum likelihood and full Bayesian posterior (using Stan) are available. The key point is to account for donor-specific variablity (and in the Bayesian model to borrow information across donor through partial pooling). Estimating donor-specific parameters is straightforward because the number of cells exceed the number of markers, whereas estimating population-level parameters is challenging because the number of markers usually exceed the number of donors.
