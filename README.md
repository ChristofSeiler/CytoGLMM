# CytoGLMM: R Package for Differential Analysis of Mass Cytometry Experiment

We use generalized linear model where the response variable is the outcome and the explanatory variables are the ~40 protein expression profiles. We build a hierarchical model of CyTOF data that allows us to estimate population level parameters and marginalize out the donor-specific parameters. Estimating donor- specific parameters is straightforward because the number of cells exceed the number of markers, whereas estimating population-level parameters is challenging because the number of markers usually exceed the number of donors. Therefore it is important to borrow information through partial pooling across donors when estimating the population-level parameters. 

Both maximum likelihood and Bayesian implementation (with Stan) are available.
