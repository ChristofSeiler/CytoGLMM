# Precompiled vignettes that depend on long computations
# Must manually move image files from CytoGLMM/figure to CytoGLMM/vignettes/figure after knit

knitr::knit("vignettes/pregnancy_dataset.Rmd.orig", "vignettes/pregnancy_dataset.Rmd")
