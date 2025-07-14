## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.retina = 2, dpi = 96, fig.width = 7.2916667, fig.asp = 0.6178571)

## ----simulated_data-----------------------------------------------------------
library("CytoGLMM")
set.seed(23)
df <- generate_data()
df[1:5,1:5]

## ----protein_names------------------------------------------------------------
protein_names <- names(df)[3:12]

## ----transform----------------------------------------------------------------
df <- dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))

## ----glm_fit------------------------------------------------------------------
glm_fit <- CytoGLMM::cytoglm(df,
                             protein_names = protein_names,
                             condition = "condition",
                             group = "donor",
                             num_cores = 1,
                             num_boot = 1000)
glm_fit

## ----glm_plot-----------------------------------------------------------------
plot(glm_fit)

## ----glm_summarize------------------------------------------------------------
summary(glm_fit)

## ----glm_p_values-------------------------------------------------------------
summary(glm_fit) |> dplyr::filter(pvalues_adj < 0.1)

## ----session_info-------------------------------------------------------------
sessionInfo()

