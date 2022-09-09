## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    dpi=200
)

library(kableExtra)

## ---- message=FALSE-----------------------------------------------------------
library(MetaPhOR)

## -----------------------------------------------------------------------------
exdegs <- read.csv(system.file("extdata/exampledegs.csv", package = "MetaPhOR"),
                    header = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
kable(head(exdegs), format="latex", booktabs = TRUE)    %>%
    kable_styling(position = "center", latex_options = c("striped", 
                                                            "hold_position"))

## ---- include=FALSE-----------------------------------------------------------
# BRCA, OVCA, PRAD 
# sampsize <- c(1095, 378, 497)

## -----------------------------------------------------------------------------
set.seed(1234)

brca <- pathwayAnalysis(system.file("extdata/BRCA_DEGS.csv", 
            package = "MetaPhOR"), "X", 1095, headers = c("logFC", "adj.P.Val"))

## ---- echo = FALSE------------------------------------------------------------
kable(head(brca), format="latex", booktabs = TRUE)    %>%
    kable_styling(position = "center", latex_options = c("striped", 
                                                            "hold_position"))

## -----------------------------------------------------------------------------
pval <- bubblePlot(brca, "Pval", .85)
plot(pval)

## -----------------------------------------------------------------------------
logfc <- bubblePlot(brca, "LogFC", .85)
plot(logfc)

## -----------------------------------------------------------------------------
##read in two additional sets of scores,
##run in the same manner as brca for comparison

ovca <- read.csv(system.file("extdata/OVCA_Scores.csv",
                            package = "MetaPhOR"), header = TRUE, row.names = 1)
prad <- read.csv(system.file("extdata/PRAD_Scores.csv",
                            package = "MetaPhOR"), header = TRUE, row.names = 1)

all.scores <- list(brca, ovca, prad)
names <- c("BRCA", "OVCA", "PRAD")
metaHeatmap(all.scores, names, 0.05)

## ---- eval = FALSE------------------------------------------------------------
#  cytoPath("Tryptophan Metabolism", "BRCA_DEGS.csv", paste(getwd(),
#  "BRCA_Tryptophan_Pathway", sep = "/"), "X", headers = c("logFC", "adj.P.Val"))

## ---- echo = F, out.width = "100%"--------------------------------------------
knitr::include_graphics(c(system.file("extdata", "BRCA_Tryptophan_Pathway.png", 
                                        package = "MetaPhOR")))

## ---- eval = FALSE------------------------------------------------------------
#  pathwayList()

## ---- echo = F----------------------------------------------------------------
kable(head(pathwayList()), format = "latex", booktabs = TRUE)    %>%
    kable_styling(position = "center", latex_options = c("striped", 
                                                            "hold_position"))

## -----------------------------------------------------------------------------
sessionInfo()

