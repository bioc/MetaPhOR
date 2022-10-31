#' @rdname MetaHeatmap
#' @title Create a Heatmap for Comparing Multiple Samples
#' @param scorelist list of outputs from pathwayAnalysis()
#' @param samplenames vector of samples names for axis labels
#' @param pvalcut numeric, the p val over which pathways will not be included
#' @return metaHeatmap() returns a heatmap of significant dysregulated pathways
#' @return for each sample included
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
#' @examples
#' brca <- read.csv(system.file("extdata/BRCA_Scores.csv",
#'             package = "MetaPhOR"), header = TRUE, row.names = 1)
#'
#' ovca <- read.csv(system.file("extdata/OVCA_Scores.csv",
#'             package = "MetaPhOR"), header = TRUE, row.names = 1)
#'
#' prad <- read.csv(system.file("extdata/PRAD_Scores.csv",
#'             package = "MetaPhOR"), header = TRUE, row.names = 1)
#'
#' all.scores <- list(brca, ovca, prad)
#' names <- c("BRCA", "OVCA", "PRAD")
#'
#' metaHeatmap(scorelist = all.scores,
#'             samplenames = names,
#'             pvalcut = 0.05)
metaHeatmap <- function(scorelist, samplenames,  pvalcut = 0.05){
    stopifnot(
        is.list(scorelist), !is.na(scorelist),
        lapply(seq_along(scorelist), function(i)
            is.data.frame(scorelist[[i]])) == TRUE,
        is.vector(samplenames), !is.na(samplenames),
        lapply(seq_along(samplenames), function(i)
            is.character(samplenames[[i]])) == TRUE,
        is.numeric(pvalcut), length(pvalcut) == 1, !is.na(pvalcut),
        length(scorelist) == length(samplenames)
    )

    #Subset by Score Type and Pval Cutoff
    matdata <- unlist(lapply(scorelist, function(x)
        replace(x$ABSScores, which(x$ABSScorePvals >= pvalcut), 0)))

    #Make Heatmap Matrix
    hmat <- matrix(data = matdata, byrow = FALSE,
                    nrow = 114,
                    ncol = length(samplenames),
                    dimnames = list(gsub("\\.", " ",
                    rownames(scorelist[[1]])), samplenames))

    #remove all zero rows
    hmat[is.na(hmat)] <- 0
    hmat <- hmat[rowSums(hmat[]) > 0, , drop = FALSE]

    #Scale by Column
    hmat <- scale(hmat, center = FALSE)

    #Plot Heatmap
    plot <- pheatmap(hmat, colorRampPalette(c("grey", "red"))(n = 299),
                    fontsize_col = 8, fontsize_row = 10, border_color = NA)
    return(plot)
}
