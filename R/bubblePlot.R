#' @rdname BubblePlot
#' @title Create a Bubble Plot for Individual Samples
#' @param scorelist dataframe(1) the output of Pathway Analysis fun
#' @param labeltext character(1) what to label points by: LogFC or Pval
#' @param labelsize numeric(1) size of text labels for points
#' @return bubblePlot() returns a bubble plot using pathway scores, pval, logfc
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples
#' brca <- read.csv(system.file("extdata/BRCA_Scores.csv",
#'                     package = "MetaPhOR"),
#'                     header = TRUE,
#'                     row.names = 1)
#'
#' #Bubble Plot Labeled By P Value
#' bubblePlot(scorelist = brca,
#'             labeltext = "Pval",
#'             labelsize = .85)
#'
#' #Bubble Plot Labeled by LogFC
#' bubblePlot(scorelist = brca,
#'             labeltext = "LogFC",
#'             labelsize = .85)
bubblePlot <- function(scorelist, labeltext, labelsize = .25){
    stopifnot(
        is.data.frame(scorelist), length(scorelist) == 4, !is.na(scorelist),
        is.character(labeltext), length(labeltext) == 1, !is.na(labeltext),
        labeltext %in% c("LogFC", "Pval")
    )

    scorelist <- scorelist[,-4]
    scorelist <- na.omit(scorelist)

    #Find Top Scores
    if (labeltext == "LogFC"){
        datalabs <- scorelist[order(-abs(scorelist$Scores))[seq_len(10)],]
    } else if (labeltext == "Pval") {
        datalabs <- scorelist[order(scorelist$ScorePvals)[seq_len(10)],]
    } else {}

    #Plot Bubble Plot
    plot <- ggplot(scorelist,
        aes(x = Scores, y = ABSScores, size = ScorePvals, color = Scores)) +
        geom_point(alpha = 0.7) +
        geom_text_repel(data = datalabs, aes(label = gsub("\\.", " ",
                rownames(datalabs)), size = labelsize), color = "black",
                force = 3, max.overlaps = Inf) +
        scale_size(trans = 'reverse') +
        scale_color_gradientn(colors =
                colorRampPalette(c("darkblue", "grey", "red"))(n = 299),
                limits = c(-max(abs(scorelist$Scores)),
                max(abs(scorelist$Scores)))) +
        theme_minimal()
        return(plot)
}
