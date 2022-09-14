#' @rdname Scoring
#' @title Metabolic Pathway Analysis of RNAseq Data
#' @param DEGpath character, the path to a txt or csv DEG file
#' @param genename character, column name with HUGO Gene Names in DEG file
#' @param sampsize numeric, the sample size of the experiment to be analyzed
#' @param headers character vector of length2 in the form c(log fold change
#' col name, adjusted p value col name)
#' @return pathwayAnalysis() returns a dataframe of pathway scores and pvals
#' @importFrom stringr str_sub
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#' @importFrom utils read.csv read.table
#' @export
#' @examples
#' set.seed(1234)
#' scores <- pathwayAnalysis(
<<<<<<< HEAD
#'                  DEGpath = system.file("extdata/BRCA_DEGS.csv",
#'                                           package = "MetaPhOR"),
#'                  genename = "X",
#'                  sampsize = 1095,
#'                  headers = c("logFC", "adj.P.Val"))
=======
#' system.file("extdata/BRCA_DEGS.csv", package = "MetaPhOR"), "X", 1095,
#' headers = c("logFC", "adj.P.Val"))
>>>>>>> upstream/master
#' View(scores)
pathwayAnalysis <- function(DEGpath, genename, sampsize,
                            headers = c("log2FoldChange", "padj")){
    stopifnot(is.character(DEGpath), length(DEGpath) == 1, !is.na(DEGpath),
<<<<<<< HEAD
                is.character(genename), length(genename) == 1, !is.na(genename),
                is.numeric(sampsize), length(sampsize) == 1, !is.na(sampsize),
                is.vector(headers), length(headers) == 2, !is.na(headers))

    if (str_sub(DEGpath, -3, -1) == "txt"){
        MP_DEGS <- read.table(DEGpath, header = TRUE, sep = "\t")
    } else if (str_sub(DEGpath, -3, -1) == "csv"){
        MP_DEGS <- read.csv(DEGpath, header = TRUE)
    } else {stop("file is not .csv or .txt")}

    MP_DEGS[,headers[2]] <- replace(MP_DEGS[,headers[2]],
                                    which(MP_DEGS[,headers[2]] == 0), 1e-314)
    MP_DEGS[,headers[2]] <- replace(MP_DEGS[,headers[2]],
                                    which(is.na(MP_DEGS[,headers[2]])), 1)
    MP_DEGS[,headers[1]] <- replace(MP_DEGS[,headers[1]],
                                    which(is.na(MP_DEGS[,headers[1]])), 0)
=======
        is.character(genename), length(genename) == 1, !is.na(genename),
        is.numeric(sampsize), length(sampsize) == 1, !is.na(sampsize),
        is.vector(headers), length(headers) == 2, !is.na(headers))

    if (str_sub(DEGpath, -3, -1) == "txt"){#Read in DEG File
        MP_DEGS <- read.table(DEGpath, header = TRUE, sep = "\t")
    } else if (str_sub(DEGpath, -3, -1) == "csv") {MP_DEGS <- read.csv(DEGpath,
        header = TRUE)} else {stop("file is not .csv or .txt")}

    MP_DEGS[,headers[2]] <- replace(MP_DEGS[,headers[2]],
                                which(MP_DEGS[,headers[2]] == 0), 1e-314)
    MP_DEGS[,headers[2]] <- replace(MP_DEGS[,headers[2]],
                                which(is.na(MP_DEGS[,headers[2]])), 1)
    MP_DEGS[,headers[1]] <- replace(MP_DEGS[,headers[1]],
                                which(is.na(MP_DEGS[,headers[1]])), 0)
>>>>>>> upstream/master

    #Calculate Scores
    MP_Scores <- as.data.frame(matrix(nrow = nrow(MP_DEGS), ncol = 2,
                    dimnames = list(rownames(MP_DEGS), c("Score", "ABSScore"))))
    MP_Scores$Score <- (-log(MP_DEGS[,headers[2]]))*MP_DEGS[,headers[1]]
    MP_Scores$ABSScore <- abs(MP_Scores$Score)

    #Read in KEGG Pathways
    AllPathwaysKEGGPS <- data.frame(lapply(KEGGdata, as.character),
                                    stringsAsFactors = FALSE)

    KEGG_Scores <- c() #Calculate Pathway Scores
    ABSScores <- c()
    for (i in seq_len(ncol(AllPathwaysKEGGPS))){
        for (j in c("Score", "ABSScore")){
            y <- MP_Scores[toupper(MP_DEGS[,genename]) %in%
                            as.character(AllPathwaysKEGGPS[,i]),]
            PathwayScore <- as.vector(sum(y[,j])/sqrt(sampsize))
            if (j == "Score"){KEGG_Scores[i] <- PathwayScore
            } else {ABSScores[i] <- PathwayScore}}}

    z <- c() #Bootstrapping for Significance
    pvals <- c()
    abspvals <- c()
    for (i in c("Score", "ABSScore")){
        for (j in seq_len(ncol(AllPathwaysKEGGPS))){
            for (k in seq_len(100000)){boot.Scores <- sample(MP_Scores[,i],
                size = length(grep(".", AllPathwaysKEGGPS[,j])), replace = TRUE)
                #Sum Scores
                z[k] <- sum(boot.Scores)/sqrt(sampsize)}
            if (i == "Score") {pvals[j] <- sum(z >= KEGG_Scores[j])/100000}
            else {abspvals[j] <- sum(z >= ABSScores[j])/100000}}}

    PipelineScores <- as.data.frame(matrix(data =    #Compile Data Frame
                    c(KEGG_Scores, ABSScores, pvals, abspvals),
                    nrow = ncol(AllPathwaysKEGGPS), ncol = 4, byrow = FALSE,
                    dimnames = list(colnames(AllPathwaysKEGGPS),
                    c("Scores", "ABSScores", "ScorePvals", "ABSScorePvals"))))
    return(PipelineScores)}
