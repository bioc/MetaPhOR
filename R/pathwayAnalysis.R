#' @rdname Scoring
#' @title Metabolic Pathway Analysis of RNAseq Data
#' @param DEGpath character, the path to a txt or csv DEG file
#' @param genename character, column name with HUGO Gene Names in DEG file
#' @param sampsize numeric, the sample size of the experiment to be analyzed
#' @param iters numeric, the number of iterations of resampling to perform in
#' bootstrapping
#' @param headers character vector of length2 in the form c(log fold change
#' col name, adjusted p value col name)
#' @return pathwayAnalysis() returns a dataframe of pathway scores and pvals
#' @importFrom stringr str_sub
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#' @importFrom utils read.csv read.table
#' @export
#' @examples
#' #iterations (iters) of resampling in bootstraping set to 50,000 for speed
#' #100,000 iterations recommended for improved power
#'
#' set.seed(1234)
#'
#' scores <- pathwayAnalysis(
#'                 DEGpath = system.file("extdata/BRCA_DEGS.csv",
#'                                         package = "MetaPhOR"),
#'                 genename = "X",
#'                 sampsize = 1095,
#'                 iters = 50000,
#'                 headers = c("logFC", "adj.P.Val"))
#' scores
pathwayAnalysis <- function(DEGpath, genename, sampsize, iters = 100000,
                            headers = c("log2FoldChange", "padj")){
    stopifnot(is.character(DEGpath), length(DEGpath) == 1, !is.na(DEGpath),
                is.character(genename), length(genename) == 1, !is.na(genename),
                is.numeric(sampsize), length(sampsize) == 1, !is.na(sampsize),
                is.vector(headers), length(headers) == 2, !is.na(headers),
                is.numeric(iters), length(iters) == 1, !is.na(iters))

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

    MP_Scores <- as.data.frame(matrix(nrow = nrow(MP_DEGS), ncol = 2,
                    dimnames = list(rownames(MP_DEGS), c("Score", "ABSScore"))))
    MP_Scores$Score <- (-log(MP_DEGS[,headers[2]]))*MP_DEGS[,headers[1]]
    MP_Scores$ABSScore <- abs(MP_Scores$Score) #Calculate Scores

    AllPathwaysKEGGPS <- data.frame(lapply(KEGGdata, as.character),
                                    stringsAsFactors = FALSE) #Read KEGG Paths

    y <- apply(AllPathwaysKEGGPS, 2, function(x) #Calculate Pathway Scores
        MP_Scores[toupper(MP_DEGS[,genename]) %in% as.character(x),])
    KEGG_Scores <- as.vector(unlist(lapply(y, function(x)
        as.vector(sum(x[,"Score"])/sqrt(sampsize)))))
    ABSScores <- as.vector(unlist(lapply(y, function(x)
        as.vector(sum(x[,"ABSScore"])/sqrt(sampsize)))))

    boot.Scores <- replicate(iters, apply(AllPathwaysKEGGPS, 2, #Bootstrapping
        function(x) sample(MP_Scores[,"Score"], size = length(grep(".", x)),
        replace = TRUE)))
    foo <- lapply(boot.Scores, function(x) sum(unlist(x))/sqrt(sampsize))
    bar <- matrix(foo, nrow = 114, ncol = iters,
        dimnames = list(rownames(boot.Scores, colnames(boot.Scores))))
    pvals <- unlist(lapply(seq_len(ncol(AllPathwaysKEGGPS)),
        function(x) sum(bar[x,] >= KEGG_Scores[x])/iters))

    boot.Scores <- replicate(iters, apply(AllPathwaysKEGGPS, 2,#Bootstrapping
        function(x) sample(MP_Scores[,"ABSScore"], size = length(grep(".", x)),
        replace = TRUE)))
    foo <- lapply(boot.Scores, function(x) sum(unlist(x))/sqrt(sampsize))
    bar <- matrix(foo, nrow = 114, ncol = iters,
        dimnames = list(rownames(boot.Scores, colnames(boot.Scores))))
    abspvals <- unlist(lapply(seq_len(ncol(AllPathwaysKEGGPS)),
        function(x) sum(bar[x,] >= ABSScores[x])/iters))

    PipelineScores <- as.data.frame(matrix(data =    #Compile Data Frame
                    c(KEGG_Scores, ABSScores, pvals, abspvals),
                    nrow = ncol(AllPathwaysKEGGPS), ncol = 4, byrow = FALSE,
                    dimnames = list(colnames(AllPathwaysKEGGPS),
                    c("Scores", "ABSScores", "ScorePvals", "ABSScorePvals"))))
    return(PipelineScores)}
