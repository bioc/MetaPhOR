##preparing internal data

#load in pathways from file
files <- list.files("/Users/em48662/Desktop/Metabolic_Pipeline/wikipathways-20220510-gpml-Homo_sapiens/")
wpid2names <- c()
wpid2genes <- c()

for (f in seq_along(files)){
    file <- paste0("/Users/em48662/Desktop/Metabolic_Pipeline/wikipathways-20220510-gpml-Homo_sapiens/", files[[f]])
    filename <- utils::tail(unlist(strsplit(file, "/")), n = 1)
    wikipId <- utils::tail(unlist(strsplit(filename, "_")), n = 2)[[1]]

    #read in the gpml file
    doc <- XML::xmlTreeParse(file)
    gpml <- XML::xmlToList(doc)
    nms <- names(gpml)

    #get genes
    tgt_nodes <- gpml[nms == "DataNode"]
    geneprods <- as.vector(unlist(lapply(tgt_nodes, function(x) x$Xref[["Database"]])))
    transcripts <- tgt_nodes[which(geneprods == "Entrez Gene" | geneprods == "Ensembl")]
    genes <- as.vector(unlist(lapply(transcripts, function(x) x$.attrs[["TextLabel"]])))

    #set up mapping
    wpid2names <- c(wpid2names, wikipId, gpml$.attrs[["Name"]])
    if (length(genes) > 0){
        for (i in seq_along(genes)){
            wpid2genes <- c(wpid2genes, wikipId, genes[i])
        }
    }
}

wpid2name <- as.data.frame(matrix(wpid2names, nrow = length(wpid2names)/2, ncol = 2, byrow = T))
colnames(wpid2name) <- c("WPID", "name")

wpid2gene <- as.data.frame(matrix(wpid2genes, nrow = length(wpid2genes)/2, ncol = 2, byrow = T))
colnames(wpid2gene) <- c("WPID", "gene")

wpid2name <- wpid2name[which(!unique(wpid2name$WPID) %in% unique(wpid2gene$WPID) == F),]

#KEGG Data
KEGGdata <- read.delim("/Users/em48662/Desktop/Metabolic_Pipeline/AllPathwaysKEGGPS.txt")
KEGGdata <- KEGGdata[,-c(grep("Glycolysis.Gluconeogenesis", colnames(KEGGdata)))]

#Saving Data
usethis::use_data(wpid2name, wpid2gene, KEGGdata, internal = T, overwrite = T)

