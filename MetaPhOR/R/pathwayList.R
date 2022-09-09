#' @rdname Pathways
#' @title List Available Metabolic rWikiPathways
#' @return pathwayList() returns a list of rWikiPathways for use in CytoPath()
#' @export
#' @examples
#' pathwayList()
pathwayList <- function(){
        as.data.frame(wpid2name$name)
}
