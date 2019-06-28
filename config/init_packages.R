# Function initializing SeqCVIBE universe (packages, persistent variables, etc.)
initPackages <- function(session) {
    # Initial page loading indicator, until all content is loaded
    ftProgress <- shiny::Progress$new()
    ftProgress$initialize(session,min=0,max=8)
    ftProgress$set(message="Starting:",value=0)
    on.exit(ftProgress$close())
    # Progress update function
    updateFtProgress <- function(value=NULL,detail=NULL) {
        if (is.null(value)) {
            value <- ftProgress$getValue()
            value <- value + 1
        }
        ftProgress$set(value=value,detail=detail)
    }
    
    # Load packages
    updateFtProgress(value=1,detail="Loading DT")
    require(DT)
    updateFtProgress(value=2,detail="Loading GenomicRanges")
    require(GenomicRanges)
    updateFtProgress(value=3,detail="Loading GenomicAlignments")
    require(GenomicAlignments)
    updateFtProgress(value=4,detail="Loading ggplot2")
    require(ggplot2)
    updateFtProgress(value=5,detail="Loading rtracklayer")
    require(rtracklayer)
    updateFtProgress(value=6,detail="Loading ggbio")
    require(ggbio)
    updateFtProgress(value=7,detail="Loading metaseqR")
    require(GenomicRanges)
    updateFtProgress(value=8,detail="Loading d3heatmap")
    require(d3heatmap)
    #require(plotly)
}
