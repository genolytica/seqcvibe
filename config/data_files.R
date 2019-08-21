dataFiles <- list(
    GEO=list(
        GSE59017="data/GSE59017/GSE59017.rda",
        GSE63485="data/GSE63485/GSE63485.rda"
    ),
    ENA=list(
        PRJEB16733="data/PRJEB16733/PRJEB16733.rda",
        PRJEB65039="data/PRJEB65039/PRJEB65039.rda"
    )
)

#config=read.delim("./config/metadata.txt")
#rownames(config) <- as.character(config$sample_id)
#buildTrackList(
#   config=config,
#   annoPath="/media/raid/tracks/hybridsuite/reference",
#   appBase="/media/raid/software/bigseqcvis",
#   urlBase="http://epigenomics.fleming.gr/tracks"
#)
