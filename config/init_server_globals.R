# Intit globals (don't want to laod these everytime a use connects and they 
# don't take much time to load)

# Load basic packages
require(shiny)
require(shinyjs)
require(colourpicker)
require(RSQLite)

# Load additional functions
source("lib/control.R")
source("lib/util.R")

# Load metadata SQLite
metadata <- dbConnect(drv=RSQLite::SQLite(),dbname='data/metadata.sqlite')

# Load Bookmarks SQLite
bookmarks <- dbConnect(drv=RSQLite::SQLite(),dbname='data/bookmarks.sqlite')

# Intialize metadata reactive content
sources <- as.character(dbGetQuery(metadata, "SELECT DISTINCT(source) FROM metadata")$source)
datasets <- as.character(dbGetQuery(metadata, paste0("SELECT DISTINCT(dataset) FROM metadata WHERE source == '",sources[1],"'"))$dataset)
classes <- as.character(dbGetQuery(metadata, paste0("SELECT DISTINCT(class) FROM metadata WHERE source == '",sources[1],"' AND dataset == '",datasets[1],"'"))$class)
genomes <- as.character(dbGetQuery(metadata, "SELECT DISTINCT(genome) FROM metadata")$genome)
genome <- genomes[1]
title <- as.character(dbGetQuery(metadata, paste0("SELECT DISTINCT(title) FROM summaries WHERE dataset == '",datasets[1],"'"))$title)
link <- as.character(dbGetQuery(metadata, paste0("SELECT DISTINCT(link) FROM summaries WHERE dataset == '",datasets[1],"'"))$link)
short_summary <- as.character(dbGetQuery(metadata, paste0("SELECT DISTINCT(short_summary) FROM summaries WHERE dataset == '",datasets[1],"'"))$short_summary)

# Load data file hash
load("data/dataFiles.RData")
# source("config/data_files.R")
allClasses <- as.character(dbGetQuery(metadata, "SELECT DISTINCT(class) FROM metadata")$class)
baseColours <- c("#B40000","#00B400","#0000B4","#B45200","#9B59B6","#21BCBF",
    "#BC4800","#135C34","#838F00","#4900B5")
baseColours <- rep(baseColours,length.out=length(allClasses))
names(baseColours) <- allClasses

# Keep track of loaded annotations. The first will be loaded when needed.
#load(file.path("genome",genome,"gene.rda"))
loadedGenomes <- vector("list",length(genomes))
names(loadedGenomes) <- genomes
for (gen in genomes) {
    loadedGenomes[[gen]] <- list(
        geneNames=NULL,
        dbGene=NULL,
        dbExon=NULL
    )
}
## Old chunk when we were loading the first genome at app start. Will be removed
## after a while.
#geneNames <- names(gene)
#names(geneNames) <- as.character(elementMetadata(gene)$gene_name)
#loadedGenomes[[genome]] <- list(
#    geneNames=geneNames,
#    dbGene=gene,
#    dbExon=NULL
#)

# Keep track of loaded data and load the first
loadedData <- vector("list",length(sources))
names(loadedData) <- sources
for (s in sources) {
    dd <- as.character(dbGetQuery(metadata, paste0("SELECT DISTINCT(dataset) FROM metadata WHERE source == '",s,"'"))$dataset)
    loadedData[[s]] <- vector("list",length(dd))
    names(loadedData[[s]]) <- dd
}

# Restrict the number of cores dedicated to SeqCVIBE
RC <- 0.25
#RC <- NULL

# Hack for p-value, FDR display in sliders and natural fold change
statScoreValues <- c(0,1e-16,1e-8,1e-4,0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.1,
    0.2,0.3,0.4,0.5,1)
fcNatValues <- c(0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,
    0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10)
