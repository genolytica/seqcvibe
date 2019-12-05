# Intit globals (don't want to laod these everytime a use connects and they 
# don't take much time to load)

# Load basic packages
library(shiny)

library(jsonlite)
library(RSQLite)
library(shinyjs)

# colourpicker must be loaded AFTER shinyjs... Sigh...
library(colourpicker)

# Load additional functions
source("lib/control.R")
source("lib/db.R")
source("lib/util.R")
# source("lib/queries.R")
# QUERIES <- initQueries()

appConfig <- fromJSON("config/app_config.json")

# If the symlink to the data directory is not in place, create it!
if (!isTRUE(nzchar(Sys.readlink("tracks/data"),keepNA=TRUE)))
    file.symlink(appConfig$paths$data,"tracks/data")

# Load metadata SQLite
metadata <- initDatabase(file.path(appConfig$paths$metadata))

## Until we merge bookmarks with main db
# Load Bookmarks SQLite
#bookmarks <- dbConnect(drv=RSQLite::SQLite(),dbname='data/bookmarks.sqlite')

# Intialize metadata reactive content
sources <- as.character(dbGetQuery(metadata,
    "SELECT DISTINCT source FROM metadata")$source)
datasets <- as.character(dbGetQuery(metadata,paste0("SELECT DISTINCT dataset ",
    "FROM metadata WHERE source='",sources[1],"'"))$dataset)
classes <- as.character(dbGetQuery(metadata, paste0("SELECT DISTINCT class ",
    "FROM metadata WHERE source='",sources[1],"' AND dataset=='",datasets[1],
    "'"))$class)
genomes <- as.character(dbGetQuery(metadata,
    "SELECT DISTINCT genome FROM metadata")$genome)
genome <- genomes[1]
title <- as.character(dbGetQuery(metadata,paste0("SELECT DISTINCT title FROM ",
    "datasets WHERE dataset='",datasets[1],"'"))$title)
link <- as.character(dbGetQuery(metadata,paste0("SELECT DISTINCT link  FROM ",
    "datasets WHERE dataset='",datasets[1],"'"))$link)
short_summary <- as.character(dbGetQuery(metadata,paste0("SELECT  DISTINCT ",
    "short_summary FROM datasets WHERE dataset='",datasets[1],
    "'"))$short_summary)

## Load data file hash
#load("data/dataFiles.RData")

# Create the data files list
tmpAll <- dbGetQuery(metadata,"SELECT DISTINCT source,dataset FROM metadata")
us <- unique(tmpAll$source)
dataFiles <- vector("list",length(us))
names(dataFiles) <- us
for (s in us) {
    tmpD <- tmpAll[tmpAll$source==s,,drop=FALSE]
    dataFiles[[s]] <- vector("list",nrow(tmpD))
    names(dataFiles[[s]]) <- tmpD$dataset
    for (d in tmpD$dataset)
        dataFiles[[s]][[d]] <- 
            file.path(appConfig$paths$data,d,paste0(d,".rda"))
}

# source("config/data_files.R")
allClasses <- as.character(dbGetQuery(metadata,
    "SELECT DISTINCT class FROM metadata")$class)

baseColours <- c("#B40000","#00B400","#0000B4","#B45200","#9B59B6","#21BCBF",
    "#BC4800","#135C34","#838F00","#4900B5")
baseColours <- rep(baseColours,length.out=length(allClasses))
names(baseColours) <- allClasses

# Keep track of loaded annotations. The first will be loaded when needed.
loadedGenomes <- vector("list",length(genomes))
names(loadedGenomes) <- genomes
for (gen in genomes) {
    loadedGenomes[[gen]] <- list(
        geneNames=NULL,
        dbGene=NULL,
        dbExon=NULL
    )
}

# Keep track of loaded data and load the first
loadedData <- vector("list",length(sources))
names(loadedData) <- sources
for (s in sources) {
    dd <- as.character(dbGetQuery(metadata,paste0("SELECT DISTINCT dataset ",
        "FROM metadata WHERE source='",s,"'"))$dataset)
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
