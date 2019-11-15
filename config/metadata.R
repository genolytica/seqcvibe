# Run in "/data" directory 
# Rscript metadata.R <install_dir>

library(readr)
library(dplyr)
library(RSQLite)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("You need to provide a full SeqCVIBE install directory as an arguement", call.=FALSE)
}

# Create metadata table
install_dir <- args[1]
files <- list.files(path=paste0(install_dir,"data/"), recursive = T, pattern = "*.bam$")

core <- t(as.data.frame(strsplit(files, split = "/")))
colnames(core) <- c("dataset", "class", "sample_id")
rownames(core)<-NULL
metadata<-as.data.frame(core)

metadata$sample_id <- gsub(pattern = ".bam", replacement = "", x = metadata$sample_id)
metadata$source[c(grep(pattern = "PRJEB", x = metadata$dataset))] <- "ENA"
metadata$source[c(grep(pattern = "GSE", x = metadata$dataset))] <- "GEO"
metadata$sample_dir <- paste0(install_dir, "data/", metadata$dataset, "/", metadata$class)
metadata$track_dir <- paste0(install_dir, "data/", metadata$dataset, "/", metadata$class)
metadata$alt_id <- paste0(metadata$sample_id, "_", metadata$class)
metadata$library_strategy <- "RNASeq"
metadata$quality <- 0

IDs <- as.list(list.files(recursive = F, pattern = "PRJEB|GSE"))

fillNormFactor <- function(ID) {
  genome <- read.table(paste0(ID,"/genome"), quote="\"", comment.char="")
  ID <- read_table2(paste0(ID, "/", ID, "_factors.txt"))
  ID<-ID[,c(1:2)]
  ID$genome <- as.character(genome[[1]])
  colnames(ID) <- c("sample_id", "norm_factor", "genome")
  ID$sample_id <- gsub(pattern = ".bedGraph", replacement = "", ID$sample_id)
  metadata <- merge(x = metadata, y = ID, by = 'sample_id')
}
final_metadata <- lapply(IDs, fillNormFactor)%>% bind_rows()

# Create summaries table
summaries <- read.delim("summaries.tsv", header=FALSE, stringsAsFactors=FALSE)
colnames(summaries) <- c("dataset", "title", "link", "short_summary")

# Write tables to metadata.sqlite
conn <- dbConnect(RSQLite::SQLite(), paste0(install_dir,"data/metadata.sqlite"))
dbWriteTable(conn, "metadata", final_metadata, overwrite = T)
dbWriteTable(conn, "summaries", summaries, overwrite = T)
dbDisconnect(conn)

# Create the dataFiles hash-type list
dataFiles <- list(
    GEO=list(),
    ENA=list()
    )
namesGEO<-c()
namesENA<-c()

for (dat in unique(final_metadata$dataset)) {
	if (unique(final_metadata$source[which(final_metadata$dataset==dat)])=="GEO") {
		dataFiles$GEO <- append(dataFiles$GEO, paste0("data/",dat,"/",dat,".rda"))
		namesGEO <- c(namesGEO,dat)
	}
	else if (unique(final_metadata$source[which(final_metadata$dataset==dat)])=="ENA") {
		dataFiles$ENA <- append(dataFiles$ENA, paste0("data/",dat,"/",dat,".rda"))
		namesENA <- c(namesENA,dat)
	}
}
names(dataFiles$GEO) <- namesGEO
names(dataFiles$ENA) <- namesENA

save(dataFiles, file=paste0(install_dir,"/data/dataFiles.RData"))

# TEST
#dbGetQuery(conn, "SELECT * FROM metadata WHERE genome = 'mm10'")
