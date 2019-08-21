library(readr)
library(dplyr)
library(RSQLite)

files_dir <- "/home/makis/elixir-RNAseq/"
files <- list.files(recursive = T, pattern = "*.bam$")

core <- t(as.data.frame(strsplit(files, split = "/")))
colnames(core) <- c("dataset", "class", "sample_id")
rownames(core)<-NULL
metadata<-as.data.frame(core)

metadata$sample_id <- gsub(pattern = ".bam", replacement = "", x = metadata$sample_id)
metadata$source[c(grep(pattern = "PRJEB", x = metadata$dataset))] <- "ENA"
metadata$source[c(grep(pattern = "GSE", x = metadata$dataset))] <- "GEO"
metadata$sample_dir <- paste0(files_dir, "data/", metadata$dataset, "/", metadata$class)
metadata$track_dir <- paste0(files_dir, "data/", metadata$dataset, "/", metadata$class)
metadata$alt_id <- paste0(metadata$sample_id, "_", metadata$class)
metadata$library_strategy <- "RNASeq"
metadata$quality <- 0


# FUNCTION 1 apply to all norm_factor files and merge into 1 DF
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


conn <- dbConnect(RSQLite::SQLite(), "/home/makis/elixir-RNAseq/data/metadata.sqlite")
dbWriteTable(conn, "metadata", final_metadata, overwrite = T)

#dbGetQuery(conn, "SELECT * FROM metadata WHERE genome = 'mm10'")

