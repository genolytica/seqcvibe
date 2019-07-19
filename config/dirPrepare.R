library(readr)
library(filesstrings)
dataset <- list.files(pattern = "PRJEB|GSE")

dirSort <- function(dataset) {
  targets <- read_delim(paste0(dataset,"/metaseqR_out/targets.txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
  conditions <- unique(targets$condition)
  
  file.move(files = paste0(dataset,"/metaseqR_out/",dataset,".rda"), destinations = dataset)
  file.move(files = paste0(dataset,"/bigwig/",dataset,"_factors.txt"), destinations = dataset)
  file.move(files = paste0(dataset,"/genome"), destinations = dataset)
  
  for (con in conditions) {
    dir.create(file.path(dataset,con))
    samples <- targets$samplename[which(targets$condition==con)]
    for (sample in samples) {
      file.move(files = paste0(dataset,"/hisat2_out/",sample,paste0("/",sample,".bam")),destinations = file.path(dataset,con))
      file.move(files = paste0(dataset,"/hisat2_out/",sample,paste0("/",sample,".bam.bai")),destinations = file.path(dataset,con))
      file.move(files = paste0(dataset,"/multiqc/",paste0(sample,"_multiqc_report.html")),destinations = file.path(dataset,con))
      file.move(files = paste0(dataset,"/bigwig/",paste0(sample,".bigWig")),destinations = file.path(dataset,con))
    }
  }
  unlink(paste0(dataset,"/fastq/"), recursive = T)
  unlink(paste0(dataset,"/multiqc/"), recursive = T)
  unlink(paste0(dataset,"/hisat2_out/"), recursive = T)
  unlink(paste0(dataset,"/metaseqR_out/"), recursive = T)
  unlink(paste0(dataset,"/bigwig/"), recursive = T)
  unlink(paste0(dataset,"/design.txt"))
}

lapply(dataset, dirSort)

