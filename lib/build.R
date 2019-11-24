getRpkm <- function(source,dataset,class,config,dbExon,activeLength=NULL,
    saveCountObj=NULL,output=c("counts","rpkm"),outBase="getRpkm_out_",
    outDir=".",rc=NULL) {
    targets <- readTargetsFromConfig(source,dataset,class,config)
    if (is.null(activeLength))
        activeLength <- rep(1,length(dbExon))
    annotation <- list(ranges=dbExon,length=activeLength)
    b2c.out <- bam2count(targets,annotation,rc=rc)
    if (!is.null(saveCountObj))
        save(b2c.out,file=file.path(outDir,paste(saveCountObj,"rda",sep=".")))
    gene.counts <- b2c.out$counts
    libsize <- b2c.out$libsize
    gene.length <- b2c.out$length
    gene.rpkms <- edgeR::rpkm(gene.counts,gene.length=gene.length,
        lib.size=unlist(libsize))
    
    if (!is.null(config$alt_id)) {
        ind <- which(config$source %in% source & config$dataset %in% dataset 
            & config$class %in% class)
        subconf <- config[ind,]
        alt_ids <- as.character(subconf$alt_id)
        names(alt_ids) <- as.character(subconf$sample_id)
        alt_ids <- alt_ids[colnames(gene.counts)]
        colnames(gene.counts) <- colnames(gene.rpkms) <- alt_ids
    }
    na <- which(is.na(colnames(gene.rpkms)))
    if (length(na)>0) {
        gene.counts <- gene.counts[,-na]
        gene.rpkms <- gene.rpkms[,-na]
    }
    
    if ("counts" %in% output)
        write.table(gene.counts,file=file.path(outDir,paste(outBase,
            "_counts.txt")),sep="\t",quote=FALSE,col.names=NA)
    if ("counts" %in% output)
        write.table(gene.rpkms,file=file.path(outDir,paste(outBase,
            "_rpkm.txt")),sep="\t",quote=FALSE,col.names=NA)
}

bam2count <- function(targets,annotation,rc=NULL) {
    if (missing(targets))
        stop("You must provide the targets argument!")
    if (missing(annotation))
        stop("You must provide an annotation data frame!")
    if (!requireNamespace("GenomicRanges"))
        stop("The Bioconductor package GenomicRanges is required to ",
            "proceed!")
    if (!requireNamespace("GenomicAlignments"))
        stop("The Bioconductor package GenomicAlignments is required to ",
            "proceed!")
    if (!requireNamespace("Rsamtools"))
        stop("The Bioconductor package Rsamtools is required to ",
            "process BAM files!")
    if (!is.list(targets)) {
        if (is.character(targets)) {
            if (file.exists(targets))
                targets <- readTargets(targets)
            else
                stop("You must provide a targets list or a valid targets file!")
        }
    }
    
    annotation.gr <- annotation$ranges
    coding.length <- annotation$length
    
    files.list <- targets$files
    sample.names <- unlist(lapply(files.list,names),use.names=FALSE)
    sample.files <- unlist(files.list,use.names=FALSE)
    names(sample.files) <- sample.names
    if (!is.null(targets$paired)) {
        paired <- unlist(targets$paired,use.names=FALSE)
        names(paired) <- sample.names
    }
    else
        paired <- NULL
    if (!is.null(targets$stranded)) {
        stranded <- unlist(targets$stranded,use.names=FALSE)
        names(stranded) <- sample.names
    }
    else
        stranded <- NULL
    
    counts <- matrix(0,nrow=length(annotation.gr),ncol=length(sample.names))
    rownames(counts) <- names(annotation.gr)
    colnames(counts) <- sample.names
    libsize <- vector("list",length(sample.names))
    names(libsize) <- sample.names
        
    ret.val <- cmclapply(sample.names,function(n,sample.files,paired,
        stranded) {
        message("Reading bam file ",basename(sample.files[n])," for sample ",
            "with name ",n,". This might take some time...")
        if (!is.null(paired)) {
            p <- tolower(paired[n])
            if (p=="single") {
                singleEnd <- TRUE
                fragments <- FALSE
                asMates <- FALSE
            }
            else if (p=="paired") {
                singleEnd <- FALSE
                fragments <- FALSE
                asMates <- TRUE
            }
            else if (p=="mixed") {
                singleEnd <- FALSE
                fragments <- TRUE
                asMates <- TRUE
            }
            else {
                warning("Information regarding single- or paired-end ",
                    "reads is not correctly provided! Assuming single...",
                    immediate.=TRUE)
                singleEnd <- TRUE
                fragments <- FALSE
                asMates <- FALSE
            }
        }
        else {
            singleEnd <- TRUE
            fragments <- FALSE
            asMates <- FALSE
        }
        if (!is.null(stranded)) {
            s <- tolower(stranded[n])
            if (s %in% c("forward","reverse"))
                ignore.strand <- FALSE
            else if (s=="no")
                ignore.strand <- TRUE
            else {
                warning("Information regarding strandedness of the reads ",
                    "is not correctly provided! Assuming unstranded...",
                    immediate.=TRUE)
                ignore.strand <- TRUE
            }
        }
        else
            ignore.strand <- TRUE
        # Check remoteness
        if (length(grep("^(http|ftp)",sample.files[n],perl=TRUE))>=1) {
            reads <- as(readGAlignments(file=sample.files[n]),"GRanges")
            libsize <- length(reads)
            is.remote <- TRUE
        }
        else {
            reads <- BamFile(sample.files[n],asMates=asMates)
            libsize <- countBam(reads,
                param=ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE)))$records
            is.remote <- FALSE
        }
        if (libsize>0) {
            message("  Counting reads overlapping with given annotation...")
            if (singleEnd & !fragments)
                message("    ...for single-end reads...")
            else if (!singleEnd & !fragments)
                message("    ...for paired-end reads...")
            else if (!singleEnd & fragments)
                message("    ...for mixed single- and paired-end reads...")
            if (ignore.strand)
                message("    ...ignoring strandedness...")
            else {
                message("    ...assuming ",s," sequenced reads...")
                if (s=="reverse")
                    strand(annotation.gr) <- ifelse(strand(
                        annotation.gr)=="+","-","+")
            }
            if (is.remote)
                message("    ...for remote BAM file... might take longer...")
            counts <- tryCatch(
                summarizeOverlaps(annotation.gr,reads),#singleEnd=singleEnd,
                    #fragments=fragments,ignore.strand=ignore.strand),
                error=function(e) {
                    message("Caught error while reading BAM file: ",
                    sample.files[n])
                    message("====================")
                    print(e)
                    message("====================")
                    return(NULL)
                },finally=""
            )
            if (!is.null(counts))
                counts <- assays(counts)$counts
        }
        else {
            warning(paste("No reads left after annotation chromosome ",
                "presence check for sample ",n,sep=""))
            counts <- NULL
        }
        gc(verbose=FALSE)
        return(list(counts=counts,libsize=libsize))
    },sample.files,paired,stranded,rc=rc)
    
    failed <- numeric(0)
    for (i in 1:length(ret.val)) {
        if (!is.null(ret.val[[i]]$counts)) {
            counts[,i] <- ret.val[[i]]$counts
            libsize[[i]] <- ret.val[[i]]$libsize
        }
        else
            failed <- c(failed,i)       
    }
    if (length(failed)>0) {
        counts <- counts[,-failed]
        libsize <- libsize[-failed]
    }
    return(list(counts=counts,libsize=libsize,length=coding.length))
}

readTargets <- function(input,path=NULL) {
    if (missing(input) || !file.exists(input))
        stop("The targets file should be a valid existing text file!")
    tab <- read.delim(input,strip.white=TRUE)
    samples <- as.character(tab[,1])
    conditions <- unique(as.character(tab[,3]))
    rawfiles <- as.character(tab[,2])
    if (!is.null(path)) {
        tmp <- dirname(rawfiles) # Test if there is already a path
        if (any(tmp=="."))
            rawfiles <- file.path(path,basename(rawfiles))
    }
    # Check if they exist!!!
    for (f in rawfiles) {
        if (!file.exists(f))
            stop("Raw reads input file ",f," does not exist! Please check!")
    }
    if (length(samples) != length(unique(samples)))
        stop("Sample names must be unique for each sample!")
    if (length(rawfiles) != length(unique(rawfiles)))
        stop("File names must be unique for each sample!")
    sample.list <- vector("list",length(conditions))
    names(sample.list) <- conditions
    for (n in conditions)
        sample.list[[n]] <- samples[which(as.character(tab[,3])==n)]
    file.list <- vector("list",length(conditions))
    names(file.list) <- conditions
    for (n in conditions) {
        file.list[[n]] <- rawfiles[which(as.character(tab[,3])==n)]
        names(file.list[[n]]) <- samples[which(as.character(tab[,3])==n)]
    }
    if (ncol(tab)>3) { # Has info about single- or paired-end reads / strand
        if (ncol(tab)==4) { # Stranded or paired
            whats <- tolower(as.character(tab[,4]))
            if (!all(whats %in% c("yes","no","forward","reverse",
                "single","paired")))
                stop("Unknown options for paired-end reads and/or ",
                    "strandedness in targets file.")
            what <- whats[1]
            if (what %in% c("single","paired")) {
                has.paired.info <- TRUE
                has.stranded.info <- FALSE
            }
            else {
                if (what %in% c("yes","no")) {
                    deprecated.warning("read.targets")
                    tmp <- as.character(tab[,4])
                    tmp[tmp=="yes"] <- "forward"
                    tab[,4] <- tmp
                    has.paired.info <- FALSE
                    has.stranded.info <- TRUE
                }
                if (what %in% c("forward","reverse","no")) {
                    has.paired.info <- FALSE
                    has.stranded.info <- TRUE
                }
            }
        }
        if (ncol(tab)==5) { # Both
            whats.paired <- tolower(as.character(tab[,4]))
            if (!all(whats.paired %in% c("single","paired","mixed")))
                stop("Unknown option for type of reads (single, paired, ",
                    "mixed) in targets file.")
            whats.strand <- tolower(as.character(tab[,5]))
            if (!all(whats.strand %in% c("yes","no","forward","reverse")))
                stop("Unknown option for read strandedness in targets file")
            has.paired.info <- TRUE
            has.stranded.info <- TRUE
        }
        if (has.paired.info && !has.stranded.info) {
            paired.list <- vector("list",length(conditions))
            names(paired.list) <- conditions
            for (n in conditions) {
                paired.list[[n]] <- character(length(sample.list[[n]]))
                names(paired.list[[n]]) <- sample.list[[n]]
                for (nn in names(paired.list[[n]]))
                    paired.list[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
        else
            paired.list <- NULL
        if (has.stranded.info && !has.paired.info) {
            stranded.list <- vector("list",length(conditions))
            names(stranded.list) <- conditions
            for (n in conditions) {
                stranded.list[[n]] <- character(length(sample.list[[n]]))
                names(stranded.list[[n]]) <- sample.list[[n]]
                for (nn in names(stranded.list[[n]]))
                    stranded.list[[n]][nn] <- as.character(tab[which(
                        as.character(tab[,1])==nn),4])
            }
        }
        else
            stranded.list <- NULL
        if (has.stranded.info && has.paired.info) {
            stranded.list <- vector("list",length(conditions))
            names(stranded.list) <- conditions
            for (n in conditions) {
                stranded.list[[n]] <- character(length(sample.list[[n]]))
                names(stranded.list[[n]]) <- sample.list[[n]]
                for (nn in names(stranded.list[[n]]))
                    stranded.list[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),5])
            }
            paired.list <- vector("list",length(conditions))
            names(paired.list) <- conditions
            for (n in conditions) {
                paired.list[[n]] <- character(length(sample.list[[n]]))
                names(paired.list[[n]]) <- sample.list[[n]]
                for (nn in names(paired.list[[n]]))
                    paired.list[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
    }
    else
        paired.list <- stranded.list <- NULL
    # Guess file type based on only one of them
    tmp <- file.list[[1]][1]
    if (length(grep("\\.bam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "bam"
    else if (length(grep("\\.sam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "sam"
    else if (length(grep("\\.bed$",tmp,ignore.case=TRUE,perl=TRUE)>0))
        type <- "bed"
    else
        type <- NULL
    return(list(samples=sample.list,files=file.list,paired=paired.list,
        stranded=stranded.list,type=type))
}

readTargetsFromConfig <- function(source,dataset,class,config) {
    ind <- which(config$source %in% source & config$dataset %in% dataset 
        & config$class %in% class)
    subconf <- config[ind,]
    samples <- as.character(subconf$sample_id)
    conditions <- unique(as.character(subconf$class))
    rawfiles <- as.character(subconf$sample_dir)
    for (f in rawfiles) {
        if (!dir.exists(f))
            stop("Raw reads input dir ",f," does not exist! Please check!")
    }
    if (length(samples) != length(unique(samples)))
        stop("Sample names must be unique for each sample!")
    if (length(rawfiles) != length(unique(rawfiles)))
        stop("File directories must be unique for each sample!")
    sample.list <- vector("list",length(conditions))
    names(sample.list) <- conditions
    for (n in conditions)
        sample.list[[n]] <- 
            as.character(subconf$sample_id[which(subconf$class==n)])
    file.list <- vector("list",length(conditions))
    names(file.list) <- conditions
    for (n in conditions) {
        dirs <- as.character(subconf$sample_dir[which(subconf$class==n)])
        ss <- as.character(subconf$sample_id[which(subconf$class==n)])
        file.list[[n]] <- dir(dirs,pattern=".bam$",full.names=TRUE)
        names(file.list[[n]]) <- ss
    }
    
    paired.list <- stranded.list <- NULL
    tmp <- file.list[[1]][1]
    if (length(grep("\\.bam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "bam"
    else if (length(grep("\\.sam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "sam"
    else if (length(grep("\\.bed$",tmp,ignore.case=TRUE,perl=TRUE)>0))
        type <- "bed"
    else
        type <- NULL
    return(list(samples=sample.list,files=file.list,paired=paired.list,
        stranded=stranded.list,type=type))
}

# annoPath : The path where JBrowse formatted reference tracks (genomic 
#            sequences and features) live
# urlBase  : The URL base (or template) to be prepended to the actual user data,
#            (not the reference, these are constructed by JBrowse control 
#            wrapper scripts) tracks. It should point to a directory served by 
#            a web-browser. This directory should be the one with the data, or 
#            with proper symlinks pointing to the data
# appBase  : The application base, default "../"
# metadata : The path to the application metadata SQLite database
buildTrackList  <- function(annoPath,urlBase,appBase="../",
    metadata="../config/metadata.sqlite") {
    #annoPath <- "/media/raid/tracks/seqcvibe/reference"
    #appBase <- "/media/raid/software/seqcvibe"
    #urlBase <- "http://epigenomics.fleming.gr/tracks"
    
    if (!requireNamespace("jsonlite"))
        stop("R package jsonlite is required to proceed!")
    if (!requireNamespace("RSQLite"))
        stop("R package RSQLite is required to proceed!")
    if (!requireNamespace("plyr"))
        stop("R package plyr is required to proceed!")
        
    if (!file.exists(metadata))
        stop("The metadata database does not exist!")
    
    # annoPath contains:
    #   names
    #   seq
    #   annotation feature tracks
    #   annotation feature tracks tracklist
    # For the above we must:
    #   create a "tracks/hg19" folder in app base
    #   create symlinks to names, seq, trakcs
    #   copy annotation feature tracks tracklist
    #   create our track list
    #   merge them
    # All these for each source and dataset
    #   e.g.: TCGA/COAD/seq, TCGA/COAD/trackList etc.
    
    # Track colors - to be recycled in case of many classes for each dataset
    posBaseColours <- c("#B40000","#00B400","#0000B4","#B45200","#9B59B6",
        "#21BCBF","#BC4800","#135C34","#838F00","#4900B5")
    negBaseColours <- c("#FF7575","#79FF79","#8484FF","#FFB77C","#EBB7FF",
        "#63E5E7","#FFB88B","#5BFFA5","#F4FF75","#C69FFF")
    
    # Open the metadata database
    con <- dbConnect(dbDriver("SQLite"),dbname=metadata)
    
    allData <- dbGetQuery(con,"SELECT * FROM metadata")
    orgs <- as.character(
        dbGetQuery(con,"SELECT DISTINCT genome FROM metadata")[,1])
    classes <- as.character(
        dbGetQuery(con,"SELECT DISTINCT class FROM metadata")[,1])
    
    #orgs <- unique(as.character(config$genome))
    #classes <- unique(as.character(config$class))
    
    # Track colors
    posBaseColours <- rep(posBaseColours,length.out=length(classes))
    negBaseColours <- rep(negBaseColours,length.out=length(classes))
    names(posBaseColours) <- names(negBaseColours) <- classes

    # Create tracks base and symlinks to annotation sources
    for (org in orgs) {
        if (!dir.exists(file.path(appBase,"tracks","reference",org)))
            dir.create(file.path(appBase,"tracks","reference",org),
                recursive=TRUE)
        if (!file.exists(file.path(appBase,"tracks","reference",org,"names")))
            file.symlink(file.path(annoPath,org,"names"),
                file.path(appBase,"tracks","reference",org,"names"))
        if (!file.exists(file.path(appBase,"reference","tracks",org,"seq")))
            file.symlink(file.path(annoPath,org,"seq"),
                file.path(appBase,"tracks","reference",org,"seq"))
        if (!file.exists(file.path(appBase,"tracks","reference",org,"tracks")))
            file.symlink(file.path(annoPath,org,"tracks"),
                file.path(appBase,"tracks","reference",org,"tracks"))
        
        # Read in JSON from annotation trackList
        annoTracks <- fromJSON(file.path(annoPath,org,"trackList.json"),
            simplifyVector=FALSE)
            
        # Construct the rest of the tracks to merge with annoTracks...
        tracksXY <- dlply(allData[allData$genome==org,,drop=FALSE],"sample_id",
            function(x,poscol,negcol) {
            list(
                metadata=list(
                    Source=as.character(x$source),
                    Dataset=as.character(x$dataset),
                    Class=as.character(x$class),
                    Organism=getFriendlyName(org),
                    Type="Signal",
                    User=NULL
                ),
                style=list(
                    height=64,
                    pos_color=poscol[as.character(x$class)],
                    neg_color=negcol[as.character(x$class)]
                ),
                autoscale="local",
                key=ifelse(is.null(x$alt_id),as.character(x$sample_id),
                    as.character(x$alt_id)),
                type="JBrowse/View/Track/Wiggle/XYPlot",
                storeClass="JBrowse/Store/SeqFeature/BigWig",
                yScalePosition="right",
                label=ifelse(is.null(x$alt_id),paste(x$sample_id,"xy",sep="_"),
                    paste(x$sample_id,x$alt_id,"xy",sep="_")),
                #urlTemplate=paste(urlBase,"/",x$source,"/",x$dataset,"/",
                #    x$class,"/",x$sample_id,".bigWig",sep="")
                urlTemplate=paste(urlBase,"/",x$dataset,"/",x$class,"/",
                    x$sample_id,".bigWig",sep="")
            )
        },posBaseColours,negBaseColours)
        tracksDen <- dlply(allData[allData$genome==org,,drop=FALSE],"sample_id",
            function(x,poscol,negcol) {
            list(
                metadata=list(
                    Source=as.character(x$source),
                    Dataset=as.character(x$dataset),
                    Class=as.character(x$class),
                    Organism=getFriendlyName(org),
                    Type="Density",
                    User=NULL
                ),
                style=list(
                    height=16,
                    pos_color=poscol[as.character(x$class)],
                    neg_color=negcol[as.character(x$class)]
                ),
                autoscale="local",
                key=ifelse(is.null(x$alt_id),as.character(x$sample_id),
                    as.character(x$alt_id)),
                type="JBrowse/View/Track/Wiggle/Density",
                storeClass="JBrowse/Store/SeqFeature/BigWig",
                yScalePosition="right",
                label=ifelse(is.null(x$alt_id),paste(x$sample_id,"density",
                    sep="_"),paste(x$sample_id,x$alt_id,"density",sep="_")),
                #urlTemplate=paste(urlBase,"/",x$source,"/",x$dataset,"/",
                #    x$class,"/",x$sample_id,".bigWig",sep="")
                urlTemplate=paste(urlBase,"/",x$dataset,"/",x$class,"/",
                    x$sample_id,".bigWig",sep="")
            )
        },posBaseColours,negBaseColours)
        tracks <- unname(c(tracksXY,tracksDen))
        
        annoTracks$tracks <- c(annoTracks$tracks,tracks)
        output <- toJSON(annoTracks,auto_unbox=TRUE,pretty=TRUE)
        write(output,file=file.path(appBase,"tracks","reference",org,
            "trackList.json"))
    }
}

modifyTrackList  <- function(trackList,what=list(annoBase=NULL,urlBase=NULL),
    appBase="") {
    # Read in existing trackList
    # If annoBase is not NULL, change the annotation base
    # If urlBase is not NULL, change the "urlTemplate" for user tracks
}

addTracks  <- function(trackList,metadata,trackMetadata=list(),appBase="") {
    # Read in existing trackList
    # Open metadata and choose tracks to add from trackMetadata
}

clearUserTracks <- function(trackList,user=NULL) {
    trackStru <- fromJSON(trackList,simplifyVector=FALSE)
    
    # Initially we clear based only on storeClass... More complex ways will be
    # required. For now, "JBrowse/Store/SeqFeature/BigWig"
    # We include the "user" argument for future use
    tracks <- trackStru$tracks
    inds <- sapply(1:length(tracks),function(i,tr) {
        if (tr[[i]]$storeClass == "JBrowse/Store/SeqFeature/BigWig")
            return(i)
        else
            return(0)
    },tracks)
    if (length(inds) > 0) {
        nouser <- which(inds==0)
        if (length(nouser) > 0)
            inds <- inds[-nouser]
    }
    if (length(inds) > 0)
        tracks <- tracks[-inds]
    
    trackStru$tracks <- tracks
    cleared <- toJSON(trackStru,auto_unbox=TRUE,pretty=TRUE)
    write(cleared,file=trackList)
}

## Build track list example
#config <- read.delim("/media/raid/software/bigseqcvis/config/metadata.txt")
#annoPath <- "/media/raid/tracks/hybridsuite/reference"
#appBase <- "/media/raid/software/bigseqcvis"
#urlBase <- "http://epigenomics.fleming.gr/tracks"
#buildTrackList(
#    config=config,
#    annoPath=annoPath,
#    urlBase=urlBase,
#    appBase=appBase
#)

#chr11:2,218,988-2,299,254
#
#refArea <- GRanges(
#   seqnames="chr11",
#   IRanges(start=2218988,end=2299254)
#)
#
#wntr <- data.frame(
#    chromosome="chr11",
#    start=2228236,
#    end=2232202,
#    strand=Rle("-"),
#    exon_id="WNTRLINC1",
#    gene_id="WNTRLINC1",
#    gene_name="WNTRLINC1",
#    biotype="lincRNA"
#)
#wntr.gr <- makeGRangesFromDataFrame(
#    df=wntr,
#    keep.extra.columns=TRUE
#)
#
#dbGene[subjectHits(findOverlaps(refArea,dbGene)]
#
#2219236
#2232100

## Build track list example
#config <- read.delim("/media/raid/software/seqcvibe/config/metadata.txt")
#annoPath <- "/media/raid/tracks/seqcvibe/reference"
#appBase <- "/media/raid/software/seqcvibe"
#urlBase <- "https://epigenomics.fleming.gr/tracks"
#buildTrackList(
#    config=config,
#    annoPath=annoPath,
#    urlBase=urlBase,
#    appBase=appBase
#)
#perl buildReferenceTracks.pl --config ../config/config_jbrowse_build.json --genomes --features --update --ncores 5


#annoPath <- "/media/raid/tracks/seqcvibe/reference"
#appBase <- "/media/data/hybridstat/elixir-RNAseq"
#metadata <- file.path(appBase,"config","metadata.sqlite")
#urlBase <- "http://epigenomics.fleming.gr/seqc_elixir"

