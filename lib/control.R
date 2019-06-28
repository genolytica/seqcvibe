getRnaGene <- function(
    gene, # A list of mixed __named__ objects
    refArea=NULL, # Flank or chr-start-end
    tumor,
    class=NULL, # If not provided, both tumor and normal
    dbGene,
    dbExon=NULL,
    config,
    exclude=list(normal=NULL,tumor=NULL),
    sumStat=c("mean","median","trimmed"),
    trim=0.1,
    plot=c("spline","track"),
    rc=NULL
) {
    # Requires kind of complex validation
    validateGeneInput(gene) # More complex validation
    
    if (!is(dbGene,"GRanges")) {
        # Might be an rda file
        tryCatch({
            load(dbGene)
            dbGene <- gene
        },error=function(e) {
            stop("Caught unexpected error while loading annotation: ",e)
        },finally="")
    }
    if (is.null(dbExon) && plot=="track")
        stop("The dbExon argument must be provided for track plots!")
    if (plot=="track") {
        if (!is(dbExon,"GRangesList")) {
            # Might be an rda file
            tryCatch({
                load(dbExon)
                dbExon <- exon
            },error=function(e) {
                stop("Caught unexpected error while loading annotation: ",e)
            },finally="")
        }
    }
    if (!is.list(config)) { # Might be an rda file
        tryCatch({
            load(config)
            config <- config
        },error=function(e) {
            stop("Caught unexpected error while loading configuration: ",e)
        },finally="")
    }
    if (!(tumor %in% names(config)))
        stop("The requested tumor code does is not yet supported!")
    class <- tolower(class)
    if (is.null(class))
        class <- c("tumor","normal")
    if (!(class %in% c("tumor","normal")))
        stop("class must be one or more of \"tumor\" and \"normal\"!")
    if (!is.null(exclude) && !(names(exclude) %in% c("tumor","normal")))
        stop("The exclude argument is not valid! Check its names.")
    if (!is.null(refArea)) {
        if (plot=="spline") {
            if (!is.numeric(refArea) || length(refArea)!=2)
                stop("The refArea argument must be a numeric vector of length ",
                    "2 denoting the flanking areas of the request genes.")
            if (any(refArea<500 || refArea>5000))
                stop("The minimum flanking area allowed is 500bp and the ",
                    "maximum is 5000bp")
        }
        if (plot=="track") {
            if (!is.list(refArea))
                stop("The refArea argument must be a list with three members, ",
                    "chr, start, end.")
            if (!all(names(refArea) %in% c("chr","start","end")))
                stop("The refArea argument must be named with the names chr, ",
                    "start, end")
            if (!is.character(refArea$chr) || !is.numeric(refArea$start)
                || !is.numeric(refArea$end))
                stop("Invalid refArea member types: must be character, numeric",
                    "numeric")
        }
    }
    
    # Define BAM paths to get coverage and check if samples are excluded
    for (status in class) {
        i <- match(exclude[[status]],names(config[[tumor]]$samples[[status]]))
        if (length(i)>0)
            config[[tumor]]$samples[[status]] <- 
                config[[tumor]]$samples[[status]][-i]
    }
    
    # Create plot objects
    if (plot=="spline") {
        ggObj <- getProfile(gene,refArea,class,sumStat,config[[tumor]],dbGene,
            trim,rc)
    }
    else if (plot=="track") {
        ggObj <- getTrack(gen,refArea,class,sumStat,config[[tumor]],dbGene,
            dbExon,trim=trim,rc=rc)
    }
}

getProfile <- function(gene,flank,source,dataset,class,sumStat,config,
    dbGene,trim=0.1,fromBam=FALSE,messageContainer=NULL,progressFun=NULL,
    rc=NULL) {
    # Get coordinates
    coords <- getGeneCoordinatesForSpline(gene,dbGene,flank)
    
    total <- length(coords)*length(class)
    
    # Calculate coverage over coordinates
    theCoverage <- vector("list",length(coords))
    names(theCoverage) <- names(coords)
    if (fromBam) {
        for (n in names(coords)) {
            theCoverage[[n]] <- vector("list",length(class))
            names(theCoverage[[n]]) <- class
            bsvMessage("Calculating coverage for ",n,
                messageContainer=messageContainer)
            for (st in class) {
                bsvMessage("  Processing ",st)
                bsvMessage("Calculating signal for ",n," in ",st,
                    messageContainer=messageContainer)
                if (is.function(progressFun)) {
                    text <- paste("Calculating signal for ",n," in ",st)
                    progressFun(detail=text)
                }
                ind <- which(config$source==source & config$dataset==dataset 
                    & config$class==st)
                subconf <- as.character(config$sample_dir[ind])
                names(subconf) <- as.character(config$sample_id[ind])
                theCoverage[[n]][[st]] <- cmclapply(subconf,function(x,coords) {
                    bam.file <- dir(x,pattern=".bam$",full.names=TRUE)
                    bam.index <- dir(x,pattern=".bai$",full.names=TRUE)
                    bp <- ScanBamParam(which=coords)                
                    reads <- unlist(grglist(readGAlignments(file=bam.file,
                        index=bam.index,param=bp,with.which_label=TRUE)))
                    seqlevels(reads) <- as.character(seqnames(coords))
                    return(calcCoverage(reads,coords,assign.names=FALSE,
                        verbose=FALSE)[[1]])
                },coords[n],rc=rc)
            }
        }
        
        # Normalize
        fac.index <- which(config$source==source & config$dataset==dataset 
            & config$class %in% class)
        factors <- config$norm_factor[fac.index]
        names(factors) <- as.character(config$sample_id[fac.index])
        if (is.function(progressFun)) {
            text <- paste("Normalizing...")
            progressFun(detail=text)
        }
        for (n in names(coords)) {
            for (status in names(theCoverage[[n]])) {
                bsvMessage("Normalizing ",n," in ",status)
                bsvMessage("Normalizing ",n," in ",status,
                    messageContainer=messageContainer)
                for (s in names(theCoverage[[n]][[status]])) {
                    theCoverage[[n]][[status]][[s]] <- theCoverage[[n]][[
                        status]][[s]]*factors[s]
                }
            }
        }
    }
    else {
        for (n in names(coords)) {
            theCoverage[[n]] <- vector("list",length(class))
            names(theCoverage[[n]]) <- class
            bsvMessage("Calculating coverage for ",n)
            for (st in class) {
                bsvMessage("  Processing ",st)
                bsvMessage("Calculating signal for ",n," in ",st,
                    messageContainer=messageContainer)
                if (is.function(progressFun)) {
                    text <- paste("Calculating signal for ",n," in ",st)
                    progressFun(detail=text)
                }
                
                ind <- which(config$source==source & config$dataset==dataset 
                    & config$class==st)
                sc <- as.character(config$track_dir[ind])
                names(sc) <- as.character(config$sample_id[ind])
                theCoverage[[n]][[st]] <- cmclapply(names(sc),
                    function(x,crds,sc) {
                        chr <- as.character(seqnames(crds))[1]
                        bw.file <- file.path(sc[x],paste(x,"bigWig",sep="."))
                        bwrle <- import.bw(bw.file,
                            selection=BigWigSelection(crds),as="RleList")
                        if (chr %in% names(bwrle))
                            return(bwrle[[chr]][start(crds):end(crds)])
                        else
                            return(NULL)
                },coords[n],sc,rc=rc)
            }
        }
    }

    # Average
    theAverage <- vector("list",length(coords))
    names(theAverage) <- names(coords)
    for (n in names(coords)) {
        theAverage[[n]] <- vector("list",length(class))
        names(theAverage[[n]]) <- class
    }
    for (n in names(theCoverage)) {
        for (status in (names(theCoverage[[n]]))) {
            bsvMessage("Averaging ",n," ",status," coverage",
                messageContainer=messageContainer)
            bsvMessage("Averaging ",n," ",status," coverage")
            if (is.function(progressFun)) {
                text <- paste("Averaging ",n," ",status," signal")
                progressFun(detail=text)
            }
            tmp <- cmclapply(theCoverage[[n]][[status]],function(x) {
                if (class(x)=="Rle")
                    return(as.numeric(x))
            },rc=rc)
            tmp <- do.call("rbind",tmp)
            if (sumStat=="trimmed")
                theAverage[[n]][[status]] <- apply(tmp,2,mean,trim=trim)
            else
                theAverage[[n]][[status]] <- apply(tmp,2,sumStat)
        }
    }
    
    # Smooth
    bsvMessage("Smoothing",messageContainer=messageContainer)
    bsvMessage("Smoothing")
    if (is.function(progressFun)) {
        text <- paste("Smoothing")
        progressFun(detail=text)
    }
    fit <- ci <- low <- high <- vector("list",length(coords))
    names(fit) <- names(ci) <- names(low) <- names(high) <- names(coords)
    for (n in names(coords)) {
        fit[[n]] <- ci[[n]] <- low[[n]] <- high[[n]] <- 
            vector("list",length(class))
        names(fit[[n]]) <- names(ci[[n]]) <- names(low[[n]]) <- 
            names(high[[n]]) <- class
    }
    for (n in names(theAverage)) {
        for (status in (names(theAverage[[n]]))) {
            fit[[n]][[status]] <- smooth.spline(theAverage[[n]][[status]],
                spar=0.5)
            ci[[n]][[status]] <- ssCI(fit[[n]][[status]])
            low[[n]][[status]] <- ci[[n]][[status]]$lower
            high[[n]][[status]] <- ci[[n]][[status]]$upper
            fit[[n]][[status]] <- fit[[n]][[status]]$y
        }
    }
    
    # Construct ggplot data frame
    y <- unlist(fit,use.names=FALSE)
    ymin <- unlist(low,use.names=FALSE)
    ymax <- unlist(high,use.names=FALSE)
    
    position=unlist(lapply(coords,function(x,h) {
        return(rep(start(x):end(x),h))
    },length(class)))
    ss <- locus <- vector("list",length(coords))
    names(ss) <- names(locus) <- names(fit)
    for (n in names(fit)) {
        ss[[n]] <- rep(class,each=length(fit[[n]][[1]]))
        locus[[n]] <- rep(n,length(class)*length(fit[[n]][[1]]))
    }
    ss <- unlist(ss,use.names=FALSE)
    locus <- unlist(locus,use.names=FALSE)
    
    bsvMessage("Creating plot structure",messageContainer=messageContainer)
    bsvMessage("Creating plot structure")
    if (is.function(progressFun)) {
        text <- paste("Creating plot structure")
        progressFun(detail=text)
    }
    ggplot.data <- data.frame(
        Position=position,
        Signal=y,
        Status=ss,
        Locus=locus,
        ymin=ymin,
        ymax=ymax
    )
    
    # Construct ggplot plot
    ggplot.plot <-
        ggplot(ggplot.data,mapping=aes(x=Position,y=Signal,colour=Status)) + 
        geom_line(size=0.8) +
        geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax,colour=Status,
            fill=Status),alpha=0.3,size=0) +
        facet_wrap(~ Locus,scales="free_x") +
        theme_bw() +
        xlab("\nPosition in bp") +
        ylab("Normalized signal\n") +
        theme(title=element_text(size=12),
            axis.title.x=element_text(size=10,face="bold"),
            axis.title.y=element_text(size=10,face="bold"),
            axis.text.x=element_text(size=9,face="bold"),
            axis.text.y=element_text(size=10,face="bold"),
            strip.text.x=element_text(size=10,face="bold"),
            strip.text.y=element_text(size=10,face="bold"),
            legend.position="bottom")
    
    # Return the gg object
    return(ggplot.plot)
}

getTrack <- function(refArea,customGene=NULL,source,dataset,class,sumStat,
    config,dbGene,dbExon,trim=0.1,fromBam=FALSE,classColours=NULL,
    messageContainer=NULL,progressFun=NULL,rc=NULL) {
    # Check colours
    if (is.null(classColours) || length(classColours)<length(class)) {
        baseColours <- c("#B40000","#00B400","#0000B4","#B45200")
        baseColours <- rep(baseColours,length.out=length(class))
        names(baseColours) <- class
    }
    else {
        classColours <- classColours[1:length(class)]
        names(classColours) <- class
    }
    
    # Get reference background area
    area <- getAreaCoordinatesForTrack(refArea$chr,refArea$start,refArea$end,
        dbGene)
    
    # Find if any genes (completely) overlap the area
    refArea <- GRanges(
        seqnames=refArea$chr,
        IRanges(start=refArea$start,end=refArea$end)
    )
    gene <- dbGene[subjectHits(findOverlaps(refArea,dbGene))]
    if (!is.null(customGene)) {
        for (i in 1:length(customGene)) {
            hh <- subjectHits(findOverlaps(refArea,customGene[[i]]))
            if (length(hh)>0) {
                cg <- customGene[[i]]
                cg <- validGRangesFromGRanges(cg)
                gene <- c(gene,cg)
                names(gene)[length(gene)] <- names(cg)
            }
        }
    }
    
    if (length(gene)>0) {
        names(gene) <- unlist(lapply(gene,function(x) return(x$gene_name)))
        # Get a GRangesList with requested gene strcuture 
        labelHelper <- getGeneCoordinatesForSpline(as.list(gene),dbGene)
        geneList <- getGeneCoordinatesForTrack(as.list(gene),dbGene,dbExon)
    }
    else
        labelHelper <- geneList <- NULL
    
    # Calculate coverage over coordinates
    theCoverage <- vector("list",length(class))
    names(theCoverage) <- class
    if (fromBam) {
        theCoverage <- vector("list",length(class))
        names(theCoverage) <- class
        for (st in class) {
            bsvMessage("Calculating coverage for ",st)
            bsvMessage("Calculating coverage for ",st,
                messageContainer=messageContainer)
            if (is.function(progressFun)) {
                text <- paste("Calculating signal for ",st)
                progressFun(detail=text)
            }
            ind <- which(config$source==source & config$dataset==dataset 
                & config$class==st)
            subconf <- as.character(config$sample_dir[ind])
            names(subconf) <- as.character(config$sample_id[ind])
            theCoverage[[st]] <- cmclapply(subconf,function(x,coords) {
                bam.file <- dir(x,pattern=".bam$",full.names=TRUE)
                bam.index <- dir(x,pattern=".bai$",full.names=TRUE)
                bp <- ScanBamParam(which=coords)
                reads <- unlist(grglist(readGAlignments(file=bam.file,
                    index=bam.index,param=bp,with.which_label=TRUE)))
                seqlevels(reads) <- as.character(seqnames(coords))
                return(calcCoverage(reads,coords,assign.names=FALSE,
                    verbose=FALSE)[[1]])
            },area,rc=rc)
        }
        
        # Normalize
        fac.index <- which(config$source==source & config$dataset==dataset 
            & config$class %in% class)
        factors <- config$norm_factor[fac.index]
        names(factors) <- as.character(config$sample_id[fac.index])
        for (status in names(theCoverage)) {
            bsvMessage("Normalizing ",status," coverage")
            bsvMessage("Normalizing ",status," coverage",
                messageContainer=messageContainer)
            if (is.function(progressFun)) {
                text <- paste("Normalizing ",status)
                progressFun(detail=text)
            }
            for (s in names(theCoverage[[status]])) {
                theCoverage[[status]][[s]] <- theCoverage[[
                    status]][[s]]*factors[s]
            }
        }
    }
    else {
        for (st in class) {
            bsvMessage("Calculating coverage for ",st)
            bsvMessage("Calculating coverage for ",st,
                messageContainer=messageContainer)
            if (is.function(progressFun)) {
                text <- paste("Calculating signal for ",st)
                progressFun(detail=text)
            }
            ind <- which(config$source==source & config$dataset==dataset 
                & config$class==st)
            sc <- as.character(config$track_dir[ind])
            names(sc) <- as.character(config$sample_id[ind])
            theCoverage[[st]] <- cmclapply(names(sc),function(x,crds,sc) {
                chr <- as.character(seqnames(crds))[1]
                    bw.file <- file.path(sc[x],paste(x,"bigWig",sep="."))
                bwrle <- import.bw(bw.file,selection=BigWigSelection(crds),
                    as="RleList")
                if (chr %in% names(bwrle))
                    return(bwrle[[chr]][start(crds):end(crds)])
                else
                    return(NULL)
            },area,sc,rc=rc)
        }
    }
    
    # 5. Construct ggplot data frames
    theAverage <- tcgaLower <- tcgaUpper <- ggData <- 
        vector("list",length(class))
    names(theAverage) <- names(tcgaLower) <- names(tcgaUpper) <- 
        names(ggData) <- class
    for (status in (names(theCoverage))) {
        bsvMessage("Averaging ",status," coverage")
        bsvMessage("Averaging ",status," coverage",
            messageContainer=messageContainer)
        if (is.function(progressFun)) {
            text <- paste("Averaging ",status," signal")
            progressFun(detail=text)
        }
        tmp <- cmclapply(theCoverage[[status]],function(x) {
        if (class(x)=="Rle")
            return(as.numeric(x))
        },rc=rc)
        tmp <- do.call("rbind",tmp)
        if (sumStat=="trimmed")
            theAverage[[status]] <- apply(tmp,2,mean,trim=trim)
        else
            theAverage[[status]] <- apply(tmp,2,sumStat)
        s <- apply(tmp,2,function(x) sd(x)/sqrt(length(x)))
        tcgaLower[[status]] <- theAverage[[status]] - 1.96*s
        tcgaLower[[status]][tcgaLower[[status]]<0] <- 0
        tcgaUpper[[status]] <- theAverage[[status]] + 1.96*s
        ggData[[status]] <- data.frame(
            Position=start(area):end(area),
            Signal=theAverage[[status]],
            ymin=tcgaLower[[status]],
            ymax=tcgaUpper[[status]]
        )
    }
    
    M <- sapply(ggData,function(x) return(max(x$ymax)))
    
    bsvMessage("Creating plot structure")
    bsvMessage("Creating plot structure",messageContainer=messageContainer)
    if (is.function(progressFun)) {
        text <- paste("Creating plot structure")
        progressFun(detail=text)
    }
    
    # Construct tracks
    # Ideogram: TCGA so hg19 by default
    pIdeo <- Ideogram(genome="hg19",subchr=as.character(seqnames(refArea))[1]) + 
        xlim(area)
    # Annotation
    if (!is.null(labelHelper) && !is.null(geneList))
        pAnn <- autoplot(geneList,gap.geom="arrow") +
            geom_text(
                data=data.frame(
                    name=names(gene),
                    x=start(labelHelper)+(end(labelHelper)-
                        start(labelHelper))/2,
                    y=rep(2,length(names(gene)))
                ),
                aes(x=x,y=y,label=name)
            ) +
            theme_bw() +
            theme(
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_blank(),
                axis.ticks.x=element_blank()
            )
    else
        pAnn <- NULL
    # Data
    ggplot.plot <- vector("list",length(class))
    names(ggplot.plot) <- class
    for (cl in class) {
        ggplot.plot[[cl]] <-
            ggplot(ggData[[cl]],mapping=aes(x=Position,y=Signal)) + 
            geom_line(size=0.5,colour=classColours[cl]) +
            geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax),
                alpha=0.3,size=0,fill=classColours[cl]) +
            ylim(0,max(M)+2) +
            ylab(paste(class,"normalized signal")) +
            theme_bw() +
            theme(
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank()
            )
    }
    
    nc <- length(class)
    if (!is.null(pAnn)) {
        if (nc==1)
            ggbio.track <- tracks(
               pIdeo,
               ggplot.plot[[1]],
               pAnn,
               heights=c(0.15,0.75,0.1),
               xlab="Position in bp",
               label.text.cex=1.3
           )
        else if (nc==2)
            ggbio.track <- tracks(
               pIdeo,
               ggplot.plot[[1]],
               ggplot.plot[[2]],
               pAnn,
               heights=c(0.11,0.36,0.36,0.07),
               xlab="Position in bp",
               label.text.cex=1.3
           )
        else if (nc==3)
            ggbio.track <- tracks(
               pIdeo,
               ggplot.plot[[1]],
               ggplot.plot[[2]],
               ggplot.plot[[3]],
               pAnn,
               heights=c(0.09,0.28,0.28,0.28,0.07),
               xlab="Position in bp",
               label.text.cex=1.2
           )
        else if (nc==4)
            ggbio.track <- tracks(
               pIdeo,
               ggplot.plot[[1]],
               ggplot.plot[[2]],
               ggplot.plot[[3]],
               ggplot.plot[[4]],
               pAnn,
               heights=c(0.07,0.22,0.22,0.22,0.22,0.05),
               xlab="Position in bp",
               label.text.cex=1.1
           )
    }
    else {
        if (nc==1)
            ggbio.track <- tracks(
               pIdeo,
               ggplot.plot[[1]],
               heights=c(0.2,0.8),
               xlab="Position in bp",
               label.text.cex=1.3
           )
        else if (nc==2)
            ggbio.track <- tracks(
               pIdeo,
               ggplot.plot[[1]],
               ggplot.plot[[2]],
               heights=c(0.14,0.43,0.43),
               xlab="Position in bp",
               label.text.cex=1.3
           )
        else if (nc==3)
            ggbio.track <- tracks(
               pIdeo,
               ggplot.plot[[1]],
               ggplot.plot[[2]],
               ggplot.plot[[3]],
               heights=c(0.1,0.3,0.3,0.3),
               xlab="Position in bp",
               label.text.cex=1.2
           )
        else if (nc==4)
            ggbio.track <- tracks(
               pIdeo,
               ggplot.plot[[1]],
               ggplot.plot[[2]],
               ggplot.plot[[3]],
               ggplot.plot[[4]],
               heights=c(0.04,0.24,0.24,0.24,0.24),
               xlab="Position in bp",
               label.text.cex=1.1
           )
    }
    
    #if (length(class)==1) {
    #   ggbio.track <- tracks(
    #       chrom=pIdeo,tumor=ggplot.plot.tumor,pAnn,
    #       heights=c(0.2,0.7,0.1),
    #       xlab="Position in bp",
    #       label.text.cex=1.3
    #   )
    #    if (class=="normal")
    #        ggbio.track <- tracks(
    #            chrom=pIdeo,tumor=ggplot.plot.normal,pAnn,
    #            heights=c(0.2,0.7,0.1),
    #            xlab="Position in bp",
    #            label.text.cex=1.3
    #        )
    #}
    #if (length(class)==2)
    #    ggbio.track <- tracks(
    #        chrom=pIdeo,tumor=ggplot.plot.tumor,normal=ggplot.plot.normal,pAnn,
    #        heights=c(0.15,0.34,0.34,0.07),
    #        xlab="Position in bp",
    #        label.text.cex=1.3
    #    )
    
    ggbio.track <- ggbio.track +
        theme(
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(size=10),
            axis.title.y=element_text(size=10),
            axis.text.x=element_text(size=9,face="bold"),
            axis.text.y=element_text(size=9,face="bold")
        )
    return(ggbio.track)
}

getCustomCounts <- function(coords,samples,config,messageContainer=NULL,
    progressFun=NULL,rc=NULL) {
    coords <- validGRangesFromGRanges(coords)
    subconf <- as.character(config[samples,"sample_dir"])
    names(subconf) <- as.character(config[samples,"sample_id"])
    customCounts <- vector("list",length(coords))
    names(customCounts) <- names(coords)
    for (n in names(coords)) {
        bsvMessage("Calculating expression for ",n)
        #bsvMessage("Calculating expression for ",n,
        #    messageContainer=messageContainer)
        if (is.function(progressFun)) {
            text <- paste("Calculating expression for ",n)
            progressFun(detail=text)
        }
        customCounts[[n]] <- unlist(cmclapply(subconf,function(x,coord) {
            bam.file <- dir(x,pattern=".bam$",full.names=TRUE)
            bam.index <- dir(x,pattern=".bai$",full.names=TRUE)
            bp <- ScanBamParam(which=coord)                
            reads <- unlist(grglist(readGAlignments(file=bam.file,
                index=bam.index,param=bp,with.which_label=TRUE)))
            seqlevels(reads) <- as.character(seqnames(coords))
            return(length(reads))
        },coords[n],rc=rc))
        names(customCounts[[n]]) <- names(subconf)
    }
    return(do.call("rbind",customCounts))
}

getGeneCoordinatesForSpline <- function(gene,anndb,flank=NULL) {
    # anndb must be GRanges in this case
    if (is(gene,"GRangesList"))
        stop("The gene argument cannot be a GRangesList in a spline plot!")
    if (!is(anndb,"GRanges"))
        stop("The anndb argument must be a GRanges object!")
    
    preCoords <- vector("list",length(gene))
    for (i in 1:length(gene)) {
        val <- gene[[i]]
        switch(class(val),
            character = {
                valCoords <- tryCatch(anndb[val],
                error=function(e) {
                    m <- match(val,as.character(anndb$gene_name))
                    if (length(m)>0) {
                        w <- width(anndb[m])
                        mm <- m[which(w==max(w))[1]]
                        #valCoords <- anndb[mm]
                        anndb[mm]
                    }
                    else
                        #valCoords <- preCoords[m]
                        preCords[m]
                },finally="")
            },
            GRanges = {
                valCoords <- validGRangesFromGRanges(val)
            }
        )
        preCoords[[i]] <- valCoords
        names(preCoords[[i]]) <- as.character(valCoords$gene_name)
    }
    coords <- do.call("c",preCoords)
    #names(coords) <- names(gene)
    if (!is.null(flank)) {
        w <- width(coords)
        coords <- promoters(coords,upstream=flank[1],downstream=0)
        coords <- resize(coords,width=w+flank[1]+flank[2])
    }
    return(coords)
}

getAreaCoordinatesForTrack <- function(chr,start,end,anndb) {
    if (!is.character(chr))
        stop("The chromosome argument must be a string!")
    if (!is.numeric(start) || !is.numeric(end))
        stop("The start and end arguments must be numeric!")
    if (start<0)
        stop("The start position must be positive!")
    if (!is(anndb,"GRanges") && !is(anndb,"GRangesList"))
        stop("The anndb argument must be a GRanges object!")
    chr <- chr[1]
    start <- start[1]
    end <- end[1]
    if (!(chr %in% seqlevels(anndb)))
        stop("Unknown chromosome: ",chr)
    sl <- seqlengths(anndb)
    if (end>sl[chr])
        stop("The end position exceeds chromosome ",chr," length!")
    # All must be OK
    area <- GRanges(
        seqnames=Rle(chr),
        ranges=IRanges(start=start,end=end),
        strand=Rle("*"),
        name=paste(chr,":",start,"-",end,sep="")
    )
    seqinfo(area) <- seqinfo(anndb)[chr]
    return(area)
}

getGeneCoordinatesForTrack <- function(gene,dbGene,dbExon) {
    if (!is(dbGene,"GRanges"))
        stop("The dbGene argument must be a GRanges object!")
    if (!is(dbExon,"GRangesList"))
        stop("The dbExon argument must be a GRangesList object!")
    
    preNames <- rep(NA,length(gene))
    for (i in 1:length(gene)) {
        val <- gene[[i]]
        switch(class(val),
            character = {
                m <- match(val,as.character(dbGene$gene_name))
                if (length(m)>0) {
                    w <- width(dbGene[m])
                    mm <- m[which(w==max(w))[1]]
                    valName <- as.character(dbGene[mm]$gene_id)
                }
                else
                    valName <- as.character(dbGene[m]$gene_id)
            },
            GRanges = {
                valRange <- validGRangesFromGRanges(val)
                valName <- as.character(valRange[1]$gene_id)
            }#,
            #GRangesList = {
            #    valRange <- validGRangesListFromGRangesList(val[[1]])
            #    valName <- as.character(valRange[[1]]$gene_id[1])
            #}
        )
        preNames[i] <- valName
    }
    na <- which(is.na(preNames))
    if (length(na)>0)
        preNames <- preNames[-na]
    if (!all(preNames %in% names(dbExon))) { 
        # Something is not annotated so we have to return it as is
        o <- which(!(preNames %in% names(dbExon)))
        theNames <- preNames[-o]
        output <- dbExon[match(theNames,names(dbExon))]
        tmp <- gene[o]
        for (k in 1:length(tmp)) {
            if (is(tmp[[k]],"GRanges") || is(tmp[[k]],"GRangesList")) {
                tmpk <- validGRangesListRangeFromGRanges(tmp[[k]])
                if (is(tmpk,"GRanges"))
                    output <- c(output,GRangesList(tmpk))
                else if (is(tmpk,"GRangesList"))
                    output <- c(output,tmpk)
                names(output)[length(output)] <- names(tmpk)
            }
        }
        return(output)
    }
    else
        return(dbExon[preNames])
}

validGRangesFromGRanges <- function(gr) {
    g <- as.data.frame(gr)
    # Must have fields chromosome, start, end, gene_id, gc_content, strand, 
    # gene_name, biotype. Strictly, the first 4 and the 6th.
    if (!all(c("seqnames","start","end","gene_id","strand") %in% names(g)))
        stop("The provided GRanges must have at least these basic fields: ",
            "seqnames, start, end, gene_id, strand")
    df <- g[,c("seqnames","start","end","gene_id","strand")]
    names(df)[1] <- "chromosome"
    if ("gc_content" %in% names(g))
        df$gc_content <- g$gc_content
    else
        df$gc_content <- 0
    if ("gene_name" %in% names(g))
        df$gene_name <- g$gene_name
    else 
        df$gene_name <- g$gene_id
    if ("biotype" %in% names(g))
        df$biotype <- g$biotype
    else 
        df$biotype <- "unknown"
    df <- df[,c("chromosome","start","end","gene_id","gc_content","strand",
        "gene_name","biotype")]
    vgr <- makeGRangesFromDataFrame(
        df=df,
        keep.extra.columns=TRUE
    )
    names(vgr) <- as.character(df$gene_id)
    return(vgr)
}

validGRangesListRangeFromGRanges <- function(gr) {
    g <- as.data.frame(gr)
    # Must have fields chromosome, start, end, gene_id, gc_content, strand, 
    # gene_name, biotype. Strictly, the first 4 and the 6th.
    if (!all(c("seqnames","start","end","gene_id","strand") %in% names(g)))
        stop("The provided GRanges must have at least these basic fields: ",
            "seqnames, start, end, gene_id, strand")
    df <- g[,c("seqnames","start","end","gene_id","strand")]
    names(df)[1] <- "chromosome"
    if ("exon_id" %in% names(g))
        df$exon_id <- g$exon_id
    else
        df$exon_id <- g$gene_id
    if ("gene_name" %in% names(g))
        df$gene_name <- g$gene_name
    else 
        df$gene_name <- g$gene_id
    if ("biotype" %in% names(g))
        df$biotype <- g$biotype
    else 
        df$biotype <- "unknown"
    df <- df[,c("chromosome","start","end","exon_id","gene_id","strand",
        "gene_name","biotype")]
    vgr <- makeGRangesFromDataFrame(
        df=df,
        keep.extra.columns=TRUE
    )
    names(vgr) <- as.character(df$gene_id)
    return(vgr)
}

validateGeneInput <- function(gene) {
    if (!is.list(gene) || (is.list(gene) && is.null(names(gene))))
        stop("gene must be a list of mixed named objects (character, ",
            "GRanges or GRangesList objects!")
    for (n in names(gene)) {
        if (!(is.character(gene[[n]]) && length(gene[[n]])==1)
            && !is(gene[[n]],"GRanges") && !is(gene[[n]],"GRangesList"))
            stop("each gene member must be a character vector of length one ",
                "or a GRanges object or a GRangesList object.")
        if (is(gene[[n]],"GRanges") && length(gene[[n]])>1)
            stop("When gene member is a GRanges, it should have length one!")
    }
}

calcCoverage <- function(input,mask,strand=NULL,ignore.strand=TRUE,
    assign.names=TRUE,verbose=TRUE) {
    if (!is(input,"GRanges"))
        stop("The input argument must be a GenomicRanges object")
    if (!is(mask,"GRanges"))
        stop("The mask argument must be a GenomicRanges object")
    if (!is.null(strand) && !is.list(strand)) {
        message("Retrieving ",strand," reads...")
        input <- input[strand(input)==strand]
    }
    if (!is.list(input))
        input <- splitBySeqname(input,verbose=verbose)
    index <- 1:length(mask)
    if (verbose) 
        message("Calculating coverage...")
    coverage <- cmclapply(index,function(i,mask,input,ignore.strand) {
        x <- mask[i]
        y<-list(
            chromosome=as.character(seqnames(x)),
            start=start(x),
            end=end(x),
            strand=as.character(strand(x)),
            reads=NULL,
            coverage=NULL
        )
        if (!is.null(input[[y$chromosome]])) {
            y$reads <- input[[y$chromosome]][
                subjectHits(findOverlaps(x,input[[y$chromosome]],
                    ignore.strand=ignore.strand))]
        }
        else {
            if (verbose) message(y$chromosome,"not found!")
            y$reads <- NULL
        }
        if (length(y$reads)>0) {
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            y$coverage <- y$coverage[[cc]][y$start:y$end]
            if (y$strand=="+")
                return(y$coverage)
            else if (y$strand=="-")
                return(rev(y$coverage))
            else
                return(y$coverage)
        }
        else
            return(NULL)
    },mask,input,ignore.strand)
    if (assign.names)
        names(coverage) <- names(mask)
    gc(verbose=FALSE)
    if (verbose) 
        message("Done!")
    return(coverage) # Rle
}

calcCoverageFromBW <- function(input,mask,assign.names=TRUE,verbose=TRUE) {
    if (!is(input,"RleList"))
        stop("The input argument must be an RleList object")
    if (!is(mask,"GRanges"))
        stop("The mask argument must be a GenomicRanges object")
    index <- 1:length(mask)
    if (verbose) 
        message("Calculating coverage...")
    coverage <- cmclapply(index,function(i,mask,input) {
        x <- mask[i]
        y<-list(
            chromosome=as.character(seqnames(x)),
            start=start(x),
            end=end(x),
            strand=as.character(strand(x)),
            coverage=NULL
        )
        if (!is.null(input[[y$chromosome]])) {
            y$coverage <- input[[y$chromosome]][y$start:y$end]
            if (y$strand=="+")
                return(y$coverage)
            else if (y$strand=="-")
                return(rev(y$coverage))
            else
                return(y$coverage)
        }
        else {
            if (verbose) message(y$chromosome,"not found!")
            y$coverage <- NULL
        }
    },mask,input)
    if (assign.names)
        names(coverage) <- names(mask)
    if (verbose) 
        message("Done!")
    return(coverage) # Rle
}

loadJBrowse <- function(source,dataset,config,org="hg19") {
    #urlBase <- paste("https://www.google.com/search?q=%",
    #    "http://epigenomics.fleming.gr/seqcbrowse/index.html?",
    #    "&btnI=Im+Feeling+Lucky",sep="")
    urlBase <- "https://epigenomics.fleming.gr/seqcbrowse/index.html?"
    tracksBase <- paste("https://epigenomics.fleming.gr/seqcvibe_tracks",
        org,sep="/")
    
    ind <- which(as.character(config$source)==source 
        & as.character(config$dataset)==dataset)
    subconf <- config[ind,]
    
    if (!is.null(subconf$alt_id))
        initTracks <- paste(as.character(subconf$sample_id),
            as.character(subconf$alt_id),"xy",sep="_")
    else
        initTracks <- paste(as.character(subconf$sample_id),"xy",sep="_")
    
    initTracks <- c(initTracks,"ucsc_human_hg19","ensGene","refGene")
    initTracks <- paste(initTracks,collapse=",")
    
    query <- paste("data=",tracksBase,"&tracks=",initTracks,sep="");
    
    return(paste(urlBase,query,sep=""))
}
