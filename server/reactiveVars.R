initReactiveVars <- function() {
    currentMetadata <- reactiveValues(
        source=sources[1],
        dataset=datasets[1],
        class=classes,
        genome=genomes[1],
        metadata=NULL,
        final=NULL
    )
    
    currentGenes <- reactiveValues(
        genes=NULL,
        coords=list(
            chromosome=NULL,
            start=NULL,
            end=NULL,
            strand=NULL,
            name=NULL
        )
    )
    
    customRegions <- reactiveValues(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    currentOpts <- reactiveValues(
        flank=c(2000,2000),
        sumStat="mean",
        trim=0.1,
        colours=NULL
    )
    
    genePlots <- reactiveValues(
        geneProfile=ggmessage("Gene profiles will\nbe displayed here"),
        rendered=TRUE
    )
    
    customArea <- reactiveValues(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    customRegionsInArea <- reactiveValues(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    areaPlots <- reactiveValues(
        areaProfile=ggmessage("Area profiles will\nbe displayed here"),
        rendered=FALSE
    )
    
    currentTables <- reactiveValues()
    
    customRnaRegions <- reactiveValues(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    currentCustomRnaTables <- reactiveValues(
        tables=list(),
        displayTables=list(),
        lengths=NULL
    )
    
    maPlots <- reactiveValues(
        #maPlot=data.frame(A=1,M=1,Status=1,
        #    Gene="Resulting MA plots will be displayed here"),
        maPlot=ggmessage("Resulting MA plots will\nbe displayed here"),
        maData=NULL,
        maColours=list(
            Up="#B40000",
            Down="#00B400",
            Neutral="#6B6B6B"
        ),
        maZoom=list(
            x=NULL,
            y=NULL
        ),
        rendered=TRUE
    )
    
    currentPipelineOutput <- reactiveValues(
        annotation=NULL,
        counts=NULL,
        flags=NULL,
        classList=NULL,
        contrastList=NULL,
        pValue=NULL,
        fdr=NULL
    )
    
    currentRnaDeTable <- reactiveValues(
        totalTable=NULL,
        tableFilters=list(
            p=0.05,
            fdr=0.05,
            scale="natural",
            fc=c(0.5,2),
            bt=NULL,
            genes=NULL,
            chr=NULL
        )
    )
    
    currentHeatmap <- reactiveValues(
        entry=ggmessage("Resulting heatmap will\nbe displayed here"),
        timeout=ggmessage(paste("Clustering operation took too long\nto ",
            "complete and was aborted.\nConsider lowering the number of genes",
            "\nand/or conditions to prevent this.",sep=""),type="error"),
        error=ggmessage(paste("Clustering operation resulted in an error\n ",
            "most probably because of memory reasons.\nConsider lowering the ",
            "number of genes\nand/or conditions to prevent this.",sep=""),
            type="error"),
        data=NULL,
        opts=list(
            dendrogram="both",
            distfun="dist",
            hclustfun="hclust",
            Rowv=TRUE,
            Colv=TRUE,
            k_row=1,
            k_col=1,
            colors="RdYlBu"
        )
    )
    
    currentCorrelation <- reactiveValues(
        entryCor=ggmessage("Resulting figures will\nbe displayed here"),
        errorCor=ggmessage(paste("Correlation analysis has failed.\nThe most ",
            "probable reason is that the\nselected gene set does not show\n",
            "enough variability to produce a\ncorrelation matrix. Select ",
            "another gene\nset and try again.",sep=""),type="error"),
        entryMds=ggmessage("MDS sample will\nbe displayed here",size="small"),
        warnMds=ggmessage("MDS failed! Try\nchanging some parameters.",
            type="warning",size="small"),
        errorMds=ggmessage("Error creating MDS!\nSee error on the right",
            type="error",size="small"),
        corMatrix=NULL,
        datMatrix=NULL,
        mdsRsq=NULL,
        refGene=NULL,
        what="samples",
        opts=list(
            method="pearson",
            symm=FALSE,
            colors=c("#FFFF00","#BEBEBE","#0000FF")
        ),
        tooManyGenesAlert=FALSE
    )
    
    currentDimRed <- reactiveValues(
        datMatrix=NULL,
        selMatrix=NULL,
        mdsObj=NULL,
        mdsData=NULL,
        pcaScoresData=NULL,
        pcaLoadingsData=NULL,
        pcaRankedLoadingsData=NULL,
        pcaObj=NULL,
        mdsGof=list(
            dist=NULL,
            rsq=NULL,
            gof=NULL
        ),
        opts=list(
            method=NULL,
            colors=NULL
        ),
        mdsPlot=ggmessage("Resulting MDS plots will\nbe displayed here"),
        pcaScreePlot=
            ggmessage("Resulting PCA scree plots\nwill be displayed here"),
        pcaScoresPlot=
            ggmessage("Resulting PCA score plots\nwill be displayed here"),
        pcaLoadingsPlot=
            ggmessage("Resulting PCA loadings plots\nwill be displayed here"),
        pcaRankedLoadingsPlot=
            ggmessage(paste("Resulting PCA ranked loadings plots\nwill be ",
                "displayed here",sep="")),
        pcaBiplotPlot=
            ggmessage("Resulting PCA biplot plots\nwill be displayed here"),
        tooManyGenes=FALSE # Safety trigger for plotting names with ggrepel
    )
    
    return(list(
        currentMetadata=currentMetadata,
        currentGenes=currentGenes,
        customRegions=customRegions,
        currentOpts=currentOpts,
        genePlots=genePlots,
        customArea=customArea,
        customRegionsInArea=customRegionsInArea,
        areaPlots=areaPlots,
        currentTables=currentTables,
        customRnaRegions=customRnaRegions,
        currentCustomRnaTables=currentCustomRnaTables,
        maPlots=maPlots,
        currentPipelineOutput=currentPipelineOutput,
        currentRnaDeTable=currentRnaDeTable,
        currentDimRed=currentDimRed,
        currentHeatmap=currentHeatmap,
        currentCorrelation=currentCorrelation
    ))
}

initReactiveMsgs <- function() {
    dataSelectorMessages <- reactiveValues(
        messages=list(
            list(
                type="INFO",
                msg=paste(getTime("INFO"),"Welcome to the data selector of ",
                    "SeqCVIBE! This is an info message. Make your selections ",
                    "on the left.")
            )
        )
    )
    
    geneExplorerMessages <- reactiveValues(
        messages=list(
            list(
                type="INFO",
                msg=paste(getTime("INFO"),"Welcome to the gene signal ",
                    "explorer of SeqCVIBE! This is an info message. Make ",
                    "your selection on the left.")
            )
        )
    )
    
    areaExplorerMessages <- reactiveValues(
        messages=list(
            list(
                type="INFO",
                msg=paste(getTime("INFO"),"Welcome to the area signal ",
                    "explorer of SeqCVIBE! This is an info message. Make ",
                    "your selection on the left.")
            )
        )
    )
    
    return(list(
        dataSelectorMessages=dataSelectorMessages,
        geneExplorerMessages=geneExplorerMessages,
        areaExplorerMessages=areaExplorerMessages
    ))
}

clearReactiveVars <- function(allReactiveVars) {
    allReactiveVars$currentGenes$genes <- NULL
    allReactiveVars$currentGenes$coords <- list(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    allReactiveVars$customRegions$chromosome <- NULL
    allReactiveVars$customRegions$start <- NULL
    allReactiveVars$customRegions$end <- NULL
    allReactiveVars$customRegions$strand <- NULL
    allReactiveVars$customRegions$name <- NULL
    
    allReactiveVars$genePlots$geneProfile <- 
        ggmessage("Gene profiles will\nbe displayed here")
    allReactiveVars$genePlots$rendered <- TRUE
    
    allReactiveVars$customArea$chromosome <- NULL
    allReactiveVars$customArea$start <- NULL
    allReactiveVars$customArea$end <- NULL
    allReactiveVars$customArea$strand <- NULL
    allReactiveVars$customArea$name <- NULL
    
    allReactiveVars$customArea$customRegionsInArea <- NULL
    allReactiveVars$customArea$chromosome <- NULL
    allReactiveVars$customArea$start <- NULL
    allReactiveVars$customArea$end <- NULL
    allReactiveVars$customArea$strand <- NULL
    allReactiveVars$customArea$name <- NULL
    
    allReactiveVars$areaPlots$areaProfile <-
        ggmessage("Area profiles will\nbe displayed here")
    allReactiveVars$areaPlots$rendered=FALSE
    
    allReactiveVars$currentTables <- reactiveValues()
    
    allReactiveVars$customRnaRegions$chromosome <- NULL
    allReactiveVars$customRnaRegions$start <- NULL
    allReactiveVars$customRnaRegions$end <- NULL
    allReactiveVars$customRnaRegions$strand <- NULL
    allReactiveVars$customRnaRegions$name <- NULL
    
    allReactiveVars$currentCustomRnaTables$tables <- list()
    allReactiveVars$currentCustomRnaTables$lengths <- NULL
    
    allReactiveVars$maPlots$maPlot <- 
        ggmessage("Resulting MA plots will\nbe displayed here")
    allReactiveVars$maPlots$maData <- NULL
    allReactiveVars$maPlots$maColours <- list(
        Up <- "#B40000",
        Down <- "#00B400",
        Neutral <- "#6B6B6B"
    )
    allReactiveVars$maPlots$maZoom <- list(
        x <- NULL,
        y <- NULL
    )
    allReactiveVars$maPlots$rendered <- TRUE
    
    allReactiveVars$currentPipelineOutput$annotation <- NULL
    allReactiveVars$currentPipelineOutput$counts <- NULL
    allReactiveVars$currentPipelineOutput$flags <- NULL
    allReactiveVars$currentPipelineOutput$classList <- NULL
    allReactiveVars$currentPipelineOutput$contrastList <- NULL
    allReactiveVars$currentPipelineOutput$pValue <- NULL
    allReactiveVars$currentPipelineOutput$fdr <- NULL
    
    allReactiveVars$currentRnaDeTable$totalTable <- NULL
    allReactiveVars$currentRnaDeTable$tableFilters <- list(
        p=0.05,
        fdr=0.05,
        scale="natural",
        fc=c(0.5,2),
        bt=NULL,
        genes=NULL,
        chr=NULL
    )
    
    allReactiveVars$currentHeatmap$entry <- 
        ggmessage("Resulting heatmap will\nbe displayed here")
    allReactiveVars$currentHeatmap$data <- NULL
    allReactiveVars$currentHeatmap$opts <- list(
        dendrogram="both",
        distfun="dist",
        hclustfun="hclust",
        Rowv=TRUE,
        Colv=TRUE,
        k_row=1,
        k_col=1,
        colors="RdYlBu"
    )
    
    allReactiveVars$currentCorrelation$entryCor <-
        ggmessage("Resulting figures will\nbe displayed here")
    allReactiveVars$currentCorrelation$entryMds <- 
        ggmessage("MDS sample will\nbe displayed here",size="small")
    allReactiveVars$currentCorrelation$corMatrix <- NULL
    allReactiveVars$currentCorrelation$datMatrix <- NULL
    allReactiveVars$currentCorrelation$what <- "samples"
    allReactiveVars$currentCorrelation$opts <- list(
        method="pearson",
        symm=FALSE,
        colors=c("#FFFF00","#BEBEBE","#0000FF")
    )
    allReactiveVars$currentCorrelation$tooManyGenesAlert=FALSE
        
    allReactiveVars$currentDimRed$datMatrix <- NULL
    allReactiveVars$currentDimRed$selMatrix <- NULL
    allReactiveVars$currentDimRed$mdsObj <- NULL
    allReactiveVars$currentDimRed$pcaObj <- NULL
    allReactiveVars$currentDimRed$mdsData <- NULL
    allReactiveVars$currentDimRed$pcaScoresData <- NULL
    allReactiveVars$currentDimRed$pcaLoadingsData <- NULL
    allReactiveVars$currentDimRed$pcaLoadingsData <- NULL
    allReactiveVars$currentDimRed$mdsGof <- list(
        dist=NULL,
        rsq=NULL,
        gof=NULL
    )
    allReactiveVars$currentDimRed$opts <- list(
        method=NULL,
        colors=NULL
    )
    allReactiveVars$currentDimRed$mdsPlot <- 
        ggmessage("Resulting MDS plots will\nbe displayed here")
    allReactiveVars$currentDimRed$pcaScreePlot <-
            ggmessage("Resulting PCA scree plots\nwill be displayed here")
    allReactiveVars$currentDimRed$pcaScoresPlot <-
            ggmessage("Resulting PCA score plots\nwill be displayed here")
    allReactiveVars$currentDimRed$pcaLoadingsPlot <-
            ggmessage("Resulting PCA loadings plots\nwill be displayed here")
    allReactiveVars$currentDimRed$pcaRankedLoadingsPlot <-
            ggmessage(paste("Resulting PCA ranked loadings plots\nwill be ",
                "displayed here",sep=""))
    allReactiveVars$currentDimRed$pcaBiplotPlot <-
            ggmessage("Resulting PCA biplots plots\nwill be displayed here")
    allReactiveVars$currentDimRed$tooManyGenes <- FALSE
    
    return(allReactiveVars)
}

ggmessage <- function(msg="",type=c("generic","info","success","warning",
    "error"),size=c("large","small")) {
    type <- tolower(type[1])
    size <- tolower(size[1])
    switch(type,
        generic = { color <- "black" },
        info = { color <- "green2" },
        success = { color <- "blue2" },
        warning = { color <- "orange" },
        error = { color <- "red2" }
    )
    switch(size,
        large = { s <- 10 },
        small = { s <- 5 }
    )
    return(
        ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,label=msg),
                aes(x=x,y=y,label=label),colour=color,size=s) +
            theme(
                axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank()
            )
    )
}

## Messages boilerplate
#dataSelectorMessages <- reactiveValues(
#   messages=list(
#       list(
#           type="INFO",
#           msg=paste(getTime("INFO"),"Example info message.")
#       ),
#       list(
#           type="SUCCESS",
#           msg=paste(getTime("SUCCESS"),"Example success message.")
#       ),
#       list(
#           type="WARNING",
#           msg=paste(getTime("WARNING"),"Example warning message.")
#       ),
#       list(
#           type="ERROR",
#           msg=paste(getTime("ERROR"),"Example error message.")
#       )
#   )
#)
