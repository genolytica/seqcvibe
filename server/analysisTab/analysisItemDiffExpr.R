diffExprTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    customRnaRegions <- allReactiveVars$customRnaRegions
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    currentPipelineOutput <- allReactiveVars$currentPipelineOutput
    maPlots <- allReactiveVars$maPlots
    
    runPipeline <- eventReactive(input$performDeAnalysis,{
        require(metaseqR)
            
        # Set up count data table with embedded annotation
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        meta <- currentMetadata$final
        samples <- as.character(meta$sample_id)
        dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
        
        # Check if we have counts from custom regions
        addAnn <- NULL
        if (input$includeCustomRegions) {
            if (!is.null(customRnaRegions$name)) {
                addAnn <- data.frame(
                    chromosome=as.character(customRnaRegions$chromosome),
                    start=as.integer(customRnaRegions$start),
                    end=as.integer(customRnaRegions$end),
                    gene_id=as.character(customRnaRegions$name),
                    gc_content=rep(0,length(customRnaRegions$name)),
                    strand=as.character(customRnaRegions$strand),
                    gene_name=as.character(customRnaRegions$name),
                    biotype=rep("unknown",length(customRnaRegions$name))
                )
                rownames(addAnn) <- addAnn$gene_id
            }
            if (!is.null(currentCustomRnaTables$lengths)) {
                addAnn <- cbind(addAnn,currentCustomRnaTables$lengths)
                names(addAnn)[ncol(addAnn)] <- "active_length"
            }
        }
        
        # Construct annotation
        ann <- data.frame(
            chromosome=as.character(seqnames(dbGene)),
            start=start(dbGene),
            end=end(dbGene),
            gene_id=as.character(dbGene$gene_id),
            gc_content=rep(0,length(dbGene)),
            strand=as.character(strand(dbGene)),
            gene_name=as.character(dbGene$gene_name),
            biotype=as.character(dbGene$biotype)
        )
        rownames(ann) <- ann$gene_id
        
        # Counts
        M <- loadedData[[s]][[d]]$counts[rownames(ann),samples]
        len <- loadedData[[s]][[d]]$length[rownames(ann)]
        ann <- cbind(ann,len)
        names(ann)[ncol(ann)] <- "active_length"
        
        # Final annotation
        ann <- rbind(ann,addAnn)
        
        # Custom counts
        A <- NULL
        if (input$includeCustomRegions) {
            if (!is.null(currentCustomRnaTables$lengths)) {
                A <- do.call("cbind",currentCustomRnaTables$tables)
                A <- A[,colnames(M),drop=FALSE]
            }
        }
        M <- rbind(M,A)
        
        # Final feed to metaseqr
        D <- cbind(ann,M)
        
        # Set up sample and contrast list for metaseqr
        cc <- unique(as.character(meta$class))
        classList <- vector("list",length(cc))
        names(classList) <- cc
        for (cl in cc)
            classList[[cl]] <- 
                as.character(meta$sample_id[which(meta$class==cl)])
        control <- input$rnaDePipelineControl
        treatments <- setdiff(cc,control)
        contrast <- c(treatments,control)
        contrast <- paste(contrast,collapse="_vs_")
        
        # Set up normalization
        if (input$rnaDePipeline=="deseq")
            norm <- "deseq"
        else
            norm <- "edger"
        
        # Set up filters
        exprFilt <- list(median=TRUE,mean=FALSE,quantile=NA,known=NA)
        switch(input$rnaDeGeneFilter,
            median = {
                exprFilt <- list(
                    median=TRUE,
                    mean=FALSE,
                    quantile=NA,
                    known=NA
                )
            },
            mean = {
                exprFilt <- list(
                    median=FALSE,
                    mean=TRUE,
                    quantile=NA,
                    known=NA
                )
            },
            quantile = {
                qq <- as.numeric(input$rnaDeQuantileFilter)
                if (!is.na(qq))
                    exprFilt <- list(
                        median=FALSE,
                        mean=FALSE,
                        quantile=qq,
                        known=NA
                    )
            },
            known = {
                if (!isEmpty(input$rnaDeKnownFilter))
                    exprFilt <- list(
                        median=FALSE,
                        mean=FALSE,
                        quantile=NA,
                        known=input$rnaDeKnownFilter
                    )
            }
        )
        btFilt <- NULL
        if (input$rnaDeBiotypeFilter) {
            bts <- getBiotypes(currentMetadata$genome)
            btFilt <- lapply(bts,function(b) {
                return(input[[b]])
            })
            names(btFilt) <- bts
        }
        geneFilters=list(
            length=list(
                length=as.numeric(input$rnaDeGeneLengthFilterValue)
            ),
            avg.reads=list(
                average.per.bp=100,
                quantile=0.25
            ),
            expression=exprFilt, 
            biotype=btFilt
        )
        
        progress <- shiny::Progress$new()
        progress$initialize(
            session,
            min=0,
            max=6
        )
        progress$set(message="",value=0)
        on.exit(progress$close())
        
        updateProgress <- function(value=NULL,detail=NULL) {
            if (is.null(value)) {
                value <- progress$getValue()
                value <- value + 1
            }
            progress$set(value=value,detail=detail)
        }
        
        pipOutput <- tryCatch(
            metaseqr(
                counts=D,
                sample.list=classList,
                contrast=contrast,
                annotation="embedded",
                id.col=4,
                gc.col=5,
                name.col=7,
                bt.col=8,
                org=currentMetadata$genome,
                count.type="gene",
                normalization=norm,
                statistics=input$rnaDePipeline,
                adjust.method=input$rnaDeMTC,
                fig.format="png",
                export.where=file.path(tempdir(),paste("analysis_",
                        format(Sys.time(),format="%Y%m%d%H%M%S"),sep="")),
                qc.plots="mds",
                exon.filters=NULL,
                gene.filters=geneFilters,
                when.apply.filter=input$rnaDeNormalizeWhen,
                export.what=c("annotation","p.value","adj.p.value","counts",
                    "flags"),
                export.scale="natural",
                export.values="normalized",
                export.stats="mean",
                save.gene.model=FALSE,
                report=FALSE,
                out.list=TRUE,
                progress.fun=updateProgress
            ),
        error=function(e) {
            #print(e)
        },
        finally="")
        
        currentPipelineOutput$annotation <- pipOutput$complete$gene.data
        currentPipelineOutput$counts <- pipOutput$complete$norm.counts
        currentPipelineOutput$flags <- pipOutput$complete$flags
        currentPipelineOutput$classList <- pipOutput$complete$sample.list
        currentPipelineOutput$contrastList <- pipOutput$complete$contrast
        currentPipelineOutput$pValue <- pipOutput$complete$p.value
        currentPipelineOutput$fdr <- pipOutput$complete$fdr
    })
    
    updateMaZoomCoords <- eventReactive(input$rnaDeMAPlotDblClick,{
        if (input$toggleRnaDeZoom) {
            brush <- input$rnaDeMAPlotBrush
            if (!is.null(brush)) {
                maPlots$maZoom$x <- c(brush$xmin,brush$xmax)
                maPlots$maZoom$y <- c(brush$ymin,brush$ymax)
            }
        }
        else {
            maPlots$maZoom$x <- NULL
            maPlots$maZoom$y <- NULL
        }
    })
    
    resetMaZoom <- eventReactive(input$resetRnaDeZoom,{
        if (!is.null(maPlots$maZoom$x))
            maPlots$maZoom$x <- NULL
        if (!is.null(maPlots$maZoom$y))
            maPlots$maZoom$y <- NULL
    })
    
    return(list(
        runPipeline=runPipeline,
        updateMaZoomCoords=updateMaZoomCoords,
        resetMaZoom=resetMaZoom
    ))
}

diffExprTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentPipelineOutput <- allReactiveVars$currentPipelineOutput
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    maPlots <- allReactiveVars$maPlots
    
    rnaDeTotalTable <- reactive({
        s <- unique(as.character(currentMetadata$source))
        d <- unique(as.character(currentMetadata$dataset))
        
        ann <- currentPipelineOutput$annotation
        counts <- currentPipelineOutput$counts
        flags <- currentPipelineOutput$flags
        classList <- currentPipelineOutput$classList
        contrastList <- currentPipelineOutput$contrastList
        pValue <- currentPipelineOutput$pValue
        fdr <- currentPipelineOutput$fdr
        
        if (!is.null(counts) && !isEmpty(input$rnaDeCurrentContrast)) {
            p <- pValue[[contrastList[1]]][,1,drop=FALSE]
            fdr <- fdr[[contrastList[1]]][,1,drop=FALSE]
            
            filters <- currentRnaDeTable$tableFilters
            filterInds <- list(
                stat=rownames(counts),
                fold=rownames(counts),
                chr=rownames(counts),
                bt=rownames(counts)
            )
            
            # If there are any genes in the gene selectize box, ignore all other
            # filters
            if (!is.null(filters$genes)) {
                g <- filters$genes
                tab <- counts[g,,drop=FALSE]
                sumTable <- data.frame(
                    pvalue=p[g,,drop=FALSE],
                    fdr=fdr[g,,drop=FALSE]
                )
                names(sumTable) <- c("pvalue","fdr")
                
                switch(input$rnaDeValueCompRadio,
                    counts = {
                        tab <- tab
                    },
                    rpkm = {
                        len <- loadedData[[s]][[d]]$length[g]
                        libsize <- 
                            unlist(loadedData[[s]][[d]]$libsize[colnames(tab)])
                        tab <- round(edgeR::rpkm(tab,gene.length=len,
                            lib.size=libsize),digits=6)
                    },
                    rpgm = {
                        len <- loadedData[[s]][[d]]$length[g]
                        for (j in 1:ncol(tab))
                            tab[,j] <- round(tab[,j]/len,digits=6)
                    }
                )
                switch(input$rnaDeValueScaleRadio,
                    natural = {
                        tab <- tab
                    },
                    log2 = {
                        tab <- round(log2(tab+1),digits=6)
                    }
                )
                
                fcMat <- round(makeFoldChange(contrastList[1],classList,tab,
                    input$rnaDeValueScaleRadio),6)
            
                avgMatrix <- do.call("cbind",lapply(classList,
                function(x,tab,s,v) {
                    makeStat(x,tab,s,v)
                },tab,input$rnaDeValueAverageRadio,input$rnaDeValueCompRadio))
                colnames(avgMatrix) <- paste(names(classList),
                    input$rnaDeValueAverageRadio,sep="_")
                
                stdMatrix <- do.call("cbind",lapply(classList,
                function(x,tab,s,v) {
                    makeStat(x,tab,s,v)
                },tab,input$rnaDeValueDeviationRadio,input$rnaDeValueCompRadio))
                colnames(stdMatrix) <- paste(names(classList),
                    input$rnaDeValueDeviationRadio,sep="_")
                
                totalTable <- cbind(
                    ann <- ann[g,,drop=FALSE],
                    sumTable,
                    fcMat,
                    as.data.frame(avgMatrix),
                    as.data.frame(stdMatrix),
                    as.data.frame(flags[g,,drop=FALSE])
                )               
                currentRnaDeTable$totalTable <- totalTable
                return()
            }
            
            # If interacting with the MA plot
            if (input$toggleRnaDeTableUpdate) {
                g <- NULL
                # If from click
                tabClick <- 
                    nearPoints(maPlots$maData,input$rnaDeMAPlotClick,
                        xvar="A",yvar="M",allRows=TRUE)
                selcl <- which(tabClick$selected_)
                if (length(selcl)>0)
                    g <- rownames(tabClick[selcl,])
                
                # If from brush
                tabBrush <- 
                    brushedPoints(maPlots$maData,input$rnaDeMAPlotBrush,
                        xvar="A",yvar="M",allRows=TRUE)
                selbr <- which(tabBrush$selected_)
                if (length(selbr)>0)
                    g <- rownames(tabBrush[selbr,])
                
                if (isEmpty(g))
                    return()
                
                tab <- counts[g,,drop=FALSE]
                sumTable <- data.frame(
                    pvalue=p[g,,drop=FALSE],
                    fdr=fdr[g,,drop=FALSE]
                )
                names(sumTable) <- c("pvalue","fdr")
                
                switch(input$rnaDeValueCompRadio,
                    counts = {
                        tab <- tab
                    },
                    rpkm = {
                        len <- loadedData[[s]][[d]]$length[g]
                        libsize <- 
                            unlist(loadedData[[s]][[d]]$libsize[colnames(tab)])
                        tab <- round(edgeR::rpkm(tab,gene.length=len,
                            lib.size=libsize),digits=6)
                    },
                    rpgm = {
                        len <- loadedData[[s]][[d]]$length[g]
                        for (j in 1:ncol(tab))
                            tab[,j] <- round(tab[,j]/len,digits=6)
                    }
                )
                switch(input$rnaDeValueScaleRadio,
                    natural = {
                        tab <- tab
                    },
                    log2 = {
                        tab <- round(log2(tab+1),digits=6)
                    }
                )
                
                fcMat <- round(makeFoldChange(contrastList[1],classList,tab,
                    input$rnaDeValueScaleRadio),6)
            
                avgMatrix <- do.call("cbind",lapply(classList,
                function(x,tab,s,v) {
                    makeStat(x,tab,s,v)
                },tab,input$rnaDeValueAverageRadio,input$rnaDeValueCompRadio))
                colnames(avgMatrix) <- paste(names(classList),
                    input$rnaDeValueAverageRadio,sep="_")
                
                stdMatrix <- do.call("cbind",lapply(classList,
                function(x,tab,s,v) {
                    makeStat(x,tab,s,v)
                },tab,input$rnaDeValueDeviationRadio,input$rnaDeValueCompRadio))
                colnames(stdMatrix) <- paste(names(classList),
                    input$rnaDeValueDeviationRadio,sep="_")
                
                totalTable <- cbind(
                    ann <- ann[g,,drop=FALSE],
                    sumTable,
                    fcMat,
                    as.data.frame(avgMatrix),
                    as.data.frame(stdMatrix),
                    as.data.frame(flags[g,,drop=FALSE])
                )               
                currentRnaDeTable$totalTable <- totalTable
                return()
            }
            
            # Gradual application of filters for faster rendering
            # 1. Chromosomes
            tmp <- rownames(counts)
            if (!is.null(filters$chr)) {
                    tmp <- which(as.character(ann$chromosome) %in% filters$chr)
                if (length(tmp)>0)
                    filterInds$chr <- rownames(counts)[tmp]
            }
            # 2. Biotypes
            if (!is.null(filters$bt)) {
                    tmp <- which(as.character(ann$biotype) %in% filters$bt)
                if (length(tmp)>0)
                    filterInds$bt <- rownames(counts)[tmp]
            }
            
            # 3. Statistical score
            if (input$statThresholdType=="pvalue")
                tmp <- names(which(p[,1]<filters$p))
            else if (input$statThresholdType=="fdr") {
                tmp <- names(which(fdr[,1]<filters$fdr))
            }
            if (length(tmp)>0)
                filterInds$stat <- tmp
                
            # Gather so far
            sofar <- Reduce("intersect",filterInds[c("stat","chr","bt")])
            sumTable <- data.frame(
                pvalue=p[sofar,,drop=FALSE],
                fdr=fdr[sofar,,drop=FALSE]
            )
            names(sumTable) <- c("pvalue","fdr")
            ann <- ann[sofar,]
            flags <- flags[sofar,]
            
            # Proceed with count table procesing
            tab <- counts[sofar,,drop=FALSE]
            switch(input$rnaDeValueCompRadio,
                counts = {
                    tab <- tab
                },
                rpkm = {
                    len <- loadedData[[s]][[d]]$length[sofar]
                    libsize=unlist(loadedData[[s]][[d]]$libsize[colnames(tab)])
                    tab <- round(edgeR::rpkm(tab,gene.length=len,
                        lib.size=libsize),digits=6)
                },
                rpgm = {
                    len <- loadedData[[s]][[d]]$length[sofar]
                    for (j in 1:ncol(tab))
                        tab[,j] <- round(tab[,j]/len,digits=6)
                }
            )
            switch(input$rnaDeValueScaleRadio,
                natural = {
                    tab <- tab
                },
                log2 = {
                    tab <- round(log2(tab+1),digits=6)
                }
            )
            
            fcMat <- round(makeFoldChange(contrastList[1],classList,tab,
                input$rnaDeValueScaleRadio),6)
            
            tmp <- names(which(apply(fcMat[,input$rnaDeCurrentContrast,
                drop=FALSE],1,function(x,f) {
                return(any(x<=f[1] | x>=f[2]))
            },filters$fc)))
            
            if (length(tmp)>0)
                filterInds$fold <- tmp
            else
                filterInds$fold <- sofar
                
            tab <- tab[filterInds$fold,,drop=FALSE]
            fcMat <- fcMat[filterInds$fold,,drop=FALSE]
            
            avgMatrix <- do.call("cbind",lapply(classList,function(x,tab,s,v) {
                makeStat(x,tab,s,v)
            },tab,input$rnaDeValueAverageRadio,input$rnaDeValueCompRadio))
            colnames(avgMatrix) <- paste(names(classList),
                input$rnaDeValueAverageRadio,sep="_")
            
            stdMatrix <- do.call("cbind",lapply(classList,function(x,tab,s,v) {
                makeStat(x,tab,s,v)
            },tab,input$rnaDeValueDeviationRadio,input$rnaDeValueCompRadio))
            colnames(stdMatrix) <- paste(names(classList),
                input$rnaDeValueDeviationRadio,input$rnaDeValueCompRadio,
                sep="_")
            
            totalTable <- cbind(
                ann[filterInds$fold,],
                sumTable[filterInds$fold,],
                fcMat,
                as.data.frame(avgMatrix),
                as.data.frame(stdMatrix),
                as.data.frame(flags[filterInds$fold,])
            )
            rownames(totalTable) <- as.character(ann[filterInds$fold,"gene_id"])
            
            currentRnaDeTable$totalTable <- totalTable
        }
    })
    
    handleRnaDeAnalysisSummarySelection <- reactive({
        observeEvent(input$clearRnaDeSummarySelection,{
            proxy <- dataTableProxy("rnaDeAnalysisSummaryTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaDeSummarySelection,{
            N <- input$rnaDeAnalysisSummaryTable_rows_all
            sel <- input$rnaDeAnalysisSummaryTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaDeAnalysisSummaryTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaDeAnalysisAnnotationSelection <- reactive({
        observeEvent(input$clearRnaDeAnnotationSelection,{
            proxy <- dataTableProxy("rnaDeAnalysisAnnotationTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaDeAnnotationSelection,{
            N <- input$rnaDeAnalysisAnnotationTable_rows_all
            sel <- input$rnaDeAnalysisAnnotationTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaDeAnalysisAnnotationTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaDeAnalysisFlagsSelection <- reactive({
        observeEvent(input$clearRnaDeFlagsSelection,{
            proxy <- dataTableProxy("rnaDeAnalysisFlagsTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaDeFlagsSelection,{
            N <- input$rnaDeAnalysisFlagsTable_rows_all
            sel <- input$rnaDeAnalysisFlagsTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaDeAnalysisFlagsTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaDeAnalysisAllSelection <- reactive({
        observeEvent(input$clearRnaDeAllSelection,{
            proxy <- dataTableProxy("rnaDeAnalysisAllTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaDeAllSelection,{
            N <- input$rnaDeAnalysisAllTable_rows_all
            sel <- input$rnaDeAnalysisAllTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaDeAnalysisAllTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaDeAnalysisSummaryDownload <- reactive({
        output$exportRnaDeSummarySelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("diffexpr_summary_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaDeAnalysisSummaryTable_rows_selected
                    if (length(sel)>0) {
                        totalTable <- currentRnaDeTable$totalTable
                        res <- totalTable[,c("gene_name","pvalue","fdr")]
                        fcInd <- grep("_vs_",names(totalTable))
                        avgInd <- grep("mean|median",colnames(totalTable),
                            perl=TRUE)
                        devInd <- grep("sd|mad|IQR",colnames(totalTable),
                            perl=TRUE)
                        res <- cbind(res,totalTable[,c(fcInd,avgInd,devInd)]) 
                        write.table(res[sel,],file=con,sep="\t",quote=FALSE,
                            row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaDeSummaryAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("diffexpr_summary_",tt,".txt", sep='')
            },
            content=function(con) {
                totalTable <- currentRnaDeTable$totalTable
                res <- totalTable[,c("gene_name","pvalue","fdr")]
                fcInd <- grep("_vs_",names(totalTable))
                avgInd <- grep("mean|median",colnames(totalTable),
                    perl=TRUE)
                devInd <- grep("sd|mad|IQR",colnames(totalTable),
                    perl=TRUE)
                res <- cbind(res,totalTable[,c(fcInd,avgInd,devInd)])
                write.table(res,file=con,sep="\t",quote=FALSE,row.names=FALSE)
            }
        )
    })
    
    handleRnaDeAnalysisAnnotationDownload <- reactive({
        output$exportRnaDeAnnotationSelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("diffexpr_annotation_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaDeAnalysisAnnotationTable_rows_selected
                    if (length(sel)>0) {
                        totalTable <- currentRnaDeTable$totalTable
                        res <- totalTable[,c("chromosome","start","end",
                            "gene_id","strand","gene_name","biotype")]
                        write.table(res[sel,],file=con,sep="\t",quote=FALSE,
                            row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaDeAnnotationAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("diffexpr_annotation_",tt,".txt", sep='')
            },
            content=function(con) {
                totalTable <- currentRnaDeTable$totalTable
                res <- totalTable[,c("chromosome","start","end","gene_id",
                    "strand","gene_name","biotype")]
                write.table(res,file=con,sep="\t",quote=FALSE,row.names=FALSE)
            }
        )
    })
    
    handleRnaDeAnalysisFlagsDownload <- reactive({
        output$exportRnaDeFlagsSelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("diffexpr_flags_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaDeAnalysisFlagsTable_rows_selected
                    if (length(sel)>0) {
                        totalTable <- currentRnaDeTable$totalTable
                        res <- totalTable[,c("gene_name","LN","MD","MN","QN",
                            "KN","BT")] 
                        write.table(res[sel,],file=con,sep="\t",quote=FALSE,
                            row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaDeFlagsAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("diffexpr_flags_",tt,".txt", sep='')
            },
            content=function(con) {
                totalTable <- currentRnaDeTable$totalTable
                res <- totalTable[,c("gene_name","LN","MD","MN","QN","KN","BT")]
                write.table(res,file=con,sep="\t",quote=FALSE,row.names=FALSE)
            }
        )
    })
    
    handleRnaDeAnalysisAllDownload <- reactive({
        output$exportRnaDeAllSelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("diffexpr_summary_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaDeAnalysisAllTable_rows_selected
                    if (length(sel)>0) {
                        totalTable <- currentRnaDeTable$totalTable
                        write.table(totalTable[sel,],file=con,sep="\t",
                            quote=FALSE,row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaDeAllAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("diffexpr_summary_",tt,".txt", sep='')
            },
            content=function(con) {
                totalTable <- currentRnaDeTable$totalTable
                write.table(totalTable,file=con,sep="\t",quote=FALSE,
                    row.names=FALSE)
            }
        )
    })
    
    # Now the sliders
    statSliderUpdate <- reactive({
        if (input$statThresholdType=="pvalue") {
            index <- as.integer(input$pvalue)
            # Dirty hack for the 1st loading...
            if (index==0) index <- 10
            currentRnaDeTable$tableFilters$p <- statScoreValues[index+1]
        }
        else if (input$statThresholdType=="fdr") {
            index <- as.integer(input$fdr)
            if (index==0) index <- 10
            currentRnaDeTable$tableFilters$fdr <- statScoreValues[index+1]
        }
    })
    
    valueScaleUpdate <- reactive({
        if (input$rnaDeValueScaleRadio=="natural")
            currentRnaDeTable$tableFilters$scale <- "natural"
        else if (input$rnaDeValueScaleRadio=="log2")
            currentRnaDeTable$tableFilters$scale <- "log2"
    })
    
    foldChangeSliderUpdate <- reactive({
        if (currentRnaDeTable$tableFilters$scale=="natural") {
            index <- as.integer(input$fcNatural)
            # Dirty hack for the 1st loading...
            if (any(index==0)) index <- c(10,20)
            currentRnaDeTable$tableFilters$fc <- fcNatValues[index+1]
        }
        else if (currentRnaDeTable$tableFilters$scale=="log2")
            currentRnaDeTable$tableFilters$fc <- input$fcLog
    })
    
    filterByGeneUpdate <- reactive({
        if (!isEmpty(input$rnaDeShowSpecificGenes))
            currentRnaDeTable$tableFilters$genes <- input$rnaDeShowSpecificGenes
        else
            currentRnaDeTable$tableFilters$genes <- NULL
    })
    
    filterByChromosomeUpdate <- reactive({
        if (!isEmpty(input$customDeChr))
            currentRnaDeTable$tableFilters$chr <- input$customDeChr
        else
            currentRnaDeTable$tableFilters$chr <- NULL
    })
    
    filterByBiotypeUpdate <- reactive({
        if (input$rnaDeAnalyzedBiotypeFilter) {
            bts <- getBiotypes(currentMetadata$genome)
            names(bts) <- bts
            wh <- names(which(sapply(bts,function(b) {
                if (!isEmpty(input[[paste(b,"asFilter",sep="_")]]))
                    return(input[[paste(b,"asFilter",sep="_")]])
                else
                    return(FALSE)
            })))
            currentRnaDeTable$tableFilters$bt <- wh
        }
        else
            currentRnaDeTable$tableFilters$bt <- NULL
    })
    
    updateMaPlot <- reactive({
        if (isEmpty(currentPipelineOutput$counts) 
            || isEmpty(input$rnaDeCurrentContrast))
            return()
        
        classList <- currentPipelineOutput$classList
        contrastList <- currentPipelineOutput$contrastList
        p <- currentPipelineOutput$pValue[[contrastList[1]]][,1,drop=FALSE]
        fdr <- currentPipelineOutput$fdr[[contrastList[1]]][,1,drop=FALSE]
        expr <- which(!is.na(p))
        ann <- currentPipelineOutput$annotation[expr,]
        tab <- currentPipelineOutput$counts[expr,]
        p <- p[expr,,drop=FALSE]
        fdr <- fdr[expr,,drop=FALSE]
        
        fct <- currentRnaDeTable$tableFilters$fc
        pcut <- currentRnaDeTable$tableFilters$p
        fdrcut <- currentRnaDeTable$tableFilters$fdr
        
        if (input$rnaDeValueScaleRadio=="natural")
            fct <- log2(fct)
        if (input$statThresholdType=="pvalue") {
            sc <- pcut
            SC <- p
        }
        else {
            sc <- fdrcut
            SC <- fdr
        }
        
        fcMat <- round(log2(makeFoldChange(contrastList[1],classList,tab)),3)
        aMat <- round(makeA(contrastList[1],classList,tab),3)
        
        #up <- which(fcMat[,1]>=fct[2] & SC<sc)
        #down <- which(fcMat[,1]<=fct[1] & SC<sc)
        #status <- rep("Neutral",nrow(tab))
        #status[up] <- "Up"
        #status[down] <- "Down"
        
        status <- rep(NA,nrow(tab))
        cat.ind <- list()
        cat.ind$up <- which(fcMat[,input$rnaDeCurrentContrast]>=fct[2] & SC<sc)
        cat.ind$down <- 
            which(fcMat[,input$rnaDeCurrentContrast]<=fct[1] & SC<sc)
        cat.ind$neutral <- setdiff(1:nrow(tab),c(cat.ind$up,cat.ind$down))
        names(cat.ind) <- c(
            paste("Up (",length(cat.ind$up)," genes)",sep=""),
            paste("Down (",length(cat.ind$down)," genes)",sep=""),
            paste("Neutral (",length(cat.ind$neutral)," genes)",sep="")
        )
        status[cat.ind[[1]]] <- names(cat.ind)[1]
        status[cat.ind[[2]]] <- names(cat.ind)[2]
        status[cat.ind[[3]]] <- names(cat.ind)[3]
        color.list <- unlist(maPlots$maColours)
        size.list <- c(1,1,0.5)
        names(color.list) <- names(size.list) <- names(cat.ind)
        #names(color.list) <- names(size.list) <- c("Up","Down","Neutral")

        maplot.data <- data.frame(
            A=aMat[,1],
            M=fcMat[,1],
            Status=status,
            Gene=as.character(ann$gene_name)
        )
        rownames(maplot.data) <- as.character(ann$gene_id)
        maPlots$maData <- maplot.data
        
        if (!is.null(currentRnaDeTable$tableFilters$genes)) {
            g <- currentRnaDeTable$tableFilters$genes
            ts <- maplot.data[g,,drop=FALSE]
            ts$Status <- rep(paste("Selected (",length(g)," genes)",sep=""),
                nrow(ts))
            size.list <- c(size.list,5)
            color.list <- c(color.list,"orange")
            names(size.list)[4] <- names(color.list)[4] <- ts$Status[1]
            maplot.data <- rbind(maplot.data,ts)
        }
        
        #maPlots$maPlot <- maplot.data
        maPlots$maPlot <- ggplot(data=maplot.data) + 
            geom_point(aes(x=A,y=M,colour=Status,fill=Status,size=Status,
                text=Gene)) +
            theme_bw() +
            xlab("\nAverage expression") +
            ylab("Fold change (log2)\n") +
            coord_cartesian(xlim=maPlots$maZoom$x,ylim=maPlots$maZoom$y) +
            theme(
                axis.title.x=element_text(size=12),
                axis.title.y=element_text(size=12),
                axis.text.x=element_text(size=10,face="bold"),
                axis.text.y=element_text(size=10,face="bold"),
                legend.key=element_blank(),
                legend.position="bottom"
            ) +
            scale_color_manual(values=color.list) +
            scale_fill_manual(values=color.list) + 
            scale_size_manual(values=size.list) + 
            guides(size=FALSE)
    })
    
    updateMaPlotColours <- reactive({
        macs <- isolate(maPlots$maColours)
        if (is.null(names(macs)))
            names(macs) <- c("Up","Down","Neutral")
        lapply(names(macs),function(x) {
            observeEvent(input[[paste("maPlotColour_",x,sep="")]],{
                newc <- input[[paste("maPlotColour_",x,sep="")]]
                if (newc!=macs[[x]])
                    macs[[x]] <- newc
                    #maPlots$maColours[[x]] <- newc
            })
        })
        maPlots$maColours <- macs
    })
    
    return(list(
        rnaDeTotalTable=rnaDeTotalTable,
        handleRnaDeAnalysisSummarySelection=handleRnaDeAnalysisSummarySelection,
        handleRnaDeAnalysisAnnotationSelection=
            handleRnaDeAnalysisAnnotationSelection,
        handleRnaDeAnalysisFlagsSelection=handleRnaDeAnalysisFlagsSelection,
        handleRnaDeAnalysisAllSelection=handleRnaDeAnalysisAllSelection,
        handleRnaDeAnalysisSummaryDownload=handleRnaDeAnalysisSummaryDownload,
        handleRnaDeAnalysisAnnotationDownload=
            handleRnaDeAnalysisAnnotationDownload,
        handleRnaDeAnalysisFlagsDownload=handleRnaDeAnalysisFlagsDownload,
        handleRnaDeAnalysisAllDownload=handleRnaDeAnalysisAllDownload,
        statSliderUpdate=statSliderUpdate,
        valueScaleUpdate=valueScaleUpdate,
        foldChangeSliderUpdate=foldChangeSliderUpdate,
        filterByGeneUpdate=filterByGeneUpdate,
        filterByChromosomeUpdate=filterByChromosomeUpdate,
        filterByBiotypeUpdate=filterByBiotypeUpdate,
        updateMaPlot=updateMaPlot,
        updateMaPlotColours=updateMaPlotColours
    ))
}

diffExprTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    maPlots <- allReactiveVars$maPlots
     
    output$rnaDeAnalysisSummary <- renderUI({
        if (is.null(currentRnaDeTable$totalTable))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please run a differential expression analysis ",
                    "first."))
            )
        else {
            totalTable <- currentRnaDeTable$totalTable
            res <- totalTable[,c("gene_name","pvalue","fdr")]
            fcInd <- grep("_vs_",names(totalTable))
            avgInd <- grep("mean|median",colnames(totalTable),perl=TRUE)
            devInd <- grep("sd|mad|IQR",colnames(totalTable),perl=TRUE)
            res <- cbind(res,totalTable[,c(fcInd,avgInd,devInd)])
            output$rnaDeAnalysisSummaryTable <- 
                DT::renderDataTable(
                    res,
                    class="display compact",
                    rownames=FALSE,
                    options=list(
                        searchHighlight=TRUE,
                        pageLength=10,
                        lengthMenu=c(10,20,50,100)
                    )
                )
            list(
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaDeAnalysisSummaryTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaDeSummarySelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaDeSummarySelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeSummarySelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeSummaryAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaDeAnalysisAnnotation <- renderUI({
        if (is.null(currentRnaDeTable$totalTable))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please run a differential expression analysis ",
                    "first."))
            )
        else {
            totalTable <- currentRnaDeTable$totalTable
            res <- totalTable[,c("chromosome","start","end","gene_id",
                "strand","gene_name","biotype")]
            output$rnaDeAnalysisAnnotationTable <- 
                DT::renderDataTable(
                    res,
                    class="display compact",
                    rownames=FALSE,
                    options=list(
                        searchHighlight=TRUE,
                        pageLength=10,
                        lengthMenu=c(10,20,50,100)
                    )
                )
            list(
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaDeAnalysisAnnotationTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaDeAnnotationSelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaDeAnnotationSelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeAnnotationSelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeAnnotationAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaDeAnalysisFlags <- renderUI({
        if (is.null(currentRnaDeTable$totalTable))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please run a differential expression analysis ",
                    "first."))
            )
        else {
            totalTable <- currentRnaDeTable$totalTable
            res <- totalTable[,c("gene_name","LN","MD","MN","QN","KN","BT")]
            output$rnaDeAnalysisFlagsTable <- 
                DT::renderDataTable(
                    res,
                    class="display compact",
                    rownames=FALSE,
                    options=list(
                        searchHighlight=TRUE,
                        pageLength=10,
                        lengthMenu=c(10,20,50,100)
                    )
                )
            list(
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaDeAnalysisFlagsTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaDeFlagsSelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaDeFlagsSelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeFlagsSelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeFlagsAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaDeAnalysisAll <- renderUI({
        if (is.null(currentRnaDeTable$totalTable))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please run a differential expression analysis ",
                    "first."))
            )
        else {
            totalTable <- currentRnaDeTable$totalTable
            ann <- totalTable[,c("chromosome","start","end","gene_id",
                "strand","gene_name","biotype")]
            sta <- totalTable[,c("pvalue","fdr")]
            fla <- totalTable[,c("gene_name","LN","MD","MN","QN","KN","BT")]
            fcInd <- grep("_vs_",names(totalTable))
            avgInd <- grep("mean|median",colnames(totalTable),perl=TRUE)
            devInd <- grep("sd|mad|IQR",colnames(totalTable),perl=TRUE)
            res <- cbind(ann,sta,totalTable[,c(fcInd,avgInd,devInd)],fla)
            output$rnaDeAnalysisAllTable <-
                DT::renderDataTable(
                    res,
                    class="display compact",
                    rownames=FALSE,
                    options=list(
                        searchHighlight=TRUE,
                        pageLength=10,
                        lengthMenu=c(10,20,50,100)
                    )
                )
            list(
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaDeAnalysisAllTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaDeAllSelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaDeAllSelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeAllSelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeAllAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaDePipelineControl <- renderUI({
        if (is.null(currentMetadata$final))
            list(
                div(style="font-weight:bold","Select control condition"),
                helpText("Please create a dataset first from the 'Data ",
                    "selector' menu on the top")
            )
        else {
            cls <- unique(as.character(currentMetadata$final$class))
            selectInput(
                inputId="rnaDePipelineControl",
                label="Select control condition",
                choices=cls
            )
        }
    })
    
    output$checkboxBiotypeListRna <- renderUI({
        bts <- getBiotypes(currentMetadata$genome)
        lapply(bts,function(b) {
            checkboxInput(
                inputId=b,
                label=b,
                value=FALSE
            )
        })
    })
    
    output$setDeChrs <- renderUI({
        if (is.null(currentMetadata$final))
            disabled(selectizeInput(
                inputId="customDeChr",
                label="Show selected chromosomes",
                multiple=TRUE,
                choices=NULL
            ))
        else
            disabled(selectizeInput(
                inputId="customDeChr",
                label="Show selected chromosomes",
                multiple=TRUE,
                choices=getValidChromosomes("hg19")
            ))
    })
    
    output$checkboxBiotypeListAnalyzedRna <- renderUI({
        bts <- getBiotypes(currentMetadata$genome)
        lapply(bts,function(b) {
            checkboxInput(
                inputId=paste(b,"asFilter",sep="_"),
                label=b,
                value=FALSE
            )
        })
    })
    
    output$rnaDeMAPlot <- renderPlot({
        maPlots$maPlot
    })
    
    output$maPlotColours <- renderUI({
        c <- maPlots$maColours
        if (is.null(names(c)))
            names(c) <- c("Up","Down","Neutral")
        lapply(names(c),function(x,c) {
            colourInput(
                inputId=paste("maPlotColour_",x,sep=""),
                label=paste(x,"colour"),
                value=c[[x]]
            )
        },c)
    })
    
    output$exportRnaDeMAPlotPDF <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("ma_plot_",tt,".pdf", sep='')
        },
        content=function(con) {
            ggsave(filename=con,plot=maPlots$maPlot,
                width=10,height=7)
        }
    )
    
    output$exportRnaDeMAPlotPNG <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("ma_plot_",tt,".png", sep='')
        },
        content=function(con) {
            ggsave(filename=con,plot=maPlots$maPlot,
                width=10,height=7)
        }
    )
    
    output$exportRnaDeMAPlotGG2 <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("ma_plot_",tt,".rda", sep='')
        },
        content=function(con) {
            gg <- maPlots$maPlot
            save(gg,file=con)
        }
    )
    
    #output$rnaDeMAPlot <- renderPlotly({
    #   if (nrow(maPlots$maPlot)==1 && maPlots$maPlot$Status==1) {
    #        ax <- list(
    #          title="",zeroline=FALSE,showline=FALSE,
    #          showticklabels=FALSE,showgrid=FALSE
    #        )
    #        te <- list(size=38,color=toRGB("grey20"))
    #        plot_ly(
    #            data=maPlots$maPlot,x=A,y=M,
    #            text=Gene,mode="text",textfont=te
    #        ) %>%
    #        layout(xaxis=ax,yaxis=ax)
    #    }
    #    else {
    #       if (!is.null(currentRnaDeTable$tableFilters$genes)) {
    #            g <- currentRnaDeTable$tableFilters$genes
    #            cols <- c("blue","grey80")
    #            dat <- maPlots$maPlot[g,]
    #           p <- plot_ly(
    #               data=dat,x=A,y=M,
    #               type="scatter",mode="markers",
    #               text=Gene,color=Status,colors=cols
    #           )
    #        }
    #        else {
    #           cols <- c("green3","grey50","red3")
    #           #size.list <- c(1,1,0.5)
    #           dat <- maPlots$maPlot
    #           dat <- dat[order(dat$Status),]
    #           p <- plot_ly(
    #               data=dat,x=A,y=M,
    #               type="scatter",mode="markers",
    #               text=Gene,color=Status,colors=cols
    #           )
    #       }
    #       
    #        p <- layout(p,
    #            xaxis=list(title="Average expression"),
    #            yaxis=list(title="Fold change (log2)"),
    #            legend=list(x=1,y=1)
    #        )
    #    }
    #})
}

diffExprTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    customRnaRegions <- allReactiveVars$customRnaRegions
    currentPipelineOutput <- allReactiveVars$currentPipelineOutput
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    
    diffExprTabPanelReactiveEvents <- 
        diffExprTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    runPipeline <- diffExprTabPanelReactiveEvents$runPipeline
    updateMaZoomCoords <- diffExprTabPanelReactiveEvents$updateMaZoomCoords
    resetMaZoom <- diffExprTabPanelReactiveEvents$resetMaZoom
    
    diffExprTabPanelReactiveExprs <- 
        diffExprTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    rnaDeTotalTable <- diffExprTabPanelReactiveExprs$rnaDeTotalTable
    handleRnaDeAnalysisSummarySelection <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisSummarySelection
    handleRnaDeAnalysisAnnotationSelection <- 
            diffExprTabPanelReactiveExprs$handleRnaDeAnalysisAnnotationSelection
    handleRnaDeAnalysisFlagsSelection <-
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisFlagsSelection
    handleRnaDeAnalysisAllSelection <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisAllSelection
    handleRnaDeAnalysisSummaryDownload <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisSummaryDownload
    handleRnaDeAnalysisAnnotationDownload <-
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisAnnotationDownload
    handleRnaDeAnalysisFlagsDownload <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisFlagsDownload
    handleRnaDeAnalysisAllDownload <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisAllDownload
    statSliderUpdate <- diffExprTabPanelReactiveExprs$statSliderUpdate
    valueScaleUpdate <- diffExprTabPanelReactiveExprs$valueScaleUpdate
    foldChangeSliderUpdate <- 
        diffExprTabPanelReactiveExprs$foldChangeSliderUpdate
    filterByGeneUpdate <- diffExprTabPanelReactiveExprs$filterByGeneUpdate
    filterByChromosomeUpdate <- 
        diffExprTabPanelReactiveExprs$filterByChromosomeUpdate
    filterByBiotypeUpdate <- 
        diffExprTabPanelReactiveExprs$filterByBiotypeUpdate
    updateMaPlot <- diffExprTabPanelReactiveExprs$updateMaPlot
    updateMaPlotColours <- 
        diffExprTabPanelReactiveExprs$updateMaPlotColours
    
    diffExprTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
    
    observe({
        rnaDeTotalTable()
    })
    
    observe({
        if (!isEmpty(currentPipelineOutput$contrastList)) {
            conds <- strsplit(currentPipelineOutput$contrastList[1],"_vs_")[[1]]
            cn <- paste(conds[1:(length(conds)-1)],"_vs_",
                conds[length(conds)],sep="")
            names(cn) <- paste(conds[1:(length(conds)-1)],"vs",
                conds[length(conds)])
            updateSelectizeInput(session,"rnaDeCurrentContrast",
                choices=cn,
                server=TRUE,
                selected=cn[1]
            )
        }
    })
    
    observe({
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"rnaDeKnownFilter",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$rnaDeKnownFilter})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"rnaDeKnownFilter",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"rnaDeShowSpecificGenes",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            if (input$includeCustomRegions) {
                ts <- as.character(customRnaRegions$name)
                names(ts) <- ts
                geneNames <- c(geneNames,ts)
            }
            g <- isolate({input$rnaDeShowSpecificGenes})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"rnaDeShowSpecificGenes",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        if (isEmpty(currentMetadata$final))
            shinyjs::disable("performDeAnalysis")
        else
            shinyjs::enable("performDeAnalysis")
    })
    
    observe({
        if (input$rnaDeGeneFilter=="quantile") {
            qq <- as.numeric(input$rnaDeQuantileFilter)
            if (is.na(qq) || qq<0 || qq>1) {
                output$rnaDeSettingsError <- renderUI({
                    div(class="error-message",paste("The quantile ",
                        "must be a number between 0 and 1!",sep=""))
                })
                shinyjs::disable("performDeAnalysis")
            }
            else {
                output$rnaDeSettingsError <- renderUI({div()})
                if (isEmpty(currentMetadata$final))
                    shinyjs::disable("performDeAnalysis")
                else
                    shinyjs::enable("performDeAnalysis")
            }
        }
        else {
            if (isEmpty(currentMetadata$final))
                shinyjs::disable("performDeAnalysis")
            else
                shinyjs::enable("performDeAnalysis")
        }           
    })
    
    observe({
        if (input$toggleRnaDeZoom)
            shinyjs::enable("resetRnaDeZoom")
        else
            shinyjs::disable("resetRnaDeZoom")
    })
    
    observe({
        if(is.null(currentRnaDeTable$totalTable)) {
            shinyjs::disable("statThresholdType")
            shinyjs::disable("pvalue")
            shinyjs::disable("fdr")
            shinyjs::disable("fcNatural")
            shinyjs::disable("fcLog")
            shinyjs::disable("rnaDeCurrentContrast")
            shinyjs::disable("rnaDeShowSpecificGenes")
            shinyjs::disable("rnaDeAnalyzedBiotypeFilter")
            shinyjs::disable("rnaDeValueCompRadio")
            shinyjs::disable("rnaDeValueScaleRadio")
            shinyjs::disable("rnaDeValueAverageRadio")
            shinyjs::disable("rnaDeValueDeviationRadio")
            shinyjs::disable("customDeChr")
        }
        else {
            shinyjs::enable("statThresholdType")
            shinyjs::enable("pvalue")
            shinyjs::enable("fdr")
            shinyjs::enable("fcNatural")
            shinyjs::enable("fcLog")
            shinyjs::enable("rnaDeCurrentContrast")
            shinyjs::enable("rnaDeShowSpecificGenes")
            shinyjs::enable("rnaDeAnalyzedBiotypeFilter")
            shinyjs::enable("rnaDeValueCompRadio")
            shinyjs::enable("rnaDeValueScaleRadio")
            shinyjs::enable("rnaDeValueAverageRadio")
            shinyjs::enable("rnaDeValueDeviationRadio")
            shinyjs::enable("customDeChr")
        }
    })
    
    observe({
        handleRnaDeAnalysisSummarySelection()
        handleRnaDeAnalysisAnnotationSelection()
        handleRnaDeAnalysisFlagsSelection()
        handleRnaDeAnalysisAllSelection()
        handleRnaDeAnalysisSummaryDownload()
        handleRnaDeAnalysisAnnotationDownload()
        handleRnaDeAnalysisFlagsDownload()
        handleRnaDeAnalysisAllDownload()
    })
    
    observe({
        statSliderUpdate()
        valueScaleUpdate()
        foldChangeSliderUpdate()
        filterByGeneUpdate()
        filterByChromosomeUpdate()
        filterByBiotypeUpdate()
    })
    
    observe({
        updateMaPlot()
        updateMaPlotColours()
        updateMaZoomCoords()
    })
    
    observe({
        resetMaZoom()
    })
    
    observe({
        tryCatch({
            shinyjs::disable("performDeAnalysis")
            # Clear interaction with MA plot
            updateCheckboxInput(session,"toggleRnaDeTableUpdate",value=FALSE)
            updateCheckboxInput(session,"toggleRnaDeZoom",value=FALSE)
            # Clear filters
            currentRnaDeTable$tableFilters <- list(
                p=0.05,
                fdr=0.05,
                scale="natural",
                fc=c(0.5,2),
                bt=NULL,
                genes=NULL,
                chr=NULL
            )
            runPipeline()
        },error=function(e) {
            #print(e)
        },
        finally={
            if (isEmpty(currentMetadata$final))
                shinyjs::disable("performDeAnalysis")
            else
                shinyjs::enable("performDeAnalysis")
        })
    })
}

